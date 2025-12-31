{
    volScalarField pCalc
    (
        IOobject
        (
            "pCalc_internal",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        meshItpPressure,
        dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0)
    );
    pCalc.internalField() = p.internalField();

    volVectorField wallShearStress
    (
        IOobject
        (
            "wallShearStress",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0, 2, -2, 0, 0, 0, 0), vector::zero)
    );

    volTensorField gradU = fvc::grad(U);
    volSymmTensorField devSymmGradU = dev(symm(gradU));
    volSymmTensorField Reff = turbulence->nuEff() * devSymmGradU * 2.0;
    
    wallShearStress.internalField() = vector::zero;
    
    
    word solidPatchName = "consoleSolid";  // Default solid patch name
    label solidPatchID = stressMesh.boundaryMesh().findPatchID(solidPatchName);
    
    if (solidPatchID < 0)
    {
        FatalErrorIn("setCellToPatchForcesSimple.H")
            << "Cannot find solid patch: " << solidPatchName
            << exit(FatalError);
    }

    const vectorField& solidPatchCenters =
        stressMesh.boundaryMesh()[solidPatchID].faceCentres();

    scalarField solidPatchPressure(solidPatchCenters.size(), 0.0);

    if (Pstream::parRun())
    {
        List<pointField> allCenters(Pstream::nProcs());
        List<scalarField> allPressures(Pstream::nProcs());
        
        allCenters[Pstream::myProcNo()] = mesh.cellCentres();
        allPressures[Pstream::myProcNo()] = pCalc.internalField();
        
        Pstream::gatherList(allCenters);
        Pstream::scatterList(allCenters);
        Pstream::gatherList(allPressures);
        Pstream::scatterList(allPressures);
        
        forAll(solidPatchCenters, faceI)
        {
            const point& targetPoint = solidPatchCenters[faceI];
            scalar minDistance = GREAT;
            scalar nearestPressure = 0.0;
            
            forAll(allCenters, procI)
            {
                const pointField& procCenters = allCenters[procI];
                const scalarField& procPressures = allPressures[procI];
                
                forAll(procCenters, cellI)
                {
                    scalar distance = mag(procCenters[cellI] - targetPoint);
                    if (distance < minDistance)
                    {
                        minDistance = distance;
                        nearestPressure = procPressures[cellI];
                    }
                }
            }
            
            solidPatchPressure[faceI] = nearestPressure;
        }
    }
    else
    {
        forAll(solidPatchCenters, faceI)
        {
            label cellI = mesh.findCell(solidPatchCenters[faceI]);
            if (cellI != -1)
            {
                solidPatchPressure[faceI] = pCalc[cellI];
            }
            else
            {
                label nearestCell = mesh.findNearestCell(solidPatchCenters[faceI]);
                if (nearestCell != -1)
                {
                    solidPatchPressure[faceI] = pCalc[nearestCell];
                }
            }
        }
    }

    vectorField solidPatchViscousForce(solidPatchCenters.size(), vector::zero);
    
    const vectorField& solidPatchNormals = 
        stressMesh.boundaryMesh()[solidPatchID].faceNormals();
    
    if (Pstream::parRun())
    {
        List<pointField> allCenters(Pstream::nProcs());
        List<symmTensorField> allStress(Pstream::nProcs());
        
        allCenters[Pstream::myProcNo()] = mesh.cellCentres();
        allStress[Pstream::myProcNo()] = Reff.internalField();
        
        Pstream::gatherList(allCenters);
        Pstream::scatterList(allCenters);
        Pstream::gatherList(allStress);
        Pstream::scatterList(allStress);
        
        forAll(solidPatchCenters, faceI)
        {
            const point& targetPoint = solidPatchCenters[faceI];
            const vector& normal = solidPatchNormals[faceI];
            
            scalar minDistance = GREAT;
            symmTensor nearestStress = symmTensor::zero;
            
            forAll(allCenters, procI)
            {
                const pointField& procCenters = allCenters[procI];
                const symmTensorField& procStress = allStress[procI];
                
                forAll(procCenters, cellI)
                {
                    scalar distance = mag(procCenters[cellI] - targetPoint);
                    if (distance < minDistance)
                    {
                        minDistance = distance;
                        nearestStress = procStress[cellI];
                    }
                }
            }
            
            solidPatchViscousForce[faceI] = normal & nearestStress;
        }
    }
    else
    {
        forAll(solidPatchCenters, faceI)
        {
            const vector& normal = solidPatchNormals[faceI];
            
            label cellI = mesh.findCell(solidPatchCenters[faceI]);
            if (cellI != -1)
            {
                solidPatchViscousForce[faceI] = normal & Reff[cellI];
            }
            else
            {
                label nearestCell = mesh.findNearestCell(solidPatchCenters[faceI]);
                if (nearestCell != -1)
                {
                    solidPatchViscousForce[faceI] = normal & Reff[nearestCell];
                }
            }
        }
    }

    tractionDisplacementFvPatchVectorField& tForce =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            Usolid.boundaryField()[solidPatchID]
        );

    vectorField totalTraction(solidPatchCenters.size());
    
    forAll(solidPatchCenters, faceI)
    {
        const vector& normal = solidPatchNormals[faceI];
        
        vector pressureTraction = solidPatchPressure[faceI] * rhoFluid.value() * (-normal);
        
        vector viscousTraction = solidPatchViscousForce[faceI] * rhoFluid.value() * (-1.0);
        
        totalTraction[faceI] = pressureTraction + viscousTraction;
    }
    
    tForce.traction() = totalTraction;
    
    tForce.pressure() = scalarField(solidPatchCenters.size(), 0.0);

    if (runTime.outputTime())
    {
        scalar maxPressure = max(mag(solidPatchPressure));
        scalar avgPressure = average(mag(solidPatchPressure));
        scalar maxViscousForce = max(mag(solidPatchViscousForce));
        scalar avgViscousForce = average(mag(solidPatchViscousForce));
        scalar maxTotalTraction = max(mag(totalTraction));
        scalar avgTotalTraction = average(mag(totalTraction));
        
        scalar viscousRatio = avgViscousForce / (avgPressure + SMALL);
        
        Info<< "FSI Complete Force Statistics:" << nl
            << "    Pressure Forces:" << nl
            << "        Max pressure: " << maxPressure << " Pa" << nl
            << "        Average pressure: " << avgPressure << " Pa" << nl
            << "    Viscous Forces:" << nl
            << "        Max viscous: " << maxViscousForce << " Pa" << nl
            << "        Average viscous: " << avgViscousForce << " Pa" << nl
            << "    Total Traction:" << nl
            << "        Max total: " << maxTotalTraction << " Pa" << nl
            << "        Average total: " << avgTotalTraction << " Pa" << nl
            << "    Force Analysis:" << nl
            << "        Viscous/Pressure ratio: " << viscousRatio << nl
            << "        Solid patch faces: " << solidPatchCenters.size() << nl
            << endl;
    }
} 
