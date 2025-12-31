{
    if (Pstream::parRun())
    {
        List<pointField> allSolidPoints(Pstream::nProcs());
        List<faceList> allSolidFaces(Pstream::nProcs());
        List<labelList> allSolidOwner(Pstream::nProcs());
        List<labelList> allSolidNeighbour(Pstream::nProcs());
        List<wordList> allPatchNames(Pstream::nProcs());
        List<labelList> allPatchSizes(Pstream::nProcs());
        List<labelList> allPatchStarts(Pstream::nProcs());
        
        allSolidPoints[Pstream::myProcNo()] = stressMesh.points();
        allSolidFaces[Pstream::myProcNo()] = stressMesh.faces();
        allSolidOwner[Pstream::myProcNo()] = stressMesh.faceOwner();
        allSolidNeighbour[Pstream::myProcNo()] = stressMesh.faceNeighbour();
        
        if (runTime.outputTime())
        {
            Pout<< "Processor " << Pstream::myProcNo() << " solid mesh: " 
                << stressMesh.points().size() << " points, " 
                << stressMesh.faces().size() << " faces" << endl;
        }
        
        wordList localPatchNames(stressMesh.boundaryMesh().size());
        labelList localPatchSizes(stressMesh.boundaryMesh().size());
        labelList localPatchStarts(stressMesh.boundaryMesh().size());
        
        label localConsoleSolidFaces = 0;
        
        forAll(stressMesh.boundaryMesh(), patchI)
        {
            localPatchNames[patchI] = stressMesh.boundaryMesh()[patchI].name();
            localPatchSizes[patchI] = stressMesh.boundaryMesh()[patchI].size();
            localPatchStarts[patchI] = stressMesh.boundaryMesh()[patchI].start();
            
            if (localPatchNames[patchI] == "consoleSolid")
            {
                localConsoleSolidFaces = localPatchSizes[patchI];
                Pout<< "Processor " << Pstream::myProcNo() << " has " 
                    << localConsoleSolidFaces << " consoleSolid faces" << endl;
            }
        }
        
        allPatchNames[Pstream::myProcNo()] = localPatchNames;
        allPatchSizes[Pstream::myProcNo()] = localPatchSizes;
        allPatchStarts[Pstream::myProcNo()] = localPatchStarts;
        
        Pstream::gatherList(allSolidPoints);
        Pstream::scatterList(allSolidPoints);
        Pstream::gatherList(allSolidFaces);
        Pstream::scatterList(allSolidFaces);
        Pstream::gatherList(allSolidOwner);
        Pstream::scatterList(allSolidOwner);
        Pstream::gatherList(allSolidNeighbour);
        Pstream::scatterList(allSolidNeighbour);
        Pstream::gatherList(allPatchNames);
        Pstream::scatterList(allPatchNames);
        Pstream::gatherList(allPatchSizes);
        Pstream::scatterList(allPatchSizes);
        Pstream::gatherList(allPatchStarts);
        Pstream::scatterList(allPatchStarts);
        
        DynamicList<point> globalPoints;
        DynamicList<face> globalFaces;
        DynamicList<label> globalOwner;
        DynamicList<label> globalNeighbour;
        
        List<labelList> pointMaps(Pstream::nProcs());
        
        forAll(allSolidPoints, procI)
        {
            const pointField& procPoints = allSolidPoints[procI];
            pointMaps[procI].setSize(procPoints.size());
            
            forAll(procPoints, pointI)
            {
                pointMaps[procI][pointI] = globalPoints.size();
                globalPoints.append(procPoints[pointI]);
            }
        }
        
        forAll(allSolidFaces, procI)
        {
            const faceList& procFaces = allSolidFaces[procI];
            const labelList& procOwner = allSolidOwner[procI];
            const labelList& procNeighbour = allSolidNeighbour[procI];
            const wordList& procPatchNames = allPatchNames[procI];
            const labelList& procPatchSizes = allPatchSizes[procI];
            const labelList& procPatchStarts = allPatchStarts[procI];
            
            label cellOffset = 0;
            for (label prevProcI = 0; prevProcI < procI; prevProcI++)
            {
                cellOffset += allSolidOwner[prevProcI].size() > 0 ? 
                    max(allSolidOwner[prevProcI]) + 1 : 0;
            }
            
            label nInternalFaces = procFaces.size();
            if (procPatchStarts.size() > 0)
            {
                nInternalFaces = procPatchStarts[0];
            }
            
            for (label faceI = 0; faceI < nInternalFaces; faceI++)
            {
                face f = procFaces[faceI];
                
                forAll(f, fpI)
                {
                    f[fpI] = pointMaps[procI][f[fpI]];
                }
                
                globalFaces.append(f);
                globalOwner.append(procOwner[faceI] + cellOffset);
                globalNeighbour.append(procNeighbour[faceI] + cellOffset);
            }
            
            forAll(procPatchNames, patchI)
            {
                const word& patchName = procPatchNames[patchI];
                bool isProcessorBoundary = 
                    (patchName == "procBoundary") ||
                    (patchName.size() >= 9 && patchName.substr(0, 9) == "processor");
                
                if (!isProcessorBoundary && patchName == "consoleSolid")
                {
                    label patchStart = procPatchStarts[patchI];
                    label patchSize = procPatchSizes[patchI];
                    
                    Pout<< "Time " << runTime.timeName() << " - Adding consoleSolid faces from processor " << procI 
                        << ": " << patchSize << " faces (start: " << patchStart << ")" << endl;
                    
                    for (label i = 0; i < patchSize; i++)
                    {
                        label faceI = patchStart + i;
                        face f = procFaces[faceI];
                        
                        forAll(f, fpI)
                        {
                            f[fpI] = pointMaps[procI][f[fpI]];
                        }
                        
                        globalFaces.append(f);
                        globalOwner.append(procOwner[faceI] + cellOffset);
                        globalNeighbour.append(-1);
                    }
                }
            }
        }
        
        if (globalPoints.size() > 0 && globalFaces.size() > 0)
        {
            Pout<< "Processor " << Pstream::myProcNo() << " time " << runTime.timeName()
                << " - Combined global solid mesh: " << globalPoints.size() 
                << " points, " << globalFaces.size() << " consoleSolid faces" << endl;
            
            DynamicList<point> triPoints;
            DynamicList<triFace> triTriFaces;
            Map<label> pointMap;
           
            pointField globalPointField(globalPoints.shrink());
            
            boundBox solidBounds(globalPointField);
            Info<< "Global solid bounds: " << solidBounds << endl;
            
            forAll(globalFaces, faceI)
            {
                if (globalNeighbour[faceI] == -1)
                {
                    const face& f = globalFaces[faceI];
                    
                    labelList facePointMap(f.size());
                    forAll(f, fpI)
                    {
                        label globalPointI = f[fpI];
                        
                        if (pointMap.found(globalPointI))
                        {
                            facePointMap[fpI] = pointMap[globalPointI];
                        }
                        else
                        {
                            label triPointI = triPoints.size();
                            triPoints.append(globalPoints[globalPointI]);
                            pointMap.insert(globalPointI, triPointI);
                            facePointMap[fpI] = triPointI;
                        }
                    }
                    
                    if (f.size() == 3)
                    {
                        triTriFaces.append
                        (
                            triFace
                            (
                                facePointMap[0],
                                facePointMap[1],
                                facePointMap[2]
                            )
                        );
                    }
                    else if (f.size() > 3)
                    {
                        for (label i = 1; i < f.size() - 1; i++)
                        {
                            triTriFaces.append
                            (
                                triFace
                                (
                                    facePointMap[0],
                                    facePointMap[i],
                                    facePointMap[i+1]
                                )
                            );
                        }
                    }
                }
            }
            
            Info<< "Generated triSurface from consoleSolid: " << triTriFaces.size() << " triangle faces" << endl;
            
            if (triPoints.size() > 0 && triTriFaces.size() > 0)
            {
                triSurface globalTriSurf
                (
                    List<triFace>(triTriFaces.shrink()),
                    pointField(triPoints.shrink())
                );
                
                Pout<< "Processor " << Pstream::myProcNo() << " time " << runTime.timeName()
                    << " - Final consoleSolid triSurface: " << globalTriSurf.size() 
                    << " triangle faces, " << globalTriSurf.points().size() << " points" << endl;
                    
                boundBox triSurfBounds(globalTriSurf.points());
                Pout<< "Processor " << Pstream::myProcNo() << " time " << runTime.timeName()
                    << " - ConsoleSolid triSurface bounds: " << triSurfBounds << endl;
                    
                if (runTime.value() > 0.001) // Skip first time step
                {
                    scalar maxPointDisp = 0;
                    forAll(globalTriSurf.points(), pointI)
                    {
                        maxPointDisp = max(maxPointDisp, mag(globalTriSurf.points()[pointI]));
                    }
                    Pout<< "Processor " << Pstream::myProcNo() << " time " << runTime.timeName()
                        << " - Max triSurface point magnitude: " << maxPointDisp << endl;
                }
                
                forAll(mesh.boundaryMesh(), patchI)
                {
                    if (isA<immersedBoundaryPolyPatch>(mesh.boundaryMesh()[patchI]))
                    {
                        const immersedBoundaryPolyPatch& ibPatchconst = refCast
                        <
                            const immersedBoundaryPolyPatch
                        >
                        (
                            mesh.boundaryMesh()[patchI]
                        );
                        
                        Pout<< "Processor " << Pstream::myProcNo() << " time " << runTime.timeName()
                            << " - Before IBM update: ibFaces=" << ibPatchconst.ibMesh().size() 
                            << ", ibPoints=" << ibPatchconst.ibMesh().points().size() << endl;
                        
                        const_cast<immersedBoundaryPolyPatch&>(ibPatchconst).updateIbMesh
                        (
                            globalTriSurf
                        );
                        
                        Pout<< "Processor " << Pstream::myProcNo() << " time " << runTime.timeName()
                            << " - After IBM update: ibFaces=" << ibPatchconst.ibMesh().size() 
                            << ", ibPoints=" << ibPatchconst.ibMesh().points().size() << endl;
                        
                        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh()).updateMesh();
                        
                        mesh.clearOut();
                        mesh.lduAddr();
                        
                        if (Pstream::parRun())
                        {
                            mesh.clearOut();
                            
                            const_cast<polyBoundaryMesh&>(mesh.boundaryMesh()).updateMesh();
                            
                            mesh.lduAddr();
                            
                            Pout<< "Processor " << Pstream::myProcNo() << " time " << runTime.timeName()
                                << " - IBM mesh update with boundary protection completed" << endl;
                                
                            Pout<< "Processor " << Pstream::myProcNo() << " time " << runTime.timeName()
                                << " - Final IBM verification: ibFaces=" << ibPatchconst.ibMesh().size() 
                                << ", ibPoints=" << ibPatchconst.ibMesh().points().size() << endl;
                        }
                        
                        phi.correctBoundaryConditions();
                        U.correctBoundaryConditions();
                        p.correctBoundaryConditions();
                        
                        forAll(U.boundaryField(), patchI)
                        {
                            if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
                            {
                                U.boundaryField()[patchI].updateCoeffs();
                                p.boundaryField()[patchI].updateCoeffs();
                                
                                Pout<< "Processor " << Pstream::myProcNo() << " time " << runTime.timeName()
                                    << " - Forced IBM boundary condition update for patch " << patchI << endl;
                            }
                        }
                        
                        const lduAddressing& lduAddr = mesh.lduAddr();
                        Pout<< "Processor " << Pstream::myProcNo() << " time " << runTime.timeName()
                            << " - Forced mesh reconstruction completed. New lduAddr size: " 
                            << lduAddr.size() << endl;
                        
                        Pout<< "Processor " << Pstream::myProcNo() << " time " << runTime.timeName()
                            << " - Final verification: ibFaces=" << ibPatchconst.ibMesh().size() 
                            << ", ibPoints=" << ibPatchconst.ibMesh().points().size() << endl;
                    }
                }
            }
            else
            {
                Info<< "Warning: No valid triSurface faces found from global solid mesh" << endl;
            }
        }
        else
        {
            Info<< "Warning: No valid global solid mesh data found" << endl;
        }
    }
    else
    {
        // Serial case: use original method
        labelHashSet includedPatches;
        
        forAll(stressMesh.boundaryMesh(), patchI)
        {
            const auto& patch = stressMesh.boundaryMesh()[patchI];
            if (patch.name() == "consoleSolid")
            {
                includedPatches.insert(patchI);
            }
        }
        
        triSurface triSurf = triSurfaceTools::triangulate
        (
            stressMesh.boundaryMesh(),
            includedPatches
        );
        
        forAll(mesh.boundaryMesh(), patchI)
        {
            if (isA<immersedBoundaryPolyPatch>(mesh.boundaryMesh()[patchI]))
            {
                const immersedBoundaryPolyPatch& ibPatchconst = refCast
                <
                    const immersedBoundaryPolyPatch
                >
                (
                    mesh.boundaryMesh()[patchI]
                );
                
                const_cast<immersedBoundaryPolyPatch&>(ibPatchconst).updateIbMesh
                (
                    triSurf
                );
                
                const_cast<polyBoundaryMesh&>(mesh.boundaryMesh()).updateMesh();
                mesh.lduAddr();
                mesh.clearOut();
                
                phi.correctBoundaryConditions();
                U.correctBoundaryConditions();
                p.correctBoundaryConditions();
            }
        }
    }

    U.correctBoundaryConditions();
} 
 
