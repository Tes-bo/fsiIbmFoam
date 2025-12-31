    Info << "Reading coupling properties" << endl;
    IOdictionary couplingProperties
    (
        IOobject
        (
            "couplingProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Read solid patch data
    word solidPatchName(couplingProperties.lookup("solidPatch"));

    label solidPatchID =
        stressMesh.boundaryMesh().findPatchID(solidPatchName);



    // Grab solid patch field
    tractionDisplacementFvPatchVectorField& tForce =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            Usolid.boundaryField()[solidPatchID]
        );


