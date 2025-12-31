ibMask.internalField() = 0.0;

forAll(mesh.boundaryMesh(), patchI)
{
    if (isA<immersedBoundaryPolyPatch>(mesh.boundaryMesh()[patchI]))
    {
        const immersedBoundaryPolyPatch& ibPatch = refCast
        <
            const immersedBoundaryPolyPatch
        >
        (
            mesh.boundaryMesh()[patchI]
        );
        
        const labelList& ibCells = ibPatch.ibCells();
        forAll(ibCells, cellI)
        {
            if (ibCells[cellI] < mesh.nCells())
            {
                ibMask[ibCells[cellI]] = 1.0;
            }
        }
        
        const labelList& deadCells = ibPatch.deadCells();
        forAll(deadCells, cellI)
        {
            if (deadCells[cellI] < mesh.nCells())
            {
                ibMask[deadCells[cellI]] = 1.0;
            }
        }
    }
}

if (Pstream::parRun())
{
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<processorFvPatch>(mesh.boundary()[patchI]))
        {
            const labelList& faceCells = mesh.boundary()[patchI].faceCells();
            forAll(faceCells, faceI)
            {
                label cellI = faceCells[faceI];
                if (cellI < mesh.nCells())
                {
                    if (ibMask[cellI] == 1.0)
                    {
                        ibMask[cellI] = 3.0;
                    }
                    else if (ibMask[cellI] == 0.0)
                    {
                        ibMask[cellI] = 2.0;
                    }
                }
            }
        }
    }
}

ibMask.correctBoundaryConditions();

if (runTime.outputTime())
{
    label ibmCells = sum(pos(ibMask.internalField() - 0.5));
    label procCells = sum(pos(ibMask.internalField() - 1.5)); 
    label intersectionCells = sum(pos(ibMask.internalField() - 2.5));
    
    Info<< "ibMask updated: IBM=" << ibmCells 
        << ", Proc=" << procCells 
        << ", Intersection=" << intersectionCells << endl;
} 
 
