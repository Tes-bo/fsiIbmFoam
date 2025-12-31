if (runTime.outputTime())
{
    Info<< "Updating IBM boundaries at time " << runTime.timeName() << endl;
}

forAll(mesh.boundaryMesh(), patchI)
{
    if (isA<immersedBoundaryPolyPatch>(mesh.boundaryMesh()[patchI]))
    {
        const immersedBoundaryPolyPatch& ibPatch = refCast
        <
            const immersedBoundaryPolyPatch
        >(mesh.boundaryMesh()[patchI]);
        
        const_cast<immersedBoundaryPolyPatch&>(ibPatch).clearOut();
        
        const_cast<immersedBoundaryPolyPatch&>(ibPatch).updateIbMesh();
        
        if (runTime.outputTime())
        {
            Info<< "  IBM patch " << ibPatch.name() 
                << ": ibFaces=" << ibPatch.ibFaces().size()
                << ", ibCells=" << ibPatch.ibCells().size() << endl;
        }
        
        forAll(U.boundaryField(), boundaryPatchI)
        {
            if (boundaryPatchI == patchI)
            {
                U.boundaryField()[boundaryPatchI].updateCoeffs();
                p.boundaryField()[boundaryPatchI].updateCoeffs();
            }
        }
    }
}

mesh.clearOut();
mesh.lduAddr();

if (runTime.outputTime())
{
    Info<< "IBM real-time update completed." << endl;
} 
 
