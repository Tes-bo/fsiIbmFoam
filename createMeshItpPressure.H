Foam::Info
    << "Create mesh for p interpolator "
    << Foam::nl << Foam::endl;

Foam::fvMesh meshItpPressure
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::NO_READ,
        Foam::IOobject::AUTO_WRITE
    )
);

meshItpPressure.removeFvBoundary();

const polyBoundaryMesh& patches = mesh.boundaryMesh();
List<polyPatch*> newPatches(patches.size());
label nNewPatches = 0;

forAll(patches, patchi)
{
    if (!isA<immersedBoundaryPolyPatch>(patches[patchi]))
    {
        newPatches[nNewPatches] = patches[patchi].clone
        (
            meshItpPressure.boundaryMesh(),
            nNewPatches,
            patches[patchi].size(),
            patches[patchi].start()
        ).ptr();
        nNewPatches++;
    }
}

newPatches.setSize(nNewPatches);

meshItpPressure.addFvPatches(newPatches);


