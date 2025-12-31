{
    pointVectorField solidPointsDispl =
        cpi.interpolate(Usolid - Usolid.oldTime());

    vectorField newPoints =
        stressMesh.points()
      + solidPointsDispl.internalField();

    stressMesh.movePoints(newPoints);


    #include "volContinuity.H"

    Info << "Motion magnitude: mean = "
        << average(mag(Usolid.boundaryField()[solidPatchID]))
        << " max = "
        << max(mag(Usolid.boundaryField()[solidPatchID])) << endl;
}
