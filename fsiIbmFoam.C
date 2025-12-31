/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    fsiIbmFoam

Description
    Transient solver for incompressible, turbulent flow with immersed boundary
    method (IBM) and fluid-structure interaction (FSI). Based on PISO algorithm
    with turbulence model support (k-epsilon, k-omega, LES, etc.).

    Features:
    - Turbulent flow simulation using OpenFOAM turbulence models
    - Immersed Boundary Method for complex geometries
    - Fluid-Structure Interaction coupling
    - Parallel processing with enhanced pressure interpolation
    - Real-time IBM boundary updates following solid deformation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "patchToPatchInterpolation.H"
#include "tractionDisplacement/tractionDisplacementFvPatchVectorField.H"
#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValueTetPolyPatchFields.H"
#include "pisoControl.H"

#include "pointMesh.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

#include "interpolationCellPoint.H"
#include "interpolationCellPointFace.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "immersedBoundaryFvPatchField.H"
#include "mixedIbFvPatchField.H"
#include "processorFvPatchField.H"

#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    //- NOTE Create solid mesh.
    // tesbo, 07/May/2025
    #include "createStressMesh.H"


    //- NOTE Create normal mesh - fluid background - including IBM boundary
    // tesbo, 07/May/2025
    #include "createMesh.H"

    //- NOTE Create interpolation mesh for pressure transfer
    // tesbo, 07/May/2025
    #include "createMeshItpPressure.H"

    pisoControl piso(mesh);

    //- NOTE Create fields with turbulence model support
    // tesbo, 07/May/2025
    #include "createFields.H"

    #include "createStressFields.H"
    #include "readMechanicalProperties.H"
    #include "readCouplingProperties.H"
    #include "createTimeControls.H"
    
    #include "initContinuityErrs.H"

    Info<< "Turbulence model type: " << turbulence->type() << endl;

    bool isLaminar = (turbulence->type() == "laminar");

    if (isLaminar)
    {
        Info<< "Using laminar flow model" << endl;
    }
    else
    {
        Info<< "Using turbulent flow model: " << turbulence->type() << endl;
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // **CRITICAL: Initialize V0 for restart calculations**
    // When restarting from a time > 0, V0 (old cell volumes) may not be available
    // This causes fatal errors in mesh motion and IB geometry calculations
    // Force initialization to current volumes for restart safety

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;


        //- NOTE Calculate force with global pressure interpolation
        // Enhanced with wall shear stress capability (moved after variable definitions)
        // tesbo, 07/May/2025
        #include "setCellToPatchForcesSimple.H"
        
        #include "solveSolid.H"
        
        #include "setMotion.H"

        //- NOTE Update IBM with global solid mesh
        // Critical for turbulent boundary layer interaction with IBM
        // tesbo, 08/May/2025
        #include "updateIbMeshGlobal.H"

        //- NOTE Update ibMask for post-processing
        // tesbo, 11/Dec/2025
        #include "updateIbMask.H"

        //- Force IBM boundary real-time update after solid mesh motion
        // This ensures IBM boundaries follow solid deformation in real-time
        #include "forceIbmUpdate.H"

        //- NOTE Solve fluid with IBM
        // tesbo, 08/May/2025
        #include "solveFluid.H"

        turbulence->correct();

        if (Pstream::parRun() && !isLaminar)
        {
            const volScalarField& k = turbulence->k();
            const volScalarField& epsilon = turbulence->epsilon();
            const volScalarField& nut = turbulence->nut();
            
            const_cast<volScalarField&>(k).boundaryField().updateCoupledPatchFields();
            const_cast<volScalarField&>(epsilon).boundaryField().updateCoupledPatchFields();
            const_cast<volScalarField&>(nut).boundaryField().updateCoupledPatchFields();
        }

        if (runTime.outputTime())
        {
            runTime.write();
        }

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
