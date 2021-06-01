/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    RSWEFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for inviscid shallow-water equations with rotation.

    If the geometry is 3D then it is assumed to be one layers of cells and the
    component of the velocity normal to gravity is removed.

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "wallFvPatch.H"
#include "fvsc.H"
//#include "twoPhaseIcoQGDThermo.H"
#include "shallowWaterQGDThermo.H"
#include "QGDInterpolate.H"
#include "fvOptions.H"

//#include "fvCFD.H"
////#include "QGD.H"
//#include "fvOptions.H"
//#include "shallowWaterQGDThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for inviscid shallow-water equations with rotation"
    );

    #define NO_CONTROL

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFaceFields.H"
    #include "createFaceFluxes.H"
    #include "createFvOptions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	//Info << Foam::qgd::QGDCoeffs::hQGD() << endl;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
		/*
         *
         * Update fields
         *
         */
        #include "updateFields.H"
        /*
         *
         * Update fluxes
         *
         */
        #include "updateFluxes.H"

        /*
         *
         * Update time step
         *
         */
	
		
        Info<< "\n Time = " << runTime.timeName() << nl << endl;

        //#include "CourantNo.H"
		
        
		solve
		(
			fvm::ddt(h)
			+
			fvc::div(phiJm)
		);		
		
		/*
		solve
		(
			fvm::ddt(hU)
			+
			fvc::div(phiJmU)
			+ 
			0.5 * magg * hStar * 
			(
				fvc::grad(ksi)
			)
			-
			magg * tau * fvc::div(hU) * fvc::grad(b)
			-
			fvc::div(phiPi)
		);*/
		
		
		solve
		(
			fvm::ddt(hU)
			+
			fvc::div(phiJmU)
			+ 
			0.5 * magg * fvc::div(phih2)
			+ 
			magg * hStar * fvc::grad(b)
			-
			fvc::div(phiPi)
			- 
			NS * divPiNS
		);
		
        if (!dryZoneCondition)
        {
            U == hU/h;
        }
        else
        {
            U == hU/max(h,dimensionedScalar("h0", dimLength, eps0));
            forAll(U,celli)
            {
                if (h[celli] <= epsilon[celli]) 
                {
                    U[celli] = Zero;
                    hU[celli] = Zero;
                }
            }
        }        
		
        
		ksi == b + h;
		
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
