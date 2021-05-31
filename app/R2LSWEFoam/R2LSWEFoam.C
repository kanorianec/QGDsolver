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
    shallowWaterFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for inviscid shallow-water equations with rotation.

    If the geometry is 3D then it is assumed to be one layers of cells and the
    component of the velocity normal to gravity is removed.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "QGD.H"
#include "fvOptions.H"

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
		
		volScalarField oldH1 = h1; 
		volScalarField oldH2 = h2; 
		Info << "I was here!" << nl;
		solve
		(
			fvm::ddt(h1)
			+
			fvc::div(phiJ1m)
		);
		
		solve
		(
			fvm::ddt(h2)
			+
			fvc::div(phiJ2m)
		);	
		
		solve
		(
			fvm::ddt(h1U1)
			+
			fvc::div(phiJ1mU1)
			+ 
			0.5 * magg * fvc::div(phih1h1)
			+ 
			//magg * h1Star * (fvc::grad(b) + r* fvc::div(phih2) - r * fvc::div(addPhih1))
			magg * h1Star * (fvc::grad(b) + r * fvc::div(phih2))
			-
			dry2 * r * magg * oldH1 * fvc::div(addPhih1)
			-
			fvc::div(phiPi1)
			- 
			NS * divPiNS1
		);
        
		
		solve
		(
			fvm::ddt(h2U2)
			+
			fvc::div(phiJ2mU2)
			+ 
			0.5 * magg * fvc::div(phih2h2)
			+ 
			//magg * h2Star * (fvc::grad(ksi1) - fvc::div(addPhih2))
			magg * h2Star * fvc::grad(ksi1)
			//magg * h2Star * fvc::div(phih1)
			-
			dry1 * magg * oldH2 * fvc::div(addPhih2)
			-
			fvc::div(phiPi2)
			- 
			NS * divPiNS2
		);
		
		//U1 == h1U1/h1;
		//U2 == h2U2/h2;
        
        if (!dryZoneCondition)
        {
            U1 == h1U1/h1;
            U2 == h2U2/h2;
        }
        else
        {
            U1 == h1U1/max(h1,dimensionedScalar("h01", dimLength, eps1));
            U2 == h2U2/max(h2,dimensionedScalar("h02", dimLength, eps2));
            forAll(U1,celli)
            {
                labelList neighbours = mesh.cellCells()[celli];
                if (h1[celli] <= epsilon1[celli]) 
                {
                    U1[celli] = Zero;
                    h1U1[celli] = Zero;
                    forAll(neighbours,cellJ)
                    {
                        U1[neighbours[cellJ]] = Zero;
                    }
                }
                
                
                if (h2[celli] <= epsilon1[celli]) 
                {
                    U2[celli] = Zero;
                    h2U2[celli] = Zero;
                    forAll(neighbours,cellJ)
                    {
                        U2[neighbours[cellJ]] = Zero;
                    }
                }
            }
        }   
        
		ksi1 == b + h1;
		ksi2 == b + h2 + h1;

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
