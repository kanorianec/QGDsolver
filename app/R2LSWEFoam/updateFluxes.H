/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd .
    Copyright (C) 2016-2019 ISP RAS (www.ispras.ru) UniCFD Group (www.unicfd.ru)
-------------------------------------------------------------------------------
License
    This file is part of QGDsolver library, based on OpenFOAM+.
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
Global
    updateFluxes
Description
    Updates fluxes for continuity equation, momentum balance equation,
    energy balance equation.
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Gradients and divergence
//---------Start---------
gradU1f = fvsc::grad(U1);

gradU2f = fvsc::grad(U2);

divU1f  = tr(gradU1f);

divU2f  = tr(gradU2f);

divh1U1f  = tr(fvsc::grad(h1U1));

divh2U2f  = tr(fvsc::grad(h2U2));

gradksi1f = fvsc::grad(ksi1);

gradksi2f = fvsc::grad(ksi2);

gradh1f = fvsc::grad(h1);

gradh2f = fvsc::grad(h2);

//---------End---------

//Continuity equation fluxes
//---------Start---------

h1W1 = tau1f *
        (
            (U1f * divh1U1f)
            +
            (h1U1f & gradU1f)
			+
			magg * h1f * (gradksi1f + r * gradh2f)
			//magg * h1f * (gradh1f + r * gradh2f)
        );
		
h2W2 = tau2f *
        (
            (U2f * divh2U2f)
            +
            (h2U2f & gradU2f)
			+
			magg * h2f * (gradh1f + gradh2f)
			//magg * h2f * gradksi2f
        );

j1m = h1f*U1f - h1W1;

j2m = h2f*U2f - h2W2;

phiJ1m = mesh.Sf() & j1m;

phiJ2m = mesh.Sf() & j2m;
//---------End---------
		
// Fluxes for momentum balance equation
//---------Start---------

phiJ1mU1 = qgdFlux(phiJ1m,U1,U1f);

phiJ2mU2 = qgdFlux(phiJ2m,U2,U2f);

phih1h1 = mesh.Sf() * h1f * h1f;

phih2h2 = mesh.Sf() * h2f * h2f;

phih1 = mesh.Sf() * h1f;

phih2 = mesh.Sf() * h2f;

addPhih1 = mesh.Sf() * tau2f * divh2U2f;
addPhih2 = mesh.Sf() * tau1f * divh1U1f;

Pi1f = 
    tau1f*
    (
        (U1h1U1f & gradU1f)
        +
        U1f * magg * h1f * (gradksi1f + r * gradh2f)
    )
    +
    tau1f *
    (
		I * ( magg * h1f * divh1U1f)
		//I * ( magg * h1f * ((U1f & gradh1f) + (h1f * divU1f)))
    );
	
Pi2f =
	tau2f*
	(
		(U2h2U2f & gradU2f)
		+
		U2f * magg * h2f * gradksi2f
	)
	+
	tau2f *
	(
		I * ( magg * h2f * divh2U2f)
		//I * ( magg * h2f * ((U2f & gradh2f) + (h2f * divU2f)))
	);
	
phiPi1 = mesh.Sf() & Pi1f;

phiPi2 = mesh.Sf() & Pi2f;
//---------End---------	
		
// ************************************************************************* //