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
gradUf = fvsc::grad(U);

divUf  = tr(gradUf);

divhUf  = tr(fvsc::grad(hU));

gradksif = fvsc::grad(ksi);

//---------End---------

//Continuity equation fluxes
//---------Start---------
/*hW = tauf *
        (
            ((Uf * gradhf) & Uf)
            +
            (hUf * divUf)
            +
            (hUf & gradUf)
			+
			magg * hf * gradksif
        );*/
		
//Info << "1: " << (Uf * divhUf) << endl;
//Info << "2: " << (((Uf * gradhf) & Uf) + (hUf * divUf)) << endl;

hW = tauf *
        (
            (Uf * divhUf)
            +
            (hUf & gradUf)
			+
			magg * hf * gradksif
        );

jm = hf*Uf - hW;

phiJm = mesh.Sf() & jm;
//---------End---------
		
// Fluxes for momentum balance equation
//---------Start---------

phiJmU = qgdFlux(phiJm,U,Uf);

phih2 = mesh.Sf() * hf * hf;

Pif = 
    //QGD diffusive fluxes
	tauf*
    (
        (UhUf & gradUf)
        +
        Uf * magg * hf * gradksif
    )
    +
    tauf *
    (
		I * ( magg * hf * divhUf)
        //I * ( magg * hf * ((Uf & gradhf) + (hf * divUf)))
    );
	
phiPi = mesh.Sf() & Pif;

//phiw = mesh.Sf() & hW;
//---------End---------	
		
// ************************************************************************* //