/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    createFaceFluxes
Description
    Create Face Fluxes for calculation
SourceFile
    RSWEFoam.C
\*---------------------------------------------------------------------------*/

//Gradients and divergence
//---------Start---------

surfaceTensorField gradUf
(
    "gradUf",
    fvsc::grad(U)
);

surfaceScalarField divUf
(
    "divUf",
    Foam::tr(gradUf)
);

surfaceScalarField divhUf
(
    "divhUf",
    Foam::tr(fvsc::grad(h*U))
);

surfaceVectorField gradksif
(
    "gradksif",
    fvsc::grad(h)
);

//---------End---------

//Continuity equation fluxes
//---------Start---------

surfaceVectorField hUf
(
    "hUf",
    linearInterpolate(h*U)
);

surfaceVectorField hW
(
    "hW",
    hUf ///!!!tauQGDf*gradPf
);

surfaceVectorField jm
(
    "jm",
    hUf - hW
);

surfaceScalarField phiJm
(
    "phiJm",
    mesh.Sf() & jm
);

//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

surfaceVectorField phiJmU
(
    "phiJmU",
    mesh.Sf() & (jm * Uf)
);



surfaceTensorField UhUf
(
    "UhUf",
    linearInterpolate(U*hU)
);

//Info << h << endl;
//
    //+
surfaceTensorField Pif
(
    "Pif",
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
    )
);

/*surfaceTensorField PiNSf
(
    "PiNSf",
    // NS Pi
	Pif
);*/

surfaceVectorField phiPi
(
    "phiPi",
    mesh.Sf() & Pif
);

/*
surfaceScalarField phiw
(
    "phiwStar",
    mesh.Sf() & hW
);
*/
//---------End---------


// ************************************************************************* //