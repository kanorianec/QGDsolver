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

surfaceTensorField gradU1f
(
    "gradU1f",
    fvsc::grad(U1)
);

surfaceTensorField gradU2f
(
    "gradU2f",
    fvsc::grad(U2)
);

surfaceScalarField divU1f
(
    "divU1f",
    Foam::tr(gradU1f)
);

surfaceScalarField divU2f
(
    "divU2f",
    Foam::tr(gradU2f)
);

surfaceScalarField divh1U1f
(
    "divh1U1f",
    Foam::tr(fvsc::grad(h1U1))
);

surfaceScalarField divh2U2f
(
    "divh2U2f",
    Foam::tr(fvsc::grad(h2U2))
);

surfaceVectorField gradh1f
(
    "gradh1f",
    fvsc::grad(h1)
);

surfaceVectorField gradh2f
(
    "gradh2f",
    fvsc::grad(h2)
);

surfaceVectorField gradksi1f
(
    "gradksi1f",
    fvsc::grad(h1)
);

surfaceVectorField gradksi2f
(
    "gradksi2f",
    fvsc::grad(h2)
);

//---------End---------

//Continuity equation fluxes
//---------Start---------

surfaceVectorField h1U1f
(
    "h1U1f",
    linearInterpolate(h1*U1)
);

surfaceVectorField h2U2f
(
    "h2U2f",
    linearInterpolate(h2*U2)
);

surfaceVectorField h1W1
(
    "h1W1",
    h1U1f 
);

surfaceVectorField h2W2
(
    "h2W2",
    h2U2f 
);

surfaceVectorField j1m
(
    "j1m",
    h1U1f - h1W1
);

surfaceVectorField j2m
(
    "j2m",
    h2U2f - h2W2
);

surfaceScalarField phiJ1m
(
    "phiJ1m",
    mesh.Sf() & j1m
);

surfaceScalarField phiJ2m
(
    "phiJ2m",
    mesh.Sf() & j2m
);

//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

surfaceVectorField phiJ1mU1
(
    "phiJ1mU1",
    mesh.Sf() & (j1m * U1f)
);

surfaceVectorField phiJ2mU2
(
    "phiJ2mU2",
    mesh.Sf() & (j2m * U2f)
);

surfaceTensorField U1h1U1f
(
    "U1h1U1f",
    linearInterpolate(U1*h1U1)
);

surfaceTensorField U2h2U2f
(
    "U2h2U2f",
    linearInterpolate(U2*h2U2)
);

//QGD diffusive fluxes
    
surfaceTensorField Pi1f
(
    "Pi1f",
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
    )
);

surfaceTensorField Pi2f
(
    "Pi2f",
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
    )
);

surfaceVectorField phiPi1
(
    "phiPi1",
    mesh.Sf() & Pi1f
);

surfaceVectorField phiPi2
(
    "phiPi2",
    mesh.Sf() & Pi2f
);


//---------End---------


// ************************************************************************* //