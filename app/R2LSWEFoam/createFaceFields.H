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
    createFaceFields
Description
    Creates the face fields: linear interpolation of fields from volumes 
    to face centers.
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*surfaceScalarField hQGDf_
(
	"hQGDf",
	1.0 / mag(mesh.surfaceInterpolation::deltaCoeffs())
);*/

// Velocity
surfaceVectorField U1f
(
    "U1f",
    linearInterpolate(U1)
);

surfaceVectorField U2f
(
    "U2f",
    linearInterpolate(U2)
);

// Height
surfaceScalarField h1f
(
    "h1f",
    linearInterpolate(h1)
);

surfaceScalarField h2f
(
    "h2f",
    linearInterpolate(h2)
);

surfaceScalarField tau1f
(
    "tau1f",
    linearInterpolate(tau1)
);

surfaceScalarField tau2f
(
    "tau2f",
    linearInterpolate(tau2)
);

surfaceVectorField phih1h1
(
    "phih1h1",
    h1f*h1f*mesh.Sf()
);

surfaceVectorField phih2h2
(
    "phih2h2",
    h2f*h2f*mesh.Sf()
);

surfaceVectorField phih2
(
    "phih2",
    h2f*mesh.Sf()
);

surfaceVectorField phih1
(
    "phih1",
    h1f*mesh.Sf()
);

surfaceVectorField addPhih1
(
    "addPhih1",
    h1f*mesh.Sf()
);

surfaceVectorField addPhih2
(
    "addPhih2",
    h2f*mesh.Sf()
);

surfaceScalarField nu1f
(
    "nu1f",
    tau1f*magg*h1f
);
surfaceScalarField nu2f
(
    "nu2f",
    tau1f*magg*h1f
);


// ************************************************************************* //