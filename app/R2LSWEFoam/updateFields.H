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
    updateFields
Description
    Updates fields.
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Height
h1f = qgdInterpolate(h1);

h2f = qgdInterpolate(h2);

// Velocity
U1f = qgdInterpolate(U1);

U2f = qgdInterpolate(U2);

// Momentum
h1U1f = qgdInterpolate(h1U1);
//h1f*U1f;

h2U2f = qgdInterpolate(h2U2);
//h2f*U2f;

// Nonlinear Tensor
U1h1U1f = U1f*h1f*U1f;
//qgdInterpolate(U1*h1U1);

U2h2U2f = U2f*h2f*U2f;
//qgdInterpolate(U2*h2U2);

if (!dryZoneCondition)
{
    tau1 = alpha1*hQGD_/sqrt(magg*h1);
    tau2 = alpha2*hQGD_/sqrt(magg*h2);
}
else
{
    tau1 = alpha1*hQGD_/(sqrt(magg*max(h1,dimensionedScalar("h01", dimLength, eps1))) + tauU*mag(U1));
    tau2 = alpha2*hQGD_/(sqrt(magg*max(h2,dimensionedScalar("h02", dimLength, eps2))) + tauU*mag(U2));
    
    forAll(tau1, celli)
    {
        if (h1[celli] <= epsilon1[celli]) 
        {
            tau1[celli] = 0;
            /*
            labelList neighbours = mesh.cellCells()[celli];
            forAll(neighbours,cellJ)
            {
                tau1[neighbours[cellJ]] = 0;
            }*/
        }
        if (h2[celli] <= epsilon1[celli]) 
        {
            tau2[celli] = 0;
            /*labelList neighbours = mesh.cellCells()[celli];
            forAll(neighbours,cellJ)
            {
                tau2[neighbours[cellJ]] = 0;
            }*/
        }
    }
} 

//tau1 = alpha1*hQGD_/sqrt(magg*h1);
//tau2 = alpha2*hQGD_/sqrt(magg*h2);

tau1f = qgdInterpolate(tau1);
tau2f = qgdInterpolate(tau2);

//tau1f = alpha1 * hQGDf_ / sqrt(magg * h1f);
//tau2f = alpha2 * hQGDf_ / sqrt(magg * h2f);

h1Star = h1;
h2Star = h2;

if (wellBalancedScheme)
{
	volScalarField h1Sum = fvc::surfaceIntegrate(h1f * mesh.Sf().component(vector::X));
	volScalarField h1h1Sum = 0.5*fvc::surfaceIntegrate(h1f * h1f * mesh.Sf().component(vector::X));
    
    volScalarField h2Sum = fvc::surfaceIntegrate(h2f * mesh.Sf().component(vector::X));
	volScalarField h2h2Sum = 0.5*fvc::surfaceIntegrate(h2f * h2f * mesh.Sf().component(vector::X));

	forAll(h1Star,celli)
	{
		if (mag(h1Sum[celli]) > 1e-15) //0.5*h[celli]
		{
			h1Star[celli] = h1h1Sum[celli] / h1Sum[celli];
		}
        
        if (mag(h2Sum[celli]) > 1e-15) //0.5*h[celli]
		{
			h2Star[celli] = h2h2Sum[celli] / h2Sum[celli];
		}
	}
}

h1Star = h1Star - tau1 * fvc::div(h1U1);
h2Star = h2Star - tau2 * fvc::div(h2U2);

nu1f = tau1f * magg * h1f;

divPiNS1 = fvc::laplacian(nu1f * h1f,U1)
        +
		fvc::div(nu1f * h1f * mesh.Sf() & qgdInterpolate(Foam::T(fvc::grad(U1))));
		
nu2f = tau2f * magg * h2f;

divPiNS2 = fvc::laplacian(nu2f * h2f,U2)
        +
		fvc::div(nu2f * h2f * mesh.Sf() & qgdInterpolate(Foam::T(fvc::grad(U2))));

//tauf = alpha * hQGDf_ / sqrt(magg * hf);

// ************************************************************************* //