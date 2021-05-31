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
hf = qgdInterpolate(h);

// Velocity
Uf = qgdInterpolate(U);

// Momentum
hUf = qgdInterpolate(hU);
//hf*Uf;

// Nonlinear Tensor
UhUf = hf*Uf*Uf;
//qgdInterpolate(U*hU);

if (!dryZoneCondition)
{
    tau = alpha*hQGD_/sqrt(magg*h);
}
else
{
    tau = alpha*hQGD_/(sqrt(magg*max(h,dimensionedScalar("h0", dimLength, eps0))) + tauU*mag(U));
    forAll(tau, celli)
    {
        if (h[celli] <= epsilon[celli]) 
        {
            tau[celli] = 0;
        }
    }
} 

tauf = qgdInterpolate(tau);
//tauf = alpha * hQGDf_ / sqrt(magg * hf);

if (wellBalancedScheme)
{
	hSum = fvc::surfaceIntegrate(hf * mesh.Sf().component(vector::X));

	hhSum = 0.5*fvc::surfaceIntegrate(hf * hf * mesh.Sf().component(vector::X));

	forAll(h_,celli)
	{

		if (mag(hSum[celli]) > 1e-15) //0.5*h[celli]
		{
			h_[celli] = hhSum[celli] / hSum[celli];
		}
		else
		{
			h_[celli] = h[celli];
		}
	}
	
	hStar = h_;
}
else
{
	hStar = h;
}

nuf = tauf * magg * hf;

divPiNS = fvc::laplacian(nuf * hf,U)
        +
		fvc::div(nuf * hf * mesh.Sf() & qgdInterpolate(Foam::T(fvc::grad(U))));

//hStar -= tau * fvc::div(hU);
//hStar -= tau * fvc::div(hUf & mesh.Sf());


// Nonlinear Tensor
//UrhoUf = qgdInterpolate(U*rhoU);

//Speed of sound
//cf = qgdInterpolate(thermo.c());

// ************************************************************************* //