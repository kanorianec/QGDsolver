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

if (wellBalancedScheme)
{
	hSumX = fvc::surfaceIntegrate(hf * mesh.Sf().component(vector::X));
	hhSumX = 0.5*fvc::surfaceIntegrate(hf * hf * mesh.Sf().component(vector::X));
    
    hSumY = fvc::surfaceIntegrate(hf * mesh.Sf().component(vector::Y));
	hhSumY = 0.5*fvc::surfaceIntegrate(hf * hf * mesh.Sf().component(vector::Y));

	forAll(h_,celli)
	{
		if (mag(hSumX[celli]) > 1e-15)
		{
			h_[celli].component(vector::X) = hhSumX[celli] / hSumX[celli];           
		}
		else
		{
            h_[celli].component(vector::X) = h[celli];
			//h_[celli] = h[celli];
		}
        if (mag(hSumY[celli]) > 1e-15)
		{
            h_[celli].component(vector::Y) = hhSumY[celli] / hSumY[celli];
		}
		else
		{
            h_[celli].component(vector::Y) = h[celli];
			//h_[celli] = h[celli];
		}
	}
    
    volVectorField gradB = fvc::grad(b);
    
    vector e1(1, 0, 0);
    vector e2(0, 1, 0);
    vector e3(0, 0, 1);
    
    ghGradB =  e1 * magg * h_.component(vector::X) * gradB.component(vector::X) + e2 * magg * h_.component(vector::Y) * gradB.component(vector::Y) + 
            e3 * magg * h * gradB.component(vector::Z);   
    
    
    //ghGradB.value().component(vector::X) = magg * h_.component(vector::X) * gradB.component(vector::X);
    //ghGradB.value().component(vector::Y) = magg * h_.component(vector::Y) * gradB.component(vector::Y);
    //ghGradB.value().component(vector::Z) = magg * h;//* gradB.component(vector::Z);
    
    //volVectorField newGradB = ghGradB / magg;
    
    //Info << "h - hx -> " << mag(h_.component(vector::X) - h) << endl << endl;
	//Info << "h - hy -> " << mag(h_.component(vector::Y) - h) << endl << endl;
	hStar = h_.component(vector::X);
    
    Info << "magg * hStar * fvc::grad(b) - ghGradB -> " << magg * hStar * fvc::grad(b) - ghGradB << endl;
}
else
{
    ghGradB = magg*h*fvc::grad(b);
	hStar = h;
}

nuf = tauf * magg * hf;

divPiNS = fvc::laplacian(nuf * hf,U)
        +
		fvc::div(nuf * hf * mesh.Sf() & qgdInterpolate(Foam::T(fvc::grad(U))));

// ************************************************************************* //