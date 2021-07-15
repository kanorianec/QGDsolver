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

bf = qgdInterpolate(b);

// Velocity
Uf = qgdInterpolate(U);

// Momentum
hUf = qgdInterpolate(hU);
//hf*Uf;

// Nonlinear Tensor
UhUf = hf*Uf*Uf;

if (!dryZoneCondition)
{
    if (gMin(h) <= 0)
    {
        FatalErrorIn("RSWEFoam.C") << "Can't calculate cases h <=0 without dryZoneCondition = true." << nl << exit(FatalError);
    }
    Info << gMin(h) << endl;
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

//volVectorField gradB = fvc::grad(b);


if (dryZoneCondition)
{
    for( label i=0; i < mesh.nInternalFaces(); i++)
    {
        if (h[own[i]] <= epsilon[own[i]])
        {
            
        }
        //if (h[own[i]] <= epsilon[own[i]] || h[nei[i]] <= epsilon[nei[i]])
        //if (h[own[i]] <= epsilon[own[i]])
        //{
        //    bf[i] =  max(ksi[nei[i]] - hf[i], bf[i]); //h[nei[i]] + b[nei[i]]
        //    //-hf[i]*hf[i]/2;
        //} 
        //if (h[nei[i]] <= epsilon[nei[i]])
        //{
        //    bf[i] =  max(ksi[own[i]] - hf[i], bf[i]); // h[own[i]] + b[own[i]]
        //    //-hf[i]*hf[i]/2;
        //}     
        // some operation of volfield[own[i]] and volfield[nei[i]]
    }
}

vector e1(1, 0, 0);
vector e2(0, 1, 0);
vector e3(0, 0, 1);

volVectorField gradB = fvc::div(bf * mesh.Sf());

if (wellBalancedScheme)
{
    //hSumX = dimensionedScalar(dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0 );
    //hhSumX = dimensionedScalar(dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0 );
    //
    //hSumY = dimensionedScalar(dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0 );
    //hhSumY = dimensionedScalar(dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0 );
    //
    //hSumZ = dimensionedScalar(dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0 );
    //hhSumZ = dimensionedScalar(dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0 );
    
    
    //for( label i=0; i < mesh.nInternalFaces(); i++)
    //{
    //    if (h[nei[i]] > epsilon[nei[i]] && h[own[i]] > epsilon[own[i]])
    //    {
    //        hSumX[nei[i]]-= hf[i] * mesh.Sf()[i].component(vector::X);
    //        hhSumX[nei[i]]-= 0.5* hf[i] * hf[i] * mesh.Sf()[i].component(vector::X);
    //        
    //        hSumY[nei[i]]-= hf[i] * mesh.Sf()[i].component(vector::Y);
    //        hhSumY[nei[i]]-= 0.5* hf[i] * hf[i] * mesh.Sf()[i].component(vector::Y);
    //        
    //        hSumZ[nei[i]]-= hf[i] * mesh.Sf()[i].component(vector::Z);
    //        hhSumZ[nei[i]]-= 0.5* hf[i] * hf[i] * mesh.Sf()[i].component(vector::Z);
    //        
    //        hSumX[own[i]]+= hf[i] * mesh.Sf()[i].component(vector::X);
    //        hhSumX[own[i]]+= 0.5* hf[i] * hf[i] * mesh.Sf()[i].component(vector::X);
    //        
    //        hSumY[own[i]]+= hf[i] * mesh.Sf()[i].component(vector::Y);
    //        hhSumY[own[i]]+= 0.5* hf[i] * hf[i] * mesh.Sf()[i].component(vector::Y);
    //        
    //        hSumZ[own[i]]+= hf[i] * mesh.Sf()[i].component(vector::Z);
    //        hhSumZ[own[i]]+= 0.5* hf[i] * hf[i] * mesh.Sf()[i].component(vector::Z);
    //    }
    //    /*if (h[own[i]]<eps0 || h[nei[i]]<eps0)
    //    {
    //        bf[i] = -hf[i]*hf[i]/2;
    //    } */   
    //    // some operation of volfield[own[i]] and volfield[nei[i]]
    //}
    
	hSumX = fvc::surfaceIntegrate(hf * mesh.Sf().component(vector::X));
	hhSumX = 0.5*fvc::surfaceIntegrate(hf * hf * mesh.Sf().component(vector::X));
    
    hSumY = fvc::surfaceIntegrate(hf * mesh.Sf().component(vector::Y));
	hhSumY = 0.5*fvc::surfaceIntegrate(hf * hf * mesh.Sf().component(vector::Y));
    
    hSumZ = fvc::surfaceIntegrate(hf * mesh.Sf().component(vector::Z));
	hhSumZ = 0.5*fvc::surfaceIntegrate(hf * hf * mesh.Sf().component(vector::Z));

	forAll(h_,celli)
	{
		if (mag(hSumX[celli]) > 1e-15 && h[celli] > epsilon[celli])
		{
			h_[celli].component(vector::X) = hhSumX[celli] / hSumX[celli];           
		}
		else
		{
            h_[celli].component(vector::X) = h[celli];
		}
        if (mag(hSumY[celli]) > 1e-15 && h[celli] > epsilon[celli])
		{
            h_[celli].component(vector::Y ) = hhSumY[celli] / hSumY[celli];
		}
		else
		{
            h_[celli].component(vector::Y) = h[celli];
		}
        if (mag(hSumZ[celli]) > 1e-15 && h[celli] > epsilon[celli])
		{
            h_[celli].component(vector::Z) = hhSumZ[celli] / hSumZ[celli];
		}
		else
		{
            h_[celli].component(vector::Z) = h[celli];
		}
	}
    
    ghGradB =  e1 * magg * h_.component(vector::X) * gradB.component(vector::X) + e2 * magg * h_.component(vector::Y) * gradB.component(vector::Y) + 
            e3 * magg * h_.component(vector::Z) * gradB.component(vector::Z);   
}
else
{
    ghGradB = magg*h*gradB;
}

nuf = tauf * magg * hf;

divPiNS = fvc::laplacian(nuf * hf,U)
        +
		fvc::div(nuf * hf * mesh.Sf() & qgdInterpolate(Foam::T(fvc::grad(U))));

// ************************************************************************* //