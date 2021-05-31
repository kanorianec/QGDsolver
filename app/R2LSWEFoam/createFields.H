#include "readGravitationalAcceleration.H"
#include "readShallowWaterProperties.H"



Info<< "Reading field h\n" << endl;
volScalarField h1
(
    IOobject
    (
        "h1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField h2
(
    IOobject
    (
        "h2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField h1Star
(
    IOobject
    (
        "h1*",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    h1
);

volScalarField h2Star
(
    IOobject
    (
        "h2*",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    h2
);

//bathymetry file

volScalarField b
(
    IOobject
    (
        "b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    h1*0
);
volScalarField ksi1
(
    IOobject
    (
        "ksi1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h1+b
);
volScalarField ksi2
(
    IOobject
    (
        "ksi2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h1 + h2 + b
);

volScalarField epsilon1
(
    IOobject
    (
        "epsilon1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h1*0
);

volScalarField epsilon2
(
    IOobject
    (
        "epsilon2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h2*0
);

volScalarField dry1
(
    IOobject
    (
        "dry1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

volScalarField dry2
(
    IOobject
    (
        "dry2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

Info<< "Creating field h1U1\n" << endl;
volVectorField h1U1
(
    IOobject
    (
        "h1U1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Creating field h2U2\n" << endl;
volVectorField h2U2
(
    IOobject
    (
        "h2U2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U1\n" << endl;
volVectorField U1
(
    IOobject
    (
        "U1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h1U1/max(h1,dimensionedScalar("h01", dimLength, 1e-6))
);

Info<< "Reading field U2\n" << endl;
volVectorField U2
(
    IOobject
    (
        "U2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h2U2/max(h2,dimensionedScalar("h02", dimLength, 1e-6))
);


volScalarField hQGD_
(
	"hQGD",
	Coeffs.hQGD()
);

surfaceScalarField hQGDf_
(
	"hQGDf",
	Coeffs.hQGDf()
);

volScalarField tau1
(
    IOobject
    (
        "tau1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*hQGD_/sqrt(magg*max(h1,dimensionedScalar("h01", dimLength, 1e-6)))
);

volScalarField tau2
(
    IOobject
    (
        "tau2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha2*hQGD_/sqrt(magg*max(h2,dimensionedScalar("h02", dimLength, 1e-6)))
);

volVectorField divPiNS1
(
    IOobject
    (
        "divPiNS1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U1*mag(U1)
);

volVectorField divPiNS2
(
    IOobject
    (
        "divPiNS2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U2*mag(U2)
);

#include "createPhi.H"

Info<< "Creating Coriolis Force" << endl;
const dimensionedVector F("F", ((2.0*Omega) & gHat)*gHat);

mesh.setFluxRequired(h1.name());

if (dryZoneCondition)
{
    Info << "Creating cutoff parameter epsilon" << nl;
    
    forAll(epsilon1,cellI)
	{
        
        if (h1[cellI] <= epsilon1[cellI])
            dry1[cellI] = 0.0;
        else
            dry1[cellI] = 1.0;
        if (h2[cellI] <= epsilon2[cellI])
            dry2[cellI] = 0.0;
        else
            dry2[cellI] = 1.0;
        
        labelList neighbours = mesh.cellCells()[cellI];
        //Info << "cell neighbours of " << b[cellI] << nl;
        scalar deltaB = 0.0;
        scalar deltah1B = 0.0;
        /*forAll(neighbours,cellJ)
        {
            deltah1B = max(deltaB,mag(b[neighbours[cellJ]] + h1[neighbours[cellJ]] - b[cellI] - h1[cellI]));
            deltaB = max(deltaB,mag(b[neighbours[cellJ]] - b[cellI]));
            //Info << b[neighbours[cellJ]] << " ";
        }*/
        epsilon1[cellI] = max(eps1,deltaB);
        epsilon2[cellI] = max(eps2,deltah1B);
        //max(eps0, deltaB);
	}
}
