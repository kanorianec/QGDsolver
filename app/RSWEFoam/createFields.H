#include "readGravitationalAcceleration.H"
#include "readShallowWaterProperties.H"



Info<< "Reading field h\n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField h_
(
    IOobject
    (
        "h_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    h
);

volScalarField hStar
(
    IOobject
    (
        "h*",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    h
);

//bathymetry
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
    h*0
);
volScalarField ksi
(
    IOobject
    (
        "ksi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h+b
);
volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h*0
);
/*
Info<< "Reading field h0 if present\n" << endl;
volScalarField h0
(
    IOobject
    (
        "h0",
        runTime.findInstance("polyMesh", "points"),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mesh,
    dimensionedScalar(dimLength, Zero)
);*/

Info<< "Creating field hU\n" << endl;
volVectorField hU
(
    IOobject
    (
        "hU",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    hU/max(h,dimensionedScalar("h0", dimLength, 1e-6))
);

volVectorField divPiNS
(
    IOobject
    (
        "divPiNS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U*mag(U)
);
/*
Info<< "Creating field hTotal for post processing\n" << endl;
volScalarField hTotal
(
    IOobject
    (
        "hTotal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    h+h0
);
hTotal.write();
*/

//Info << hQGDf_ << endl;

//Info << hQGD_ << endl;

volScalarField hQGD_
(
	"hQGD",
	Coeffs.hQGD()
);

volScalarField tau
(
    IOobject
    (
        "tau",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha*hQGD_/sqrt(magg*max(h,dimensionedScalar("h0", dimLength, 1e-6)))
);
volScalarField hSum
(
    IOobject
    (
        "hSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("hSum", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0 )
);

volScalarField hhSum
(
    IOobject
    (
        "hhSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("hhSum", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0 )
);

/*surfaceScalarField tau
(
    "rhof",
    linearInterpolate(rho)
);*/

#include "createPhi.H"

Info<< "Creating Coriolis Force" << endl;
const dimensionedVector F("F", ((2.0*Omega) & gHat)*gHat);

mesh.setFluxRequired(h.name());

if (dryZoneCondition)
{
    Info << "Creating cutoff parameter epsilon" << nl;
    
    forAll(epsilon,cellI)
	{
        
        labelList neighbours = mesh.cellCells()[cellI];
        //Info << "cell neighbours of " << b[cellI] << nl;
        scalar deltaB = 0.0;
        forAll(neighbours,cellJ)
        {
            deltaB = max(deltaB,mag(b[neighbours[cellJ]] - b[cellI]));
            //Info << b[neighbours[cellJ]] << " ";
        }
        epsilon[cellI] = max(eps0,deltaB);
        //max(eps0, deltaB);
	}
}
