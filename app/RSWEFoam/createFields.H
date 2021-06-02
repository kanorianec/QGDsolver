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

volVectorField h_
(
    IOobject
    (
        "h_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    h*fvc::grad(b)
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

volScalarField hQGD_ = thermo.hQGD();

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
volScalarField hSumX
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

volScalarField hhSumX
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

volScalarField hSumY
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

volScalarField hhSumY
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

volVectorField ghGradB
(
    IOobject
    (
        "ghgradB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    magg*h*fvc::grad(b)
);


if (dryZoneCondition)
{
    Info << "Creating cutoff parameter epsilon" << nl;
    
    forAll(epsilon,cellI)
	{
        
        labelList neighbours = mesh.cellCells()[cellI];
        scalar deltaB = 0.0;
        forAll(neighbours,cellJ)
        {
            deltaB = max(deltaB,mag(b[neighbours[cellJ]] - b[cellI]));
        }
        epsilon[cellI] = max(eps0,deltaB);
	}
}
