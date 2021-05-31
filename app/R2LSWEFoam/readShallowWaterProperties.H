Info<< "\nReading shallowWaterProperties" << endl;

IOdictionary shallowWaterProperties
(
    IOobject
    (
        "shallowWaterProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

scalar alpha1(readScalar(shallowWaterProperties.lookup("alpha1")));
scalar alpha2(readScalar(shallowWaterProperties.lookup("alpha2")));

Switch dryZoneCondition = false;
scalar eps1 = 1e-6;
scalar eps2 = 1e-6;
scalar tauU = 0;

scalar r(readScalar(shallowWaterProperties.lookup("r")));

 autoPtr<Foam::qgd::QGDCoeffs>  pCoeffs
(
    Foam::qgd::QGDCoeffs::New(shallowWaterProperties.subDict("QGD").get<word>("QGDCoeffs"), mesh, shallowWaterProperties.subDict("QGD"))
);
Foam::qgd::QGDCoeffs& Coeffs = pCoeffs();

Switch wellBalancedScheme = false;
scalar NS = 0.0;

if (shallowWaterProperties.found("wellBalancedScheme"))
{
	shallowWaterProperties.lookup("wellBalancedScheme") >> wellBalancedScheme;
}

if (shallowWaterProperties.found("NS"))
{
	shallowWaterProperties.lookup("NS") >> NS;
}

if (shallowWaterProperties.found("dryZoneCondition"))
{
	shallowWaterProperties.lookup("dryZoneCondition") >> dryZoneCondition;
}

if (shallowWaterProperties.found("eps1"))
{
	shallowWaterProperties.lookup("eps1") >> eps1;
}

if (shallowWaterProperties.found("eps2"))
{
	shallowWaterProperties.lookup("eps2") >> eps2;
}
if (shallowWaterProperties.found("tauU"))
{
	shallowWaterProperties.lookup("tauU") >> tauU;
}

//Coeffs.correct(NULL);

//Info << Coeffs.alphauQGD() << endl;


//const scalar alpha("alpha", shallowWaterProperties);

