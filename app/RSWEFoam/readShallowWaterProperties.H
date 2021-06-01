Info<< "\nReading shallowWaterProperties" << endl;

shallowWaterQGDThermo thermo(mesh);
thermo.correct();

scalar alpha = thermo.alpha();

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

//scalar alpha(readScalar(shallowWaterProperties.lookup("alpha")));
Switch wellBalancedScheme = false;
Switch dryZoneCondition = false;
scalar NS = 0.0;
scalar eps0 = 1e-6;
scalar tauU = 0;

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

if (shallowWaterProperties.found("eps0"))
{
	shallowWaterProperties.lookup("eps0") >> eps0;
}

if (shallowWaterProperties.found("tauU"))
{
	shallowWaterProperties.lookup("tauU") >> tauU;
}

/* autoPtr<Foam::qgd::QGDCoeffs>  pCoeffs
(
    Foam::qgd::QGDCoeffs::New(shallowWaterProperties.subDict("QGD").get<word>("QGDCoeffs"), mesh, shallowWaterProperties.subDict("QGD"))
);
Foam::qgd::QGDCoeffs& Coeffs = pCoeffs();
*/

//Coeffs.correct(NULL);

//Info << Coeffs.alphauQGD() << endl;


//const scalar alpha("alpha", shallowWaterProperties);

