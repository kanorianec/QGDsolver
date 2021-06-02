Info<< "\nReading shallowWaterProperties" << endl;

shallowWaterQGDThermo thermo(mesh);
thermo.correct();

scalar alpha = thermo.alpha();

Switch wellBalancedScheme = thermo.wellBalancedScheme();
Switch dryZoneCondition = thermo.dryZoneCondition();
scalar NS = thermo.NS();
scalar eps0 = thermo.eps0();
scalar tauU = thermo.tauU();