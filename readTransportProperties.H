    singlePhaseTransportModel laminarTransport(U, phi);

    // Thermal expansion coefficient [1/K] 1/Theta0
    dimensionedScalar beta(laminarTransport.lookup("beta"));

    // Laminar Prandtl number
    dimensionedScalar Pr(laminarTransport.lookup("Pr"));

    // Turbulent Prandtl number
    dimensionedScalar Prt(laminarTransport.lookup("Prt"));


	
