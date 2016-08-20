{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn(
		(nonlinear ? fvm::ddt(T) + fvm::div(phi, T) : fvm::ddt(T))
	      + fvc::div(phi,Tmean)
	      - fvm::laplacian(alphaEff, T)
	     ==
	        fvOptions(T)
    ); 


    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    rhok_tag = - beta*T;
}
