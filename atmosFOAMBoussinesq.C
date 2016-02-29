/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Application
	atmosFOAMBuossinesq
	 
	0.4.0 
	-----
		* Reading the coriolis force, and the geostrphic wind. if not found use 0. 
	 
	0.3.0 
	------ 
		* Include the geostrophic pressures. 

Description
	A basic buossinesq solver for the atmosphere. 
>>>>>>> laptop_OF2.2.2/masterLaptop

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "radiationModel.H"
#include "turbulentTransportModel.H"

#include "fvIOoptionList.H"
#include "pimpleControl.H"
#include "interpolation.H"
#include "Random.H"
#include "meshSearch.H"

#include <sstream>
#include <fstream>



//#include "stdlib.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createIncompressibleRadiationModel.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    Info << " \n\nAtmospheric solver  :  0.4.0" << endl; 
    Info << " -------------------------- " << endl; 

    const bool nonlinear = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault("nonlinearSolver", true);
	Info << "\t\t " << (nonlinear ? "Nonlinear " : "Linear ") << " solver. " << endl;
	
	const bool forcepreference = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault("forcepreference", false);
	Info << "\t\t " << (forcepreference ? "force p reference " : "") << endl;
	
	Info << "\t\t " << "Using coriolis vector " << fCoriolis << endl; 
	
    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
	

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
		
           #include "UEqn.H"
	       #include "TEqn.H"
	    
	        // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
	
		Ttotal = T+Tmean;

		runTime.write();

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;
            
    }

    Info<< "End\n" << endl;

    return 0;
}

