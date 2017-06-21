    Info<< "Reading field p_rgh\n" << endl;

    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

// add by zhijie
    Info<< "Reading field c_CO2\n" << endl;
    volScalarField c_CO2
    (
        IOobject
        (
            "c_CO2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field c_RNH2\n" << endl;
    volScalarField c_RNH2
    (
        IOobject
        (
            "c_RNH2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field c_RNHCOO\n" << endl;
    volScalarField c_RNHCOO
    (
        IOobject
        (
            "c_RNHCOO",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading physicalProperties\n" << endl;
    
    IOdictionary physicalProperties
    (
	IOobject
	(
	    "physicalProperties",
	    runTime.constant(),
	    mesh,
	    IOobject::MUST_READ_IF_MODIFIED,
	    IOobject::NO_WRITE
	)
    );

   Info<< "Reading diffusivity D1L\n" << endl;

   dimensionedScalar D1L
   (
//diffusivity of CO2 in solvent phase
	physicalProperties.lookup("D1L")
   );

   dimensionedScalar D1G
   (
//diffusivity of CO2 in gas phase	
	physicalProperties.lookup("D1G")
   );


   Info<< "Reading diffusivity D2L\n" << endl;
//for MEA RNH2
   dimensionedScalar D2L
   (
	physicalProperties.lookup("D2L")
   );


   dimensionedScalar D2G
   (
	physicalProperties.lookup("D2G")
   );


   Info<< "Reading diffusivity D2L\n" << endl;
//for product RNHCOO-
   dimensionedScalar D3L
   (
	physicalProperties.lookup("D3L")
   );


   dimensionedScalar D3G
   (
	physicalProperties.lookup("D3G")
   );


   dimensionedScalar RConst1
//for reaction kinetic constant R
   (
	physicalProperties.lookup("RConst1")
   );

   dimensionedScalar Henry_const
   (
	physicalProperties.lookup("Henry_const")
   );


//   dimensionedScalar Const2
//for reaction 
//   (
//	physicalProperties.lookup("Const2")
//   );


	

//add by zhijie


    #include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi);

    volScalarField& alpha1(twoPhaseProperties.alpha1());

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();


    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + (scalar(1) - alpha1)*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();


    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );


    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);


    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    #include "readGravitationalAcceleration.H"

    /*
    dimensionedVector g0(g);

    // Read the data file and initialise the interpolation table
    interpolationTable<vector> timeSeriesAcceleration
    (
        runTime.path()/runTime.caseConstant()/"acceleration.dat"
    );
    */

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    fv::IOoptionList fvOptions(mesh);

// add by zhijie
    Info<< "Initializing D1 \n" << endl;	
    volScalarField D1
    (
	IOobject
	(
	    "D1",
	    runTime.timeName(),
	    mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
//	(scalar(1)-alpha1)*D1G+alpha1*D1L
 	D1L*D1G/(D1G*alpha1+D1L*(1-alpha1)) 
   );	

    Info<< "Initializing D2 \n" << endl;	
    volScalarField D2
    (
	IOobject
	(
	    "D2",
	    runTime.timeName(),
	    mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	alpha1*D2L
    );	


    Info<< "Initializing D3 \n" << endl;	
    volScalarField D3
    (
	IOobject
	(
	    "D3",
	    runTime.timeName(),
	    mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
        alpha1*D3L
    );	

    volScalarField Gama
    (
        IOobject
        (
            "Gama",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
//	D1*c_CO2*(1-Henry_const)/(alpha1+Henry_const*(1-alpha1)),
	dimensionedScalar("Gama", dimensionSet(0,-1,-1,0,1,0,0), scalar(0.0))
//	dimensionedScalar("Gama", dimensionSet(0,-1,-1,0,1,0,0))
    );

/*
    Info<< "Calculating field alphaFlux1\n" << endl;	
    surfaceScalarField alphaFlux
    (
	IOobject
	(
	    "alphaFlux",
	    runTime.timeName(),
	    mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	mesh.magSf()*fvc::snGrad(alpha1)
    );	
*/


/*
    Info<< "Calculating field ephiFlux1\n" << endl;	
    surfaceScalarField ephiFlux1
    (
	IOobject
	(
	    "ephiFlux1",
	    runTime.timeName(),
	    mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	-(linearInterpolate(Dc1))*mesh.magSf()*fvc::snGrad(ephi)
    );	

    Info<< "Calculating field ephiFlux2\n" << endl;	
    surfaceScalarField ephiFlux2
    (
	IOobject
	(
	    "ephiFlux2",
	    runTime.timeName(),
	    mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	-(linearInterpolate(Dc2))*mesh.magSf()*fvc::snGrad(ephi)
    );	

    Info<< "Initializing surface charge\n" << endl;	
    volScalarField c_surf
    (
	IOobject
	(
	    "c_surf",
	    runTime.timeName(),
	    mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
	),
	mesh
    );	

    Info<< "Initializing d1_c1\n" << endl;	
    volScalarField d1_c1
    (
	IOobject
	(
	    "d1_c1",
	    runTime.timeName(),
	    mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("d1_c1", dimensionSet(0,-1,-1,0,1,0,0), scalar(0.0))
    );	

    Info<< "Initializing d2_c2\n" << endl;	
    volScalarField d2_c2
    (
	IOobject
	(
	    "d2_c2",
	    runTime.timeName(),
	    mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("d2_c2", dimensionSet(0,-1,-1,0,1,0,0), scalar(0.0))
    );	


    Info<< "Initializing field c1_c2\n" << endl;
    volScalarField c1_c2
    (
        IOobject
        (
            "c1_c2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("c1_c2", dimensionSet(0,-3,0,0,1,0,0), scalar(0.0))
    );

    Info<< "Initializing Random Concentration Flux Field J1 for ion 1\n" << endl;
    volVectorField J1
    (
        IOobject
        (
            "J1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
       mesh,
       dimensionedVector("J1", dimensionSet(0,-2,-1,0,1,0,0), vector(0.0,0.0,0.0))
    );

    Info<< "Initializing Random Concentration Flux Field J2 for ion 1\n" << endl;
    volVectorField J2
    (
        IOobject
        (
            "J2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("J2", dimensionSet(0,-2,-1,0,1,0,0), vector(0.0,0.0,0.0))
    );
*/
