{
        solve
        (

//           fvm::ddt(c_CO2) + fvm::div(phi,c_CO2) - fvm::laplacian(D1,c_CO2) + fvm::div(alpha_flux,Gama) + RConst1*c_CO2*c_RNH2

           fvm::ddt(c_CO2) + fvm::div(phi,c_CO2) - fvm::laplacian(D1,c_CO2) + fvc::laplacian(Gama,alpha1) + RConst1*c_CO2*c_RNH2

//           fvc::ddt(c1) + fvc::div(phi,c1) - fvc::laplacian(Dc1,c1) - Const1*fvc::div(ephiFlux1,c1) - fvc::div(J1)

//           fvm::ddt(c1) + fvm::div(phi,c1) - fvm::laplacian(Dc1,c1) - Const1*fvm::laplacian(d1_c1,ephi)
//           fvm::ddt(c1) + fvm::div(phi,c1) - fvm::laplacian(Dc1,c1) 
	);
}
       
