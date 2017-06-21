{
	 solve
        (
            fvm::ddt(c_RNHCOO) + fvm::div(phi,c_RNHCOO) - fvm::laplacian(D3,c_RNHCOO) - RConst1*c_CO2*c_RNH2

//            fvc::ddt(c2) + fvc::div(phi,c2) - fvc::laplacian(Dc2,c2) + Const1*fvc::div(ephiFlux2,c2) - fvc::div(J2)

//            fvm::ddt(c2) + fvm::div(phi,c2) - fvm::laplacian(Dc2,c2) + Const1*fvm::laplacian(c2,ephi)
//            fvm::ddt(c2) + fvm::div(phi,c2) - fvm::laplacian(Dc2,c2) 
	);
}
       
