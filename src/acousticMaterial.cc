#include "acousticMaterial.h"
#include "xComplexUtils.h"
#include "heParser.h"


namespace xfem {

//##############################Liquid_material##########################//
xAcousticLiquid::xAcousticLiquid()
{
    properties_signature.register_string("MATERIAL_CLASS");
    properties_signature.register_string("NAME");
    properties_signature.register_scalar("FLUID_DENSITY");
    properties_signature.register_scalar("SOUND_CELERITY");
    properties.setSignature(&properties_signature);
    properties.astring("MATERIAL_CLASS") = "MATERIAL_ACOUSTIC_FLUID";
}

void xAcousticLiquid::checkProperties()
{
    //    throw;
    //check the properties
    fluid_density = properties.scalar("FLUID_DENSITY");
    sound_celerity = properties.scalar("SOUND_CELERITY");

    Z_a = fluid_density * sound_celerity;

    assert(fluid_density > 0.);
    assert(sound_celerity > 0.);
}

//c^2 grad p grad q + w^2 q p
void  xAcousticLiquid::sensitivityTo(const std::string& phys_token,
                                    std::complex<double>& sensitivity) {

    if (phys_token == "density") sensitivity = fluid_density;
    else if (phys_token == "inv_celerity") sensitivity = 1./sound_celerity;
    else if (phys_token == "inv_density"){
        sensitivity = 1. / fluid_density;
    }
    else if (phys_token == "abs_inv_density"){
        sensitivity = std::abs(1. / (fluid_density));
    }
    else if (phys_token == "impedance"){
        sensitivity = 1. / (fluid_density*fluid_density);
    }
    else if (phys_token == "celerity_square_density") {
        sensitivity = sound_celerity * sound_celerity / fluid_density;
    }
    else if (phys_token == "inv_celerity_square_density") {
        sensitivity = 1. / (sound_celerity * sound_celerity * fluid_density);
    }
    else {
        fprintf(stderr, "No sensitivity coded\n"); assert(0); }
    return;
}


//##############################gas_material##########################//
xAcousticAir::xAcousticAir()
{
    properties_signature.register_string("MATERIAL_CLASS");
    properties_signature.register_string("NAME");

    properties_signature.register_scalar("FLUID_DENSITY");
    properties_signature.register_scalar("ATMOS_PRESSURE");
    properties_signature.register_scalar("POLYTROPIC_COEFFICIENT");
    properties_signature.register_scalar("DYNAMIC_VISCOSITY");
    properties.setSignature(&properties_signature);
    properties.astring("MATERIAL_CLASS") = "MATERIAL_ACOUSTIC_FLUID";
}

void xAcousticAir::checkProperties()
{
    //    throw;
    //check the properties
    rho_a = properties.scalar("FLUID_DENSITY");
    P = properties.scalar("ATMOS_PRESSURE");
    gamma = properties.scalar("POLYTROPIC_COEFFICIENT");

    K_a = gamma * P;
    c_a = sqrt(K_a/rho_a);
    Z_a = rho_a * c_a;

    assert(rho_a > 0.);
    // assert(sound_celerity > 0.);
}

//c^2 grad p grad q + w^2 q p
void  xAcousticAir::sensitivityTo(const std::string& phys_token,
                                    std::complex<double>& sensitivity) {

    if (phys_token == "density") sensitivity = rho_a;
    else if (phys_token == "celerity") sensitivity = c_a;
    else if (phys_token == "inv_density"){
        sensitivity = 1. / rho_a;
    //    std::cout<<"Air density : "<<sensitivity<<std::endl;
    }
        else if (phys_token == "abs_inv_density"){
        sensitivity = std::abs(1. / rho_a);
    }
    else if (phys_token == "inv_celerity") {
        sensitivity = 1. / c_a;
    }
    else if (phys_token == "inv_bulk") {
        sensitivity = 1. / K_a;
    }
    else if (phys_token == "impedance") {
        sensitivity = 1. / (c_a*rho_a);
    }
    else if (phys_token == "celerity_square_density") {
        sensitivity = c_a * c_a / rho_a;
    }
    else if (phys_token == "inv_celerity_square_density") {
        sensitivity = 1./ (c_a * c_a * rho_a);
    }
    else {
        fprintf(stderr, "No sensitivity coded\n"); assert(0); }
    return;
}


void xAcousticAir::sensitivityTo(const std::string& phys_token, xtensor::xTensor4<valtype>  &sensitivity)
{
    // const bool debug = xdebug_flag;
//    if (debug) cout << "inside  xElastic::sensitivityTo" << endl;
   { fprintf(stderr, "No sensitivity coded IN AIR\n"); assert(0); }
   return; 

}


//############################## Equivalent Fluid (Porous) (JCA) ##########################//
double xAcousticEquivalentFluid::omega = 0;

xAcousticEquivalentFluid::xAcousticEquivalentFluid()
{
    properties_signature.register_string("MATERIAL_CLASS");
    properties_signature.register_string("NAME");

    properties_signature.register_scalar("THICKNESS");
    properties_signature.register_scalar("AIR_DENSITY");
    properties_signature.register_scalar("ATMOS_PRESSURE");
    properties_signature.register_scalar("POLYTROPIC_COEFFICIENT");
    properties_signature.register_scalar("DYNAMIC_VISCOSITY");
    properties_signature.register_scalar("PRANDTL_NUMBER");
    properties_signature.register_scalar("POROSITY");
    properties_signature.register_scalar("RESISTIVITY");
    properties_signature.register_scalar("TORTUOSITY");
    properties_signature.register_scalar("VISCOUS_CHARAC_LENGTH");
    properties_signature.register_scalar("VISCOUS_LENGTH");
    properties.setSignature(&properties_signature);
    properties.astring("MATERIAL_CLASS") = "MATERIAL_ACOUSTIC_FLUID";
}

void xAcousticEquivalentFluid::checkProperties()
{
    //    throw;
    //check the properties
//    omega = properties.scalar("ANG_FREQUENCY");
//    fluid = properties.scalar("FLUID");
//    porous = properties.scalar("POROUS");
    rho_a = properties.scalar("AIR_DENSITY");
    P = properties.scalar("ATMOS_PRESSURE");
    gamma = properties.scalar("POLYTROPIC_COEFFICIENT");
    mu = properties.scalar("DYNAMIC_VISCOSITY");
    Pr = properties.scalar("PRANDTL_NUMBER");

    d = properties.scalar("THICKNESS");
    phi = properties.scalar("POROSITY");
    sigma = properties.scalar("RESISTIVITY");
    alpha = properties.scalar("TORTUOSITY");
    lambda_prime = properties.scalar("VISCOUS_CHARAC_LENGTH");
    lambda = properties.scalar("VISCOUS_LENGTH");

    nu = mu / rho_a;
    nu_prime = nu / Pr;

    omega_0 = sigma * phi / (rho_a * alpha);
    omega_inf = pow(sigma, 2)*pow(phi,2)*pow(lambda,2)/(4.*mu*rho_a*pow(alpha,2));
    rho_eq  = ((rho_a*alpha/phi)*(1.+((omega_0)/(j*omega))*sqrt(1.+j*omega/omega_inf)));

    omega_prime_infty = (16.*nu_prime)/(pow(lambda_prime, 2));
    F_prime_CA = sqrt(1.+j*omega/omega_prime_infty);
    alpha_prime_til = 1.+omega_prime_infty*F_prime_CA/(2.*j*omega);
    K = ((gamma*P/phi)/(gamma-(gamma-1.)/alpha_prime_til));
    c_eq = sqrt(K/rho_eq);
//    assert(rho_a > 0.);
//    assert(c_eq > 0.);


}

void  xAcousticEquivalentFluid::sensitivityTo(const std::string& phys_token,
                                    std::complex<double> &sensitivity) {
//    std::cout<<fluid_density<<" "<<sound_celerity<<std::endl;
//    std::cout<<rho_eq<<" "<<c_eq<<std::endl;
    if (phys_token == "inv_density") sensitivity = 1./rho_eq;
        else if (phys_token == "abs_inv_density"){
        sensitivity = std::abs(1./rho_eq);
    }
    else if (phys_token == "celerity_square_density"){
        sensitivity = c_eq * c_eq/rho_eq;
    }
    else if (phys_token == "impedance"){
        sensitivity = 1./(c_eq*rho_eq);
    }
    else if (phys_token == "inv_celerity") {
        sensitivity = 1./c_eq ;
    }
    else if (phys_token == "inv_bulk") {
        sensitivity = 1. / K;
    }
    else if (phys_token == "inv_celerity_square_density") {
        sensitivity = 1./(c_eq * c_eq*rho_eq);
    }
    else if (phys_token == "inv_pressure_jump") {
        sensitivity = 1./(sigma * d);
    }
    else {
        fprintf(stderr, "No sensitivity coded\n"); assert(0); }
    return;
}


//############################## Limp Fluid (Porous) (JCA) ##########################//
double xAcousticLimpFluid::omega = 0;

xAcousticLimpFluid::xAcousticLimpFluid()
{
    properties_signature.register_string("MATERIAL_CLASS");
    properties_signature.register_string("NAME");

    properties_signature.register_scalar("THICKNESS");
    properties_signature.register_scalar("AIR_DENSITY");
    properties_signature.register_scalar("ATMOS_PRESSURE");
    properties_signature.register_scalar("POLYTROPIC_COEFFICIENT");
    properties_signature.register_scalar("DYNAMIC_VISCOSITY");
    properties_signature.register_scalar("PRANDTL_NUMBER");
    properties_signature.register_scalar("POROSITY");
    properties_signature.register_scalar("RESISTIVITY");
    properties_signature.register_scalar("TORTUOSITY");
    properties_signature.register_scalar("THERMAL_LENGTH");
    properties_signature.register_scalar("VISCOUS_LENGTH");

    properties_signature.register_scalar("AGGREGATE_DENSITY");
    properties_signature.register_scalar("POISSON_RATIO");
    properties_signature.register_scalar("YOUNG_MODULUS");
    properties_signature.register_scalar("LOSS_FACTOR");
    properties_signature.register_scalar("LOSS_TYPE");
    properties.setSignature(&properties_signature);
    properties.astring("MATERIAL_CLASS") = "MATERIAL_ACOUSTIC_FLUID";

}

void xAcousticLimpFluid::checkProperties()
{
    //    throw;
    //check the properties
//    omega = properties.scalar("ANG_FREQUENCY");
//    fluid = properties.scalar("FLUID");
//    porous = properties.scalar("POROUS");
    rho_a = properties.scalar("AIR_DENSITY");
    P = properties.scalar("ATMOS_PRESSURE");
    gamma = properties.scalar("POLYTROPIC_COEFFICIENT");
    mu = properties.scalar("DYNAMIC_VISCOSITY");
    Pr = properties.scalar("PRANDTL_NUMBER");

    d = properties.scalar("THICKNESS");
    phi = properties.scalar("POROSITY");
    sigma = properties.scalar("RESISTIVITY");
    alpha = properties.scalar("TORTUOSITY");
    lambda_prime = properties.scalar("THERMAL_LENGTH");
    lambda = properties.scalar("VISCOUS_LENGTH");

    rho_1 = properties.scalar("AGGREGATE_DENSITY");
    nu_p = properties.scalar("POISSON_RATIO");
    E = properties.scalar("YOUNG_MODULUS");
    eta = properties.scalar("LOSS_FACTOR");
    loss_type = properties.scalar("LOSS_TYPE");

    nu = mu / rho_a;
    nu_prime = nu / Pr;

    omega_0 = sigma * phi / (rho_a * alpha);
    omega_inf = pow(sigma, 2)*pow(phi,2)*pow(lambda,2)/(4.*mu*rho_a*pow(alpha,2));
    rho_eq_til  = ((rho_a*alpha/phi)*(1.+((omega_0)/(j*omega))*sqrt(1.+j*omega/omega_inf)));

    omega_prime_infty = (16.*nu_prime)/(pow(lambda_prime, 2));
    F_prime_CA = sqrt(1.+j*omega/omega_prime_infty);
    alpha_prime_til = 1.+omega_prime_infty*F_prime_CA/(2.*j*omega);
    K = ((gamma*P/phi)/(gamma-(gamma-1.)/alpha_prime_til));
    c_eq = sqrt(K/rho_eq_til);
//    assert(rho_a > 0.);
//    assert(c_eq > 0.);
   loss = 1. + j*eta ;
    // loss = 1.;
   rho_12 = -phi*rho_a*(alpha-1.);
   rho_11 = rho_1-rho_12;
   rho_2 = phi*rho_a;
   rho_22 = rho_2-rho_12;

   rho_22_til = pow(phi,2)*rho_eq_til;
   rho_12_til = rho_2-rho_22_til;
   rho_11_til = rho_1-rho_12_til;
   rho_til = rho_11_til-(pow(rho_12_til, 2)/rho_22_til);

   gamma_til = phi*(rho_12_til/rho_22_til-(1.-phi)/phi);
   rho_s_til =  rho_til+pow(gamma_til,2)*rho_eq_til;

   rho_limp_til = rho_til*rho_eq_til/(rho_til+rho_eq_til*pow(gamma_til,2));
   


}

void  xAcousticLimpFluid::sensitivityTo(const std::string& phys_token,
                                    std::complex<double> &sensitivity) {
//    std::cout<<fluid_density<<" "<<sound_celerity<<std::endl;
//    std::cout<<rho_eq<<" "<<c_eq<<std::endl;
    if (phys_token == "inv_density") {sensitivity = 1./rho_limp_til;}
    // std::cout<<"rho limp til: "<<rho_limp_til<<std::endl;}
    else if (phys_token == "abs_inv_density"){
        sensitivity = std::abs(1./rho_limp_til);
        
    }
    else if (phys_token == "celerity_square_density"){
        sensitivity = c_eq * c_eq/rho_limp_til;
    }
    else if (phys_token == "impedance"){
        sensitivity = 1./(c_eq*rho_limp_til);
    }
    else if (phys_token == "inv_celerity") {
        sensitivity = 1./rho_limp_til ;
    }
    else if (phys_token == "inv_bulk") {
        sensitivity = 1. / K;
    }
    else if (phys_token == "inv_celerity_square_density") {
        sensitivity = 1./(c_eq * c_eq*rho_limp_til);
    }
    else if (phys_token == "inv_pressure_jump") {
        sensitivity = 1./(sigma * d);
    }
    else {
        fprintf(stderr, "No sensitivity coded\n"); assert(0); }
    return;
}


//############################## Porous_material (BIOT) ##########################//
double xAcousticPEMBiot::omega = 0;
xAcousticPEMBiot::xAcousticPEMBiot()

{
   variables_signature.register_tensor2("strain");
   properties_signature.register_string("MATERIAL_CLASS");
   properties_signature.register_string("NAME");

   properties_signature.register_scalar("AIR_DENSITY");
   properties_signature.register_scalar("ATMOS_PRESSURE");
   properties_signature.register_scalar("POLYTROPIC_COEFFICIENT");
   properties_signature.register_scalar("DYNAMIC_VISCOSITY");
   properties_signature.register_scalar("PRANDTL_NUMBER");

   properties_signature.register_scalar("POROSITY");
   properties_signature.register_scalar("RESISTIVITY");
   properties_signature.register_scalar("TORTUOSITY");
   properties_signature.register_scalar("THERMAL_LENGTH");
   properties_signature.register_scalar("VISCOUS_LENGTH");

   properties_signature.register_scalar("AGGREGATE_DENSITY");
   properties_signature.register_scalar("POISSON_RATIO");
   properties_signature.register_scalar("YOUNG_MODULUS");
   properties_signature.register_scalar("LOSS_FACTOR");
   properties_signature.register_scalar("LOSS_TYPE");
   properties.setSignature(&properties_signature);
   properties.astring("MATERIAL_CLASS") = "MATERIAL_ACOUSTIC_FLUID";
}

void xAcousticPEMBiot::checkProperties()
{

    rho_a = properties.scalar("AIR_DENSITY");
    P = properties.scalar("ATMOS_PRESSURE");
    gamma = properties.scalar("POLYTROPIC_COEFFICIENT");
    mu = properties.scalar("DYNAMIC_VISCOSITY");
    Pr = properties.scalar("PRANDTL_NUMBER");

    //d = properties.scalar("THICKNESS");
    phi = properties.scalar("POROSITY");
    sigma = properties.scalar("RESISTIVITY");
    alpha = properties.scalar("TORTUOSITY");
    lambda_prime = properties.scalar("THERMAL_LENGTH");
    lambda = properties.scalar("VISCOUS_LENGTH");


    nu = mu / rho_a;
    nu_prime = nu / Pr;

    omega_0 = sigma * phi / (rho_a * alpha);
    omega_inf = pow(sigma, 2)*pow(phi,2)*pow(lambda,2)/(4.*mu*rho_a*pow(alpha,2));
    rho_eq_til  = ((rho_a*alpha/phi)*(1.+((omega_0)/(j*omega))*sqrt(1.+j*omega/omega_inf)));

    omega_prime_infty = (16.*nu_prime)/(pow(lambda_prime, 2));
    F_prime_CA = sqrt(1.+j*omega/omega_prime_infty);
    alpha_prime_til = 1.+omega_prime_infty*F_prime_CA/(2.*j*omega);
    K_eq_til = ((gamma*P/phi)/(gamma-(gamma-1.)/alpha_prime_til));
    c_eq = sqrt(K_eq_til/rho_eq_til);
    
   rho_1 = properties.scalar("AGGREGATE_DENSITY");
   nu_p = properties.scalar("POISSON_RATIO");
   E = properties.scalar("YOUNG_MODULUS");
   eta = properties.scalar("LOSS_FACTOR");
   loss_type = properties.scalar("LOSS_TYPE");

   loss = 1. + j*eta ;
    // loss = 1.;
   rho_12 = -phi*rho_a*(alpha-1.);
   rho_11 = rho_1-rho_12;
   rho_2 = phi*rho_a;
   rho_22 = rho_2-rho_12;

   rho_22_til = pow(phi,2)*rho_eq_til;
   rho_12_til = rho_2-rho_22_til;
   rho_11_til = rho_1-rho_12_til;
   rho_til = rho_11_til-(pow(rho_12_til, 2)/rho_22_til);

   gamma_til = phi*(rho_12_til/rho_22_til-(1.-phi)/phi);
   rho_s_til =  rho_til+pow(gamma_til,2)*rho_eq_til;

    SetElasticStiffnessIsotropic (E, nu_p, loss, elastic_stiffness);
}

void  xAcousticPEMBiot::sensitivityTo(const std::string& phys_token, std::complex<double> &sensitivity)
{
    if (phys_token == "inv_density") {sensitivity = 1./rho_eq_til;
    // std::cout<<"density : "<<sensitivity<<std::endl;
    }
    else if (phys_token == "abs_inv_density"){
        sensitivity = std::abs(1./rho_eq_til);
    }
    else if (phys_token == "celerity_square_density"){
        sensitivity = c_eq * c_eq/rho_eq_til;
    }
    else if (phys_token == "inv_bulk"){
        sensitivity = 1./K_eq_til;
    }
    else if (phys_token == "inv_celerity") {
        sensitivity = 1./c_eq ;
    }
    else if (phys_token == "inv_celerity_square_density") {
        sensitivity = 1./(c_eq * c_eq*rho_eq_til);
    }
    else if (phys_token == "rho_til") {
        sensitivity = rho_til;
    }
    else if (phys_token == "gamma_til") {
        sensitivity = gamma_til;
    }
    else {
        fprintf(stderr, "No sensitivity coded\n"); assert(0); }
    return;


}

void xAcousticPEMBiot::sensitivityTo(const std::string& phys_token, xtensor::xTensor4<valtype>  &sensitivity)
{
    // const bool debug = xdebug_flag;
//    if (debug) cout << "inside  xElastic::sensitivityTo" << endl;
   if (phys_token == "strain")       sensitivity = elastic_stiffness;
   else { fprintf(stderr, "No sensitivity coded\n"); assert(0); }
   return; 

}


void xAcousticPEMBiot::SetElasticStiffnessIsotropic(double E, double nu_p, std::complex<double> loss,
					    xtensor::xTensor4< valtype >& stiffness) {

//nu_p*=0.;
const valtype lam = nu_p*E*loss/((1.+nu_p)*(1.-2.*nu_p));
const valtype mu  = E*loss/(2.*(1.+nu_p));


  int i, j, k, l;
  for (i = 0; i < 3; ++i){
    for (j = 0; j < 3; ++j){
      for (k = 0; k < 3; ++k){
	for (l = 0; l < 3; ++l){
	  stiffness(i,j,k,l) = lam * delta(i,j) * delta(k,l) 
	    + mu * (delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k));
	}
      }
    }
  }
  return;
}

}

