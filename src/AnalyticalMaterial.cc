#include "AnalyticalMaterial.h"
// #include "xEval.h"
#include <iomanip>

using namespace std;
// using namespace xfem;


Analytical_MatBiot::Analytical_MatBiot(double omega_, double beta_, double phi_, double sigma_, double alpha_, double lambda_prime_, double lambda_, 
double rho_1_, double nu_, double E_, double eta_) : omega(omega_), beta(beta_),
phi(phi_), sigma(sigma_), alpha(alpha_), lambda_prime(lambda_prime_),
lambda(lambda_), rho_1(rho_1_), nu(nu_), E(E_), eta(eta_)

{
  const double P = 101325., gamma = 1.4, rho_a = 1.213, mu = 0.1839e-4, Pr = 0.710;

  double K_a = gamma * P;
  double c_a = sqrt(K_a/rho_a);
  k_a = omega/c_a;
  double Z_a = rho_a * c_a;
  double nu_p = mu / rho_a;
  double nu_prime = nu_p / Pr;

  valtype omega_0 = sigma * phi / (rho_a * alpha);
  valtype omega_inf = pow(sigma, 2)*pow(phi,2)*pow(lambda,2)/(4.*mu*rho_a*pow(alpha,2));
  rho_eq_til  = (rho_a*alpha/phi)*(1.+((omega_0)/(j*omega))*sqrt(1.+j*omega/omega_inf));

  valtype omega_prime_infty = (16.*nu_prime)/(pow(lambda_prime, 2));
  valtype F_prime_CA = sqrt(1.+j*omega/omega_prime_infty);
  valtype alpha_prime_til = 1.+omega_prime_infty*F_prime_CA/(2.*j*omega);
  K_eq_til = ((gamma*P/phi)/(gamma-(gamma-1.)/alpha_prime_til));
  valtype c_eq = sqrt(K_eq_til/rho_eq_til);
  k_eq_til = omega / c_eq;
//   cout<<"Wave number"<<k_eq_til<<endl;

  valtype loss = 1. + j*eta;
//   valtype loss = 1.;
  valtype rho_12 = -phi*rho_a*(alpha-1.);
  valtype rho_11 = rho_1-rho_12;
  valtype rho_2 = phi*rho_a;
  valtype rho_22 = rho_2-rho_12;

  valtype rho_22_til = pow(phi,2)*rho_eq_til;
  valtype rho_12_til = rho_2-rho_22_til;
  valtype rho_11_til = rho_1-rho_12_til;
  rho_til = rho_11_til-(pow(rho_12_til, 2)/rho_22_til);
  gamma_til = phi*(rho_12_til/rho_22_til-(1.-phi)/phi);
  rho_s_til =  rho_til+pow(gamma_til,2)*rho_eq_til;

  N = E*loss/(2.*(1.+nu));
  A_hat = nu*E*loss/((1.+nu)*(1.-2.*nu));
  P_hat = A_hat + 2.*N;

  valtype R_til = K_eq_til * pow(phi,2);
  valtype Q_til = ((1 - phi)/phi)*R_til;
  valtype P_til = P_hat + pow(Q_til,2)/R_til;

  valtype delta_eq = omega * sqrt(rho_eq_til/K_eq_til);
  valtype delta_s1 = omega * sqrt(rho_til/P_hat);
  valtype delta_s2 = omega * sqrt(rho_s_til/P_hat);

  valtype Psi = pow((pow(delta_s2,2)+pow(delta_eq,2)),2)-4.*pow(delta_eq,2)*pow(delta_s1,2);
  valtype sdelta_total = sqrt(Psi);
  delta_1 = sqrt(0.5*(pow(delta_s2,2)+pow(delta_eq,2)+sdelta_total));
  delta_2 = sqrt(0.5*(pow(delta_s2,2)+pow(delta_eq,2)-sdelta_total));
  delta_3 = omega*sqrt(rho_til/N);
    // cout<<"wave number is: "<<delta_1<<delta_2<<delta_3<<endl;

  mu_1 = gamma_til*pow(delta_eq,2)/(pow(delta_1,2)-pow(delta_eq,2));
  mu_2 = gamma_til*pow(delta_eq,2)/(pow(delta_2,2)-pow(delta_eq,2));
  mu_3 = -gamma_til;

}


Analytical_MatLimp::Analytical_MatLimp(double omega_, double beta_, double phi_, double sigma_, double alpha_, double lambda_prime_, double lambda_, 
double rho_1_, double nu_, double E_, double eta_) : omega(omega_), beta(beta_),
phi(phi_), sigma(sigma_), alpha(alpha_), lambda_prime(lambda_prime_),
lambda(lambda_), rho_1(rho_1_), nu(nu_), E(E_), eta(eta_)

{
  const double P = 101325., gamma = 1.4, rho_a = 1.213, mu = 0.1839e-4, Pr = 0.710;

  double K_a = gamma * P;
  double c_a = sqrt(K_a/rho_a);
  k_a = omega/c_a;
  double Z_a = rho_a * c_a;
  double nu_p = mu / rho_a;
  double nu_prime = nu_p / Pr;

  valtype omega_0 = sigma * phi / (rho_a * alpha);
  valtype omega_inf = pow(sigma, 2)*pow(phi,2)*pow(lambda,2)/(4.*mu*rho_a*pow(alpha,2));
  rho_eq_til  = (rho_a*alpha/phi)*(1.+((omega_0)/(j*omega))*sqrt(1.+j*omega/omega_inf));

  valtype omega_prime_infty = (16.*nu_prime)/(pow(lambda_prime, 2));
  valtype F_prime_CA = sqrt(1.+j*omega/omega_prime_infty);
  valtype alpha_prime_til = 1.+omega_prime_infty*F_prime_CA/(2.*j*omega);
  K_eq_til = ((gamma*P/phi)/(gamma-(gamma-1.)/alpha_prime_til));
  c_eq_til = sqrt(K_eq_til/rho_eq_til);
  k_eq_til = omega / c_eq_til;
//   cout<<"Wave number"<<k_eq_til<<endl;

  valtype loss = 1. + j*eta;
//   valtype loss = 1.;
  valtype rho_12 = -phi*rho_a*(alpha-1.);
  valtype rho_11 = rho_1-rho_12;
  valtype rho_2 = phi*rho_a;
  valtype rho_22 = rho_2-rho_12;

  valtype rho_22_til = pow(phi,2)*rho_eq_til;
  valtype rho_12_til = rho_2-rho_22_til;
  valtype rho_11_til = rho_1-rho_12_til;
  rho_til = rho_11_til-(pow(rho_12_til, 2)/rho_22_til);
  gamma_til = phi*(rho_12_til/rho_22_til-(1.-phi)/phi);
  rho_s_til =  rho_til+pow(gamma_til,2)*rho_eq_til;

  N = E*loss/(2.*(1.+nu));
  A_hat = nu*E*loss/((1.+nu)*(1.-2.*nu));
  P_hat = A_hat + 2.*N;

  valtype R_til = K_eq_til * pow(phi,2);
  valtype Q_til = ((1 - phi)/phi)*R_til;
  valtype P_til = P_hat + pow(Q_til,2)/R_til;

  valtype delta_eq = omega * sqrt(rho_eq_til/K_eq_til);
  valtype delta_s1 = omega * sqrt(rho_til/P_hat);
  valtype delta_s2 = omega * sqrt(rho_s_til/P_hat);

  valtype Psi = pow((pow(delta_s2,2)+pow(delta_eq,2)),2)-4.*pow(delta_eq,2)*pow(delta_s1,2);
  valtype sdelta_total = sqrt(Psi);
  delta_1 = sqrt(0.5*(pow(delta_s2,2)+pow(delta_eq,2)+sdelta_total));
  delta_2 = sqrt(0.5*(pow(delta_s2,2)+pow(delta_eq,2)-sdelta_total));
  delta_3 = omega*sqrt(rho_til/N);
    // cout<<"wave number is: "<<delta_1<<delta_2<<delta_3<<endl;

  mu_1 = gamma_til*pow(delta_eq,2)/(pow(delta_1,2)-pow(delta_eq,2));
  mu_2 = gamma_til*pow(delta_eq,2)/(pow(delta_2,2)-pow(delta_eq,2));
  mu_3 = -gamma_til;

  rho_limp_til = rho_til*rho_eq_til/(rho_til+rho_eq_til*pow(gamma_til,2));


}


Analytical_MatEF::Analytical_MatEF(double omega_, double beta_, double phi_, double sigma_, double alpha_, double lambda_prime_, double lambda_) : omega(omega_), beta(beta_),
phi(phi_), sigma(sigma_), alpha(alpha_), lambda_prime(lambda_prime_),
lambda(lambda_)

{
  const double P = 101325., gamma = 1.4, rho_a = 1.213, mu = 0.1839e-4, Pr = 0.710;

  double K_a = gamma * P;
  double c_a = sqrt(K_a/rho_a);
  k_a = omega/c_a;
  double Z_a = rho_a * c_a;
  double nu_p = mu / rho_a;
  double nu_prime = nu_p / Pr;

  valtype omega_0 = sigma * phi / (rho_a * alpha);
  valtype omega_inf = pow(sigma, 2)*pow(phi,2)*pow(lambda,2)/(4.*mu*rho_a*pow(alpha,2));
  rho_eq_til  = (rho_a*alpha/phi)*(1.+((omega_0)/(j*omega))*sqrt(1.+j*omega/omega_inf));

  valtype omega_prime_infty = (16.*nu_prime)/(pow(lambda_prime, 2));
  valtype F_prime_CA = sqrt(1.+j*omega/omega_prime_infty);
  valtype alpha_prime_til = 1.+omega_prime_infty*F_prime_CA/(2.*j*omega);
  K_eq_til = ((gamma*P/phi)/(gamma-(gamma-1.)/alpha_prime_til));
  c_eq_til = sqrt(K_eq_til/rho_eq_til);
  k_eq_til = omega / c_eq_til;
//   cout<<"Wave number"<<k_eq_til<<endl;

}

