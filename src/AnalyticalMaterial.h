#pragma once
// #include "xEval.h"
#include <iomanip>
#include <complex>
using namespace std;



class Analytical_MatBiot
{
public:
  Analytical_MatBiot(double omega_, double beta_,
  double phi_, double sigma_, double alpha_, double lambda_prime_, double lambda_, double rho_1_,
  double nu_, double E_, double eta_);

  using valtype = complex<double>;
  const valtype j{0.,1.};
  double omega, beta, k_a, ky;
  double phi, sigma, alpha, lambda_prime, lambda, rho_1, nu, E, eta;
  valtype gamma_til, N, P_hat, A_hat, delta_1, delta_2, delta_3, mu_1, mu_2, mu_3;
  valtype k_eq_til, K_eq_til, rho_eq_til, rho_til, rho_s_til;

};

class Analytical_MatLimp
{
public:
  Analytical_MatLimp(double omega_, double beta_,
  double phi_, double sigma_, double alpha_, double lambda_prime_, double lambda_, double rho_1_,
  double nu_, double E_, double eta_);

  using valtype = complex<double>;
  const valtype j{0.,1.};
  double omega, beta, k_a, ky;
  double phi, sigma, alpha, lambda_prime, lambda, rho_1, nu, E, eta;
  valtype gamma_til, N, P_hat, A_hat, delta_1, delta_2, delta_3, mu_1, mu_2, mu_3;
  valtype k_eq_til, K_eq_til, rho_eq_til,  c_eq_til, rho_til, rho_s_til, rho_limp_til;

};

class Analytical_MatEF
{
public:
  Analytical_MatEF(double omega_, double beta_,
  double phi_, double sigma_, double alpha_, double lambda_prime_, double lambda_);

  using valtype = complex<double>;
  const valtype j{0.,1.};
  double omega, beta, k_a, ky;
  double phi, sigma, alpha, lambda_prime, lambda;
  valtype k_eq_til, K_eq_til, rho_eq_til, c_eq_til;

};
