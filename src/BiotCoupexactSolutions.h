#pragma once

#include "xEval.h"
#include "xGeomElem.h"

class PEM_Biot
{
public:
  PEM_Biot(double omega_, double beta_,
  double phi_, double sigma_, double alpha_, double lambda_prime_, double lambda_, double rho_1_,
  double nu_, double E_, double eta_);

  using valtype = std::complex<double>;
  const valtype j{0.,1.};
  double omega, beta, k_a, ky;
  double phi, sigma, alpha, lambda_prime, lambda, rho_1, nu, E, eta;
  valtype N, P_hat, A_hat, delta_1, delta_2, delta_3, mu_1, mu_2, mu_3;
  valtype K_eq_til, rho_eq_til;


};


class xEvalCoupCoefficient
{
public:
    xEvalCoupCoefficient(double omega_, double beta_, double sigmaF_, double d_);
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;

    valtype R1, R2, R3, T1, T2, T3, A0, A11, A12, A21, A22;
    valtype kx_11, kx_12, kx_13, kx_21, kx_22, kx_23;
    valtype mat1_N, mat1_A_hat, mat2_N, mat2_A_hat;
    valtype mat1_delta_1, mat1_delta_2, mat1_delta_3;
    valtype mat2_delta_1, mat2_delta_2, mat2_delta_3;
    valtype mat1_mu_1, mat1_mu_2, mat1_mu_3;
    valtype mat2_mu_1, mat2_mu_2, mat2_mu_3;
    valtype mat1_rho_eq_til, mat1_K_eq_til, mat2_rho_eq_til, mat2_K_eq_til;
    double ky;
    
};

class xEvalExactCoupBiotUs : public xfem::xEval<xtensor::xVector<std::complex<double> > > {

public:
    xEvalExactCoupBiotUs(double omega_, double beta_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;
    
private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;;
    valtype A0, A11, A12, A21, A22, R3, T3;
    valtype kx_11, kx_12, kx_13, kx_21, kx_22, kx_23;
    valtype delta_1, delta_2, delta_3;
    double ky;
    
};


class xEvalExactCoupBiotTraction : public xfem::xEval<xtensor::xTensor2<std::complex<double>> > {

public:
    xEvalExactCoupBiotTraction(double omega_, double beta_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xTensor2<std::complex<double> >& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;;
    valtype A0, A11, A12, R3, A21, A22, T3;
    valtype kx_11, kx_12, kx_13, kx_21, kx_22, kx_23;
    valtype mat1_N, mat1_A_hat, mat2_N, mat2_A_hat;
    double ky;
};


class xEvalExactCoupBiotUtotal : public xfem::xEval<xtensor::xVector<std::complex<double> > > {

public:
    xEvalExactCoupBiotUtotal(double omega_, double beta_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;;
    valtype A0, A11, A12, A21, A22, R3, T3;
    valtype mat1_mu_1, mat1_mu_2, mat1_mu_3;
    valtype mat2_mu_1, mat2_mu_2, mat2_mu_3;
    valtype kx_11, kx_12, kx_13, kx_21, kx_22, kx_23;
    double ky;
};


class xEvalExactCoupBiotPressure : public xfem::xEval<std::complex<double> > {

public:
    xEvalExactCoupBiotPressure(double omega_, double beta_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;;
    double ky;
    valtype R1, R2, R3, T1, T2, T3;
    valtype kx_11, kx_12, kx_13, kx_21, kx_22, kx_23;
    valtype mat1_delta_1, mat1_delta_2, mat1_delta_3;
    valtype mat2_delta_1, mat2_delta_2, mat2_delta_3;
    valtype mat1_mu_1, mat1_mu_2, mat1_mu_3;
    valtype mat2_mu_1, mat2_mu_2, mat2_mu_3;
    valtype mat1_K_eq_til, mat2_K_eq_til;

    
};


