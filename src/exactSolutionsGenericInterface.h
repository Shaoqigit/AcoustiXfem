#pragma once

#include <iostream>

#include "xEval.h"
#include "xGeomElem.h"
#include "AnalyticalMaterial.h"

// constexpr unsigned int max_modes = 40;
using valtype = std::complex<double>;

class xEvalInterfaceCoefficient
{
public:
    xEvalInterfaceCoefficient(double omega_, double beta_,  valtype A_, valtype B_, valtype C_, valtype D_);
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta;
    valtype A, B, C, D;
    double A_inc, rho_a, k_a;
    valtype rho_eq, k_eq;
    valtype R, T;
    
};

class xSolObliquePlaneGenericInterface : public xfem::xEval<std::complex<double> > {

public:
    xSolObliquePlaneGenericInterface(double omega_, double beta_, valtype A_, valtype B_, valtype C_, valtype D_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const;

private:
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta;
    valtype A, B, C, D;
    valtype rho_a, rho_eq, k_a, k_eq;
    valtype R, T;
};


class xGradSolObliquePlaneGenericInterface : public xfem::xEval<xtensor::xVector<std::complex<double> > > {

public:
    xGradSolObliquePlaneGenericInterface(double omega_, double beta_, valtype A_, valtype B_, valtype C_, valtype D_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;

private:
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta;
    valtype A, B, C, D;
    valtype rho_a, rho_eq, k_a, k_eq;
    valtype R, T;

};


class xEvalTubeCoeff
{
    public:
    xEvalTubeCoeff(double omega_, double v0_, valtype M11_, valtype M12_, valtype M21_, valtype M22_);

    const valtype j{0.,1.};
    double omega, v0;
    valtype A, B, C, D;
    valtype rho_a, k_a;
    valtype rho_eq, k_eq;
    valtype M11, M12, M21, M22;
};


class xPressureTubeNormal : public xfem::xEval<std::complex<double> > {

public:
    xPressureTubeNormal(double omega_, double v0_, valtype M11_, valtype M12_, valtype M21_, valtype M22_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const;

private:
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, v0;
    valtype A, B, C, D;
    valtype rho_a, k_a;
    valtype rho_eq, k_eq;
    valtype M11, M12, M21, M22;
};


class xGradPressTubeNormal : public xfem::xEval<xtensor::xVector<std::complex<double> > > {

public:
    xGradPressTubeNormal(double omega_, double v0_, valtype M11_, valtype M12_, valtype M21_, valtype M22_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;

private:
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, V0;
    valtype A, B, C, D;
    valtype rho_a, k_a;
    valtype rho_eq, k_eq;
    valtype M11, M12, M21, M22;

};


class xEvalCoeffObliqueRigidWall
{
    public:
    xEvalCoeffObliqueRigidWall(double omega_, double beta_, valtype M11_, valtype M12_, valtype M21_, valtype M22_);

    const valtype j{0.,1.};
    double omega, beta;
    valtype A, B, C;
    valtype rho_a, k_a;
    valtype rho_eq, k_eq;
    valtype M11, M12, M21, M22;
};


class xEvalPressureObliqueRigidWall : public xfem::xEval<std::complex<double> > {

public:
    xEvalPressureObliqueRigidWall(double omega_, double beta_, valtype M11_, valtype M12_, valtype M21_, valtype M22_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const;

private:
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta;
    valtype A, B, C;
    valtype rho_a, k_a;
    valtype rho_eq, k_eq;
    valtype M11, M12, M21, M22;
};


class xEvalGradPressObliqueRigidWall : public xfem::xEval<xtensor::xVector<std::complex<double> > > {

public:
    xEvalGradPressObliqueRigidWall(double omega_, double beta_, valtype M11_, valtype M12_, valtype M21_, valtype M22_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;

private:
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta;
    valtype A, B, C;
    valtype rho_a, k_a;
    valtype rho_eq, k_eq;
    valtype M11, M12, M21, M22;

};


class xEvalInterfaceCoefficientAirBiot
{
public:
    xEvalInterfaceCoefficientAirBiot(double omega_, double beta_,  Analytical_MatEF mat_F);
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta;
    valtype A, B, C, D;
    double A_inc, rho_a, k_a;
    valtype rho_eq_til, K_eq_til;
    valtype R, T1, T2, T3, A_1, A_2;
    valtype N, A_hat, P_hat;
    valtype delta_1, delta_2, delta_3;
};

class xPressureAirBiot
{
public:
    xPressureAirBiot(double omega_, double beta_,  Analytical_MatEF mat_F);
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta;
    valtype A, B, C, D;
    double A_inc, rho_a, k_a;
    valtype rho_eq_til, K_eq_til;
    valtype R, T1, T2, T3, A_1, A_2;
    valtype N, A_hat, P_hat;
    valtype delta_1, delta_2, delta_3;
};


class xGradPAirBiot
{
public:
    xGradPAirBiot(double omega_, double beta_,  Analytical_MatEF& mat_F);
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta;
    valtype A, B, C, D;
    double A_inc, rho_a, k_a;
    valtype rho_eq_til, K_eq_til;
    valtype R, T1, T2, T3, A_1, A_2;
    valtype N, A_hat, P_hat;
    valtype delta_1, delta_2, delta_3;

private:
    Analytical_MatEF& mat_F;
};


class xUtotalAirBiot
{
public:
    xUtotalAirBiot(double omega_, double beta_,  Analytical_MatEF& mat_F_);
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta;
    valtype A, B, C, D;
    double A_inc, rho_a, k_a;
    valtype rho_eq_til, K_eq_til;
    valtype R, T1, T2, T3, A_1, A_2;
    valtype N, A_hat, P_hat;
    valtype delta_1, delta_2, delta_3;

private:
    Analytical_MatEF& mat_F;

};


class xTractionAirBiot
{
public:
    xTractionAirBiot(double omega_, double beta_,  Analytical_MatEF& mat_F_);
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta;
    valtype A, B, C, D;
    double A_inc, rho_a, k_a;
    valtype rho_eq_til, K_eq_til;
    valtype R, T1, T2, T3, A_1, A_2;
    valtype N, A_hat, P_hat;
    valtype delta_1, delta_2, delta_3;

private:
    Analytical_MatEF& mat_F;
};


class xUsAirBiot
{
public:
    xUsAirBiot(double omega_, double beta_,  Analytical_MatEF& mat_F_);
    // using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta;
    valtype A, B, C, D;
    double A_inc, rho_a, k_a;
    valtype rho_eq_til, K_eq_til;
    valtype R, T1, T2, T3, A_1, A_2;
    valtype N, A_hat, P_hat;
    valtype delta_1, delta_2, delta_3;

private:
    Analytical_MatEF& mat_F;
};



