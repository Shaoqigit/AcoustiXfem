#pragma once
#include <iostream>
#include <vector>

#include "xEval.h"
#include "xGeomElem.h"
constexpr unsigned int max_mode = 35;

class xEvalCylinderCoefficient
{
public:
    xEvalCylinderCoefficient(double omega_, double radius_, double sigmaF_, double d_);
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;
    double A_inc, radius, rho_a, k_0;
    int inf_mode;
    valtype A[4][max_mode];
    valtype k_1, k_2, k_3;
    valtype mu_1, mu_2, mu_3;
    valtype N, A_hat;
    valtype rho_eq_til, K_eq_til;
    
};

class xEvalExactCylinderUs : public xfem::xEval<xtensor::xVector<std::complex<double> > > {

public:
    xEvalExactCylinderUs(double omega_, double radius_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;
    
private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, radius, sigmaF, d;
    int inf_mode;
    valtype A[4][max_mode];
    valtype k_1, k_2, k_3;
    
};


class xEvalExactCylinderTraction : public xfem::xEval<xtensor::xTensor2<std::complex<double>> > {

public:
    xEvalExactCylinderTraction(double omega_, double radius_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xTensor2<std::complex<double> >& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, radius, sigmaF, d;
    int inf_mode;
    valtype A[4][max_mode];
    valtype k_1, k_2, k_3;
    valtype N, A_hat ;

};


class xEvalExactCylinderUtotal : public xfem::xEval<xtensor::xVector<std::complex<double> > > {

public:
    xEvalExactCylinderUtotal(double omega_, double radius_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, sigmaF, d;
    double radius, k_0, rho_a, A_inc;
    int inf_mode;
    valtype A[4][max_mode];
    valtype mu_1, mu_2, mu_3;
    valtype k_1, k_2, k_3;
};


class xEvalExactCylinderPressure : public xfem::xEval<std::complex<double> > {

public:
    xEvalExactCylinderPressure(double omega_, double radius_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;
    double A_inc, rho_a, k_0;
    double radius;
    int inf_mode ;
    valtype A[4][max_mode];
    valtype mu_1, mu_2, mu_3;
    valtype k_1, k_2, k_3;
    valtype K_eq_til;

};


class xEvalExactCylinderVelocity : public xfem::xEval<xtensor::xVector<std::complex<double> > > {

public:
    xEvalExactCylinderVelocity(double omega_, double radius_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;
    double A_inc, rho_a, k_0;
    double radius;
    int inf_mode ;
    valtype A[4][max_mode];
    valtype mu_1, mu_2, mu_3;
    valtype k_1, k_2, k_3;
    valtype K_eq_til;
};


class xEvalExactCylinderVelocityS : public xfem::xEval<std::complex<double> > {

public:
    xEvalExactCylinderVelocityS(double omega_, double radius_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, sigmaF, d;
    double A_inc, rho_a, k_0;
    double radius;
    int inf_mode ;
    valtype A[4][max_mode];
    valtype mu_1, mu_2, mu_3;
    valtype k_1, k_2, k_3;
    valtype K_eq_til;
};

