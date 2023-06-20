#pragma once

#include "xEval.h"
#include "xGeomElem.h"
using namespace xfem;

class xEvalCoefficient
{
public:
    xEvalCoefficient(double omega_, double beta_, double sigmaF_, double d_);
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;

    valtype R, T1, T2, T3, A_1, A_2;
    valtype N, A_hat, P_hat;
    valtype delta_1, delta_2, delta_3;
    valtype mu_1, mu_2, mu_3;
    valtype rho_eq_til, K_eq_til;
    double k_a;
    
    
};

class xEvalExactBiotUs : public xEval<xtensor::xVector<std::complex<double> > > {

public:
    xEvalExactBiotUs(double omega_, double beta_, double sigmaF_, double d_);
    void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;
    
private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;
    valtype A_1, A_2, T3;
    valtype delta_1, delta_2, delta_3;
    double k_a;
    
};


class xEvalExactBiotTraction : public xEval<xtensor::xTensor2<std::complex<double>> > {

public:
    xEvalExactBiotTraction(double omega_, double beta_, double sigmaF_, double d_);
    void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<std::complex<double> >& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;
    valtype A_1, A_2, T3;
    valtype delta_1, delta_2, delta_3;
    valtype N, A_hat, P_hat;
    double k_a;
};


class xEvalExactBiotUtotal : public xEval<xtensor::xVector<std::complex<double> > > {

public:
    xEvalExactBiotUtotal(double omega_, double beta_, double sigmaF_, double d_);
    void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;
    valtype A_1, A_2, T3;
    valtype mu_1, mu_2, mu_3;
    valtype delta_1, delta_2, delta_3;
    double k_a;
};


class xEvalExactBiotPressure : public xEval<std::complex<double> > {

public:
    xEvalExactBiotPressure(double omega_, double beta_, double sigmaF_, double d_);
    void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double>& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;
    double k_a;
    valtype R, T1, T2;
    valtype delta_1, delta_2, delta_3;
    
};

class xEvalExactBiotVelocity : public xEval<xtensor::xVector<std::complex<double> > > {

public:
    xEvalExactBiotVelocity(double omega_, double beta_, double sigmaF_, double d_);
    void operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;
    double k_a;
    valtype R, T1, T2, T3;
    valtype delta_1, delta_2, delta_3;
};
