#pragma once

#include <iostream>

#include "xEval.h"
#include "xGeomElem.h"

constexpr unsigned int max_modes = 30;


class xEvalCylinderFluidCoefficient
{
public:
    xEvalCylinderFluidCoefficient(double omega_, double radius_, double sigmaF_, double d_);
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, sigmaF, d;
    double A_inc, radius, rho_a, k_0;
    int inf_mode;
    valtype A[2][max_modes];
    valtype rho_eq, k_eq;
    
};

class xEvalObliqueExactSolutionBouillard2DCplx : public xfem::xEval<std::complex<double> > {

public:
    xEvalObliqueExactSolutionBouillard2DCplx(double omega_, double beta_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;
    valtype rho_a, rho_eq, k_a, k_eq;

};


class xEvalObliqueSpeedBouillard2DCplx : public xfem::xEval<xtensor::xVector<std::complex<double> > > {

public:
    xEvalObliqueSpeedBouillard2DCplx(double omega_, double beta_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, sigmaF, d;
    valtype rho_a, rho_eq, k_a, k_eq;

};

class xEvalCylinderExactSolution : public xfem::xEval<std::complex<double> > {

public:
    xEvalCylinderExactSolution(double omega_, double radius_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, radius, k_0, rho_a, A_inc, sigmaF, d;
    int inf_mode;
    valtype A[2][max_modes];
    valtype rho_eq, k_eq;

};


class xEvalCylinderExactspeed : public xfem::xEval<xtensor::xVector<std::complex<double> > > {

public:
    xEvalCylinderExactspeed(double omega_, double radius_, double sigmaF, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, sigmaF, d;
    double radius, k_0, rho_a, A_inc;
    int inf_mode;
    valtype A[2][max_modes];
    valtype rho_eq, k_eq;

};

class xEvalCylinderExactspeedS : public xfem::xEval<std::complex<double> > {

public:
    xEvalCylinderExactspeedS(double omega_, double radius_, double sigmaF_, double d_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, sigmaF, d;
    double radius, k_0, rho_a, A_inc;
    int inf_mode;
    valtype A[2][max_modes];
    valtype rho_eq, k_eq;

};


class xEvalExactSolutionBouillard2DCplxTwo : public xfem::xEval<std::complex<double> > {

public:
    xEvalExactSolutionBouillard2DCplxTwo(double omega_, double beta_, double d_, std::complex<double> v0_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const;

private:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};
    double omega, beta, d;
    valtype v0;

};


class xEvalGradExactSolutionBouillard2DCplx : public xfem::xEval<xtensor::xVector<std::complex<double> > > {

public :
    xEvalGradExactSolutionBouillard2DCplx(double omega_, double beta_, double d_, std::complex<double> v0_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> > & res) const ;

private:
    using valtype = std::complex<double>;
    const std::complex<double> j{0.,1.};
    double omega, beta, d;
    valtype v0;

};


class xEvalExactSpeedBouillard2DCplx : public xfem::xEval<xtensor::xVector<std::complex<double> > > {

public :
    xEvalExactSpeedBouillard2DCplx(double omega_, double beta_, double d_, std::complex<double> v0_);
    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> > & res) const ;

private:
    using valtype = std::complex<double>;
    const std::complex<double> j{0.,1.};
    double omega, beta, d;
    valtype v0;

};

class xEvalPointSourceCplx : public xfem::xEval<std::complex<double> > {

public :
     xEvalPointSourceCplx(double k_, xtensor::xPoint s_);
     void operator()(const  xfem::xGeomElem*  geo_appro, const xfem::xGeomElem*
geo_integ, std::complex<double>& res) const ;

private:
     double k;
     const std::complex<double> j{0.,1.};
     xtensor::xVector<>  betaDir;
     xtensor::xPoint O,s;
};
