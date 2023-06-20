#include "exactSolutions.h"
#include "xMesh.h"
#include "complex_bessel.h"
#include "Eigen/Dense"
#include "AnalyticalMaterial.h"
#include <boost/math/special_functions/bessel.hpp> 

using namespace std;
using namespace xfem;
using namespace Eigen;
using namespace sp_bessel;
using namespace boost::math;


xEvalCylinderFluidCoefficient::xEvalCylinderFluidCoefficient(double omega_, double radius_, double sigmaF_, double d_) : omega(omega_), radius(radius_), sigmaF(sigmaF_), d(d_)
{
    double phi = 0.98, sigma = 45e3, alpha = 1.0, lambda_prime = 250e-6, lambda = 110e-6;
    rho_a = 1.213;
    Analytical_MatEF mat(omega, 0., phi, sigma, alpha, lambda_prime, lambda);
    k_0 = mat.k_a;
    k_eq = mat.k_eq_til;
    // k_eq = k_a;
    rho_eq = mat.rho_eq_til;
    // rho_eq = k_eq
    A_inc = 1.e-7;
    // double sigmaF = 775e3;
    // double d = 0.45e-3;
    // double a = 0.1;
    double a = radius;
    inf_mode = max_modes;
    for (int m=0; m<inf_mode; m++) {
        // double mm = 1.0*m;
        Matrix2cd M;
        Vector2cd b;

        // pressure 
        M(0,0) = rho_a*hankelH2(m,k_0*a);
        M(0,1) = -rho_eq*besselJ(m,k_eq*a)-sigmaF*d*k_eq/(j*omega)*besselJp(m,k_eq*a,1);
        // velocity
        M(1,0) = k_0*hankelH2p(m,k_0*a,1);
        M(1,1) = -k_eq*besselJp(m,k_eq*a,1);

        b(0) = -2.*A_inc*rho_a*pow((-j),m)*besselJ(m,k_0*a);
        b(1) = -2*A_inc*pow((-j),m)*k_0*besselJp(m,k_0*a,1);



        Matrix2cd M0 = M;
        Matrix2cd M1 = M;

        M0(0,0) = b(0);
        M0(1,0) = b(1);

        M1(0,1) = b(0);
        M1(1,1) = b(1);

        valtype det_M = M.determinant();
        valtype det_M0 = M0.determinant();
        valtype det_M1 = M1.determinant();
        // Vector2cd c;
        // c(0) = det_M0/det_M;
        // c(1) = det_M1/det_M;
        Vector2cd c = M.partialPivLu().solve(b);

        if (m==0){
            A[0][m] = 0.5*c(0);
            A[1][m] = 0.5*c(1);
        }
        else {
            A[0][m] = c(0);
            A[1][m] = c(1);

        }

    }
    

}




xEvalObliqueExactSolutionBouillard2DCplx::xEvalObliqueExactSolutionBouillard2DCplx(double omega_, double beta_, double sigmaF_, double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_) {
    rho_a = 1.213;
    double P = 101325., gamma = 1.4, mu = 0.1839e-4, Pr = 0.710;
    double K_a = gamma * P;
    valtype c_a = sqrt(K_a/rho_a);
    k_a = omega/c_a;
    double phi = 0.97, sigma = 57e3, alpha = 1.54, lambda_prime = 73.8e-6, lambda = 24.6e-6;
    Analytical_MatEF mat(omega, beta, phi, sigma, alpha, lambda_prime, lambda);
    // k_a = mat.k_a;
    k_eq = mat.k_eq_til;
    rho_eq = mat.rho_eq_til;
    // rho_eq = rho_a*10.;
    // valtype c_eq = sqrt(K_a/rho_eq);
    // k_eq = omega/c_eq;

    cout<<"wave number of EF: "<<k_eq<<endl;
    cout<<"wave number of air: "<<k_a<<endl;
    cout<<"EF density: "<<rho_eq<<endl;
}

void xEvalObliqueExactSolutionBouillard2DCplx::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double>& res) const {

    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    valtype k_ax = k_a*cos(beta);
    valtype k_ay = k_a*sin(beta);
    valtype k_eqx = sqrt(pow(k_eq,2)-pow(k_a,2)*pow(sin(beta),2));
    valtype k_eqy = k_ay;

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    //Check if side was prescribed
    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation < 0) cdgx = -1.e-8;
    else if(orientation > 0) cdgx = 1.e-8;
    valtype T = 2./(1.+rho_a*k_eqx/(rho_eq*k_ax)+sigmaF*d*k_eqx/(rho_eq*omega));
    valtype R = 1.-rho_a*k_eqx/(rho_eq*k_ax)*T;

    if (cdgx<0.){
        res = exp(-j*k_ax*x-j*k_ay*y)+R*exp(j*k_ax*x-j*k_ay*y);
    }
    else
    {
        res = T*exp(-j*k_eqx*x-j*k_eqy*y);
    }

    return;
}


xEvalObliqueSpeedBouillard2DCplx::xEvalObliqueSpeedBouillard2DCplx(double omega_, double beta_, double sigmaF_, double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_) {
    rho_a = 1.213;
    double P = 101325., gamma = 1.4, mu = 0.1839e-4, Pr = 0.710;
    double K_a = gamma * P;
    valtype c_a = sqrt(K_a/rho_a);
    k_a = omega/c_a;
    double phi = 0.97, sigma = 57e3, alpha = 1.54, lambda_prime = 73.8e-6, lambda = 24.6e-6;
    Analytical_MatEF mat(omega, beta, phi, sigma, alpha, lambda_prime, lambda);
    k_a = mat.k_a;
    k_eq = mat.k_eq_til;
    rho_eq = mat.rho_eq_til;
    // rho_eq = rho_a*10.;
    // valtype c_eq = sqrt(K_a/rho_eq);
    // k_eq = omega/c_eq;
}

void xEvalObliqueSpeedBouillard2DCplx::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const {

    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);
    valtype k_ax = k_a*cos(beta);
    valtype k_ay = k_a*sin(beta);
    valtype k_eqx = sqrt(pow(k_eq,2)-pow(k_a,2)*pow(sin(beta),2));
    valtype k_eqy = k_ay;

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    //Check if side was prescribed
    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation < 0) cdgx = -1.e-8;
    else if(orientation > 0) cdgx = 1.e-8;
    
    valtype T = 2./(1.+rho_a*k_eqx/(rho_eq*k_ax)+sigmaF*d*k_eqx/(rho_eq*omega));
    valtype R = 1.-rho_a*k_eqx/(rho_eq*k_ax)*T;

    if (cdgx<0.){
        res(0) = (-j*k_ax*exp(-j*k_ax*x-j*k_ay*y)+j*k_ax*R*exp(j*k_ax*x-j*k_ay*y));
        res(1) = (-j*k_ay*exp(-j*k_ax*x-j*k_ay*y)-j*k_ay*R*exp(j*k_ax*x-j*k_ay*y));
    }
    else
    {
        res(0) = (-j*k_eqx*T*exp(-j*k_eqx*x-j*k_eqy*y));
        res(1) = (-j*k_eqy*T*exp(-j*k_eqx*x-j*k_eqy*y));
    }
    
    return;
}

xEvalCylinderExactSolution::xEvalCylinderExactSolution(double omega_, double radius_, double sigmaF_, double d_) : omega(omega_), radius(radius_), sigmaF(sigmaF_), d(d_) {
    xEvalCylinderFluidCoefficient Coeff(omega, radius, sigmaF, d);
    inf_mode  = Coeff.inf_mode;
    for (int i=0; i<2; i++) {
        for (int j=0; j<inf_mode; j++) {
            A[i][j] = Coeff.A[i][j];
        }
    }
    radius  = Coeff.radius;
    inf_mode = Coeff.inf_mode;
    k_0 = Coeff.k_0;
    k_eq = Coeff.k_eq;
    rho_a = Coeff.rho_a;
    rho_eq =  Coeff.rho_eq;
    A_inc = Coeff.A_inc;
    cout<<"equivalent density"<<rho_eq<<endl;
    cout<<"equivalent wave number"<<k_eq<<endl;
    // cout<<"reflection coefficient"<<A[0][10]<<endl;
}

void xEvalCylinderExactSolution::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double>& res) const{
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    double r = sqrt(x * x + y * y);
    double th = atan2(y, x);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    double cdgy = cdgxyz(1);

    double cdgr = sqrt(cdgx * cdgx + cdgy * cdgy);

    // 2D computation numerically
    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);
    if(orientation < 0) cdgr = radius-1.e-8;
    else if(orientation > 0) cdgr = radius+1.e-8;

    if (cdgr < radius){
        res = 0.;
        // valtype res1 = 0.;
        for (int m=0; m<inf_mode; m++)
        {
            double mm = 1.0*m;
          
                res += rho_eq*pow(omega,2)*A[1][m]*besselJ(m,k_eq*r)*cos(mm*th);
            
        }
    }
    else {
        res = 0.;
        valtype res1 = 0.;
        valtype res2 = 0.;
        for (int m=0; m<inf_mode; m++)
        {
            double mm = 1.0*m;
            if (m==0){
                res1 += 0.; 
            }
            else {
                res1 += pow((-j),m)*besselJ(m,k_0*r)*cos(mm*th);
            }
        res2 += A[0][m]*hankelH2(m, k_0*r)*cos(mm*th);
        }
        res = rho_a*pow(omega,2)*A_inc*besselJ(0,k_0*r) + A_inc*rho_a*pow(omega,2)*2.*res1+ rho_a*pow(omega,2)*res2;
    }
    return;
    
}

xEvalCylinderExactspeed::xEvalCylinderExactspeed(double omega_, double radius_, double sigmaF_, double d_) : omega(omega_), radius(radius_), sigmaF(sigmaF_), d(d_) {
    xEvalCylinderFluidCoefficient Coeff(omega, radius, sigmaF, d);
    inf_mode  = Coeff.inf_mode;
    for (int i=0; i<2; i++) {
        for (int j=0; j<inf_mode; j++) {
            A[i][j] = Coeff.A[i][j];
        }
    }
    radius  = Coeff.radius;
    inf_mode = Coeff.inf_mode;
    k_0 = Coeff.k_0;
    k_eq = Coeff.k_eq;
    rho_a = Coeff.rho_a;
    rho_eq =  Coeff.rho_eq;
    A_inc = Coeff.A_inc;

}

void xEvalCylinderExactspeed::operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const {
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    double r = sqrt(x * x + y * y);
    double th = atan2(y, x);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    double cdgy = cdgxyz(1);

    double cdgr = sqrt(cdgx * cdgx + cdgy * cdgy);

    // 2D computation numerically
    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);
    if(orientation < 0) cdgr = radius-1.e-8;
    else if(orientation > 0) cdgr = radius+1.e-8;

    if (cdgr < radius){
        res(0) = 0.;
        res(1) = 0.;
        // res(1) = 0.;
        // valtype res1 = 0.;
        // for (int m=0; m<inf_mode; m++)
        // {
        //     double mm = 1.0*m;
        //         res += k_eq*rho_eq*pow(omega,2)*A[1][m]*besselJp(0, k_eq*r,1)*cos(mm*th);
            
        // }
    }
    else {
        // res(0) = 0.;
        // res(1) = 0.;
        
        valtype res11 = 0.;
        valtype res12 = 0.;
        valtype res21 = 0.;
        valtype res22 = 0.;
        valtype vrr = 0.;
        valtype vrt = 0.;
        for (int m=0; m<inf_mode; m++)
        {
            double mm = 1.0*m;
            if (m==0){
                // res1 += rho_a*pow(omega,2)*A_inc*besselJ(0,k_0*r); 
                res11 += k_0*rho_a*pow(omega,2)*A_inc*besselJp(0,k_0*r,1); 
                res21 += 0.;
            }
            else {
                // res1 += A_inc*rho_a*pow(omega,2)*2.*pow((-j),m)*besselJ(m,k_0*r)*cos(mm*th);
                res11 += k_0*A_inc*rho_a*pow(omega,2)*2.*pow((-j),m)*besselJp(m,k_0*r,1)*cos(mm*th);
                res21 += -mm/r*A_inc*rho_a*pow(omega,2)*2.*pow((-j),m)*besselJ(m,k_0*r)*sin(mm*th);
            }
        // res12 += rho_a*pow(omega,2)*A[0][m]*hankelH2(m, k_0*r)*cos(mm*th);
        res12 += k_0*rho_a*pow(omega,2)*A[0][m]*hankelH2p(m, k_0*r,1)*cos(mm*th);
        res22 += -mm/r*rho_a*pow(omega,2)*A[0][m]*hankelH2(m, k_0*r)*sin(mm*th);
        }
        vrr += res11 + res12;
        vrt += res21 + res22;
        res(0) = vrr*cos(th)-vrt*sin(th);
        res(1) = vrr*sin(th)+vrt*cos(th);
    }
    return;
}


xEvalCylinderExactspeedS::xEvalCylinderExactspeedS(double omega_, double radius_, double sigmaF_, double d_) : omega(omega_), radius(radius_), sigmaF(sigmaF_), d(d_) {
    xEvalCylinderFluidCoefficient Coeff(omega, radius, sigmaF, d);
    inf_mode  = Coeff.inf_mode;
    for (int i=0; i<2; i++) {
        for (int j=0; j<inf_mode; j++) {
            A[i][j] = Coeff.A[i][j];
        }
    }
    radius  = Coeff.radius;
    inf_mode = Coeff.inf_mode;
    k_0 = Coeff.k_0;
    k_eq = Coeff.k_eq;
    rho_a = Coeff.rho_a;
    rho_eq =  Coeff.rho_eq;
    A_inc = Coeff.A_inc;

}

void xEvalCylinderExactspeedS::operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, std::complex<double>& res) const {
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    double r = sqrt(x * x + y * y);
    double th = atan2(y, x);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    double cdgy = cdgxyz(1);

    double cdgr = sqrt(cdgx * cdgx + cdgy * cdgy);

    // 2D computation numerically
    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);
    if(orientation < 0) cdgr = radius-1.e-8;
    else if(orientation > 0) cdgr = radius+1.e-8;

    if (cdgr < radius){
        res = 0.;
        // res(1) = 0.;
        // valtype res1 = 0.;
        // for (int m=0; m<inf_mode; m++)
        // {
        //     double mm = 1.0*m;
        //         res += k_eq*rho_eq*pow(omega,2)*A[1][m]*besselJp(0, k_eq*r,1)*cos(mm*th);
            
        // }
    }
    else {
        // res(0) = 0.;
        // res(1) = 0.;
        
        valtype res11 = 0.;
        valtype res12 = 0.;
        valtype res21 = 0.;
        valtype res22 = 0.;
        valtype vrr = 0.;
        valtype vrt = 0.;
        for (int m=0; m<inf_mode; m++)
        {
            double mm = 1.0*m;
            if (m==0){
                // res1 += rho_a*pow(omega,2)*A_inc*besselJ(0,k_0*r); 
                res11 += 0.;
                res21 += 0.;
            }
            else {
                // res1 += A_inc*rho_a*pow(omega,2)*2.*pow((-j),m)*besselJ(m,k_0*r)*cos(mm*th);
                res11 += pow((-j),mm)*k_0*besselJp(m,k_0*r,1)*cos(mm*th);
                res21 += -mm/r*A_inc*rho_a*pow(omega,2)*2.*pow((-j),m)*besselJ(m,k_0*r)*sin(mm*th);
            }
        // res12 += rho_a*pow(omega,2)*A[0][m]*hankelH2(m, k_0*r)*cos(mm*th);
        res12 += k_0*A[0][m]*hankelH2p(m, k_0*r,1)*cos(mm*th);
        res22 += -mm/r*rho_a*pow(omega,2)*A[0][m]*hankelH2(m, k_0*r)*sin(mm*th);
        }
        vrr += res11 + res12;
        vrt += res21 + res22;
        res = k_0*rho_a*pow(omega,2)*A_inc*besselJp(0,k_0*r,1)+A_inc*rho_a*pow(omega,2)*2.*res11+rho_a*pow(omega,2)*res12;
        // res(0) = vrr*cos(th)-vrt*sin(th);
        // res(1) = vrr*sin(th)+vrt*cos(th);
    }
    return;
}

xEvalExactSolutionBouillard2DCplxTwo::xEvalExactSolutionBouillard2DCplxTwo(double omega_, double beta_, double d_, std::complex<double> v0_) : omega(omega_), beta(beta_), d(d_), v0(v0_) {}

void xEvalExactSolutionBouillard2DCplxTwo::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double>& res) const{

    const double l_a=1., l_eq=1.;

    const double P = 101325., gamma = 1.4, rho_a = 1.213, mu = 0.1839e-4, Pr = 0.710, phi = 0.994, sigma = 9045, alpha = 1.02, lambda_prime = 197e-6, lambda = 103e-6;
    double K_a = gamma * P;
    double c_a = sqrt(K_a/rho_a);
    double k_a = omega/c_a;
    double Z_a = rho_a * c_a;
    double nu = mu / rho_a;
    double nu_prime = nu / Pr;

    //// fluid equivalent part######################
    valtype omega_0 = sigma * phi / (rho_a * alpha);
    valtype omega_inf = pow(sigma, 2)*pow(phi,2)*pow(lambda,2)/(4.*mu*rho_a*alpha);
//    valtype rho_eq  = (rho_a*alpha/phi)*(1.+((omega_0)/(j*omega))*sqrt(1.+j*omega/omega_inf));
    valtype rho_eq = rho_a;
    valtype omega_prime_infty = (16.*nu_prime)/(pow(lambda_prime, 2));
    valtype F_prime_CA = sqrt(1.+j*omega/omega_prime_infty);
    valtype alpha_prime_til = 1.+omega_prime_infty*F_prime_CA/(2.*j*omega);
    valtype K = (gamma*P/phi)/(gamma-(gamma-1.)/alpha_prime_til);
    valtype c_eq = sqrt(K/rho_eq);
//    valtype k_eq = omega/c_eq;
    valtype k_eq = k_a;

    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);
    valtype k_ax = k_a*cos(beta);
    valtype k_ay = k_a*sin(beta);
    valtype k_eqx = k_eq*cos(beta);
    valtype k_eqy = k_eq*sin(beta);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    valtype C = -pow(k_ax,2)*k_eq*sin(k_ax*l_a)*sin(k_eqx*l_eq)*sigma*d*exp(-j*k_eqy*y)+j*omega*pow(k_ax,2)*rho_eq*sin(k_ax*l_a)*cos(k_eqx*l_eq)*exp(-j*k_ay*y)+j*omega*k_ax*k_eqx*rho_a*sin(k_eqx*l_eq)*cos(k_ax*l_a)*exp(-j*k_eqy*y);


    if (cdgx<0.){
        res = (exp(-j*k_ay*y)*exp(-j*k_eqy*y)*j*v0*sigma*d*omega*k_ax*k_eqx*rho_a*sin(k_eqx*l_eq)+exp(-j*k_ay*y)*exp(-j*k_eqy*y)*v0*pow(omega,2)*rho_a*rho_eq*k_ax*cos(k_eqx*l_eq))*cos(k_ax*x)/C+exp(-j*k_ay*y)*exp(-j*k_eqy*y)*v0*pow(omega,2)*pow(rho_a, 2)*k_eqx*sin(k_eqx*l_eq)*sin(k_ax*x)/C;
    }
    else
    {
        res = exp(-j*k_eqy*y)*exp(-j*k_eqy*y)*v0*pow(omega,2)*rho_a*rho_eq*k_ax*cos(k_eqx*l_eq)*cos(k_eqx*x)/C+exp(-j*k_eqy*y)*exp(-j*k_eqy*y)*v0*pow(omega,2)*rho_a*rho_eq*k_ax*sin(k_eqx*l_eq)*sin(k_eqx*x)/C;
    }
    return;
}


xEvalGradExactSolutionBouillard2DCplx::xEvalGradExactSolutionBouillard2DCplx(double omega_, double beta_, double d_, std::complex<double> v0_) : omega(omega_), beta(beta_), d(d_), v0(v0_) {
    if(beta > 2 * M_PI){
        cout<<"xEvalGradExactSolutionBouillard2DCplx has to be updated to accept beta in degree !!!\n";
        throw;
    }
}

void xEvalGradExactSolutionBouillard2DCplx::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> > & res) const {
    const double l_a=1., l_eq=1.;

    const double P = 101325., gamma = 1.4, rho_a = 1.213, mu = 0.1839e-4, Pr = 0.710, phi = 0.994, sigma = 9045, alpha = 1.02, lambda_prime = 197e-6, lambda = 103e-6;
    double K_a = gamma * P;
    double c_a = sqrt(K_a/rho_a);
    double k_a = omega/c_a;
    double Z_a = rho_a * c_a;
    double nu = mu / rho_a;
    double nu_prime = nu / Pr;

    //// fluid equivalent part######################
    valtype omega_0 = sigma * phi / (rho_a * alpha);
    valtype omega_inf = pow(sigma, 2)*pow(phi,2)*pow(lambda,2)/(4.*mu*rho_a*alpha);
//    valtype rho_eq  = (rho_a*alpha/phi)*(1.+((omega_0)/(j*omega))*sqrt(1.+j*omega/omega_inf));
    valtype rho_eq = rho_a;
    valtype omega_prime_infty = (16.*nu_prime)/(pow(lambda_prime, 2));
    valtype F_prime_CA = sqrt(1.+j*omega/omega_prime_infty);
    valtype alpha_prime_til = 1.+omega_prime_infty*F_prime_CA/(2.*j*omega);
    valtype K = (gamma*P/phi)/(gamma-(gamma-1.)/alpha_prime_til);
    valtype c_eq = sqrt(K/rho_eq);
//    valtype k_eq = omega/c_eq;
    valtype k_eq = k_a;
    valtype C = -pow(k_a,2)*k_eq*sin(k_a*l_a)*sin(k_eq*l_eq)*sigma*d+j*omega*pow(k_a,2)*rho_eq*sin(k_a*l_a)*cos(k_eq*l_eq)+j*omega*k_a*k_eq*rho_a*sin(k_eq*l_eq)*cos(k_a*l_a);
    xtensor::xPoint xyz = geo_appro->getXYZ();
    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    double x = xyz(0);
    double y = xyz(1);
    double b = beta;
    //     double rc=1.21*340;

//    res(0) = {-k*sin(k*(x*cos(b) + y*sin(b)))*cos(b),  k*cos(b)*cos(k*(x*cos(b) + y*sin(b)))};
//    res(1) = {-k*sin(b)*sin(k*(x*cos(b) + y*sin(b))), k*sin(b)*cos(k*(x*cos(b) + y*sin(b)))};

    if (cdgx<0.){
        res(0) = (-k_a*cos(b)*(j*v0*sigma*d*omega*k_a*k_eq*rho_a*sin(k_eq*l_eq)+v0*pow(omega,2)*rho_a*rho_eq*k_a*cos(k_eq*l_eq))*sin(k_a*(x*cos(b)+y*sin(b)))/C+k_a*cos(b)*v0*pow(omega,2)*pow(rho_a, 2)*k_eq*sin(k_eq*l_eq)*cos(k_a*(x*cos(b)+y*sin(b)))/C)/(-j*rho_a*omega);
        res(1) = (-k_a*sin(b)*(j*v0*sigma*d*omega*k_a*k_eq*rho_a*sin(k_eq*l_eq)+v0*pow(omega,2)*rho_a*rho_eq*k_a*cos(k_eq*l_eq))*sin(k_a*(x*cos(b)+y*sin(b)))/C+k_a*sin(b)*v0*pow(omega,2)*pow(rho_a, 2)*k_eq*sin(k_eq*l_eq)*cos(k_a*(x*cos(b)+y*sin(b)))/C)/(-j*rho_a*omega);
    }
    else
    {
        res(0) = (-k_eq*cos(b)*v0*pow(omega,2)*rho_a*rho_eq*k_a*cos(k_eq*l_eq)*sin(k_eq*(x*cos(b)+y*sin(b)))/C+k_eq*cos(b)*v0*pow(omega,2)*rho_a*rho_eq*k_a*sin(k_eq*l_eq)*cos(k_eq*(x*cos(b)+y*sin(b)))/C)/(-j*rho_eq*omega);
        res(1) = (-k_eq*sin(b)*v0*pow(omega,2)*rho_a*rho_eq*k_a*cos(k_eq*l_eq)*sin(k_eq*(x*cos(b)+y*sin(b)))/C+k_eq*sin(b)*v0*pow(omega,2)*rho_a*rho_eq*k_a*sin(k_eq*l_eq)*cos(k_eq*(x*cos(b)+y*sin(b)))/C)/(-j*rho_eq*omega);
    }

    return;
}


xEvalExactSpeedBouillard2DCplx::xEvalExactSpeedBouillard2DCplx(double omega_, double beta_, double d_, std::complex<double> v0_) : omega(omega_), beta(beta_), d(d_), v0(v0_) {
    if(beta > 2 * M_PI){
        cout<<"xEvalGradExactSolutionBouillard2DCplx has to be updated to accept beta in degree !!!\n";
        throw;
    }
}

void xEvalExactSpeedBouillard2DCplx::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> > & res) const {
    const double l_a=1., l_eq=1.;

    const double P = 101325., gamma = 1.4, rho_a = 1.213, mu = 0.1839e-4, Pr = 0.710, phi = 0.994, sigma = 9045, alpha = 1.02, lambda_prime = 197e-6, lambda = 103e-6;
    double K_a = gamma * P;
    double c_a = sqrt(K_a/rho_a);
    double k_a = omega/c_a;
    double Z_a = rho_a * c_a;
    double nu = mu / rho_a;
    double nu_prime = nu / Pr;

    //// fluid equivalent part######################
    valtype omega_0 = sigma * phi / (rho_a * alpha);
    valtype omega_inf = pow(sigma, 2)*pow(phi,2)*pow(lambda,2)/(4.*mu*rho_a*alpha);
//    valtype rho_eq  = (rho_a*alpha/phi)*(1.+((omega_0)/(j*omega))*sqrt(1.+j*omega/omega_inf));
    valtype rho_eq = rho_a;
    valtype omega_prime_infty = (16.*nu_prime)/(pow(lambda_prime, 2));
    valtype F_prime_CA = sqrt(1.+j*omega/omega_prime_infty);
    valtype alpha_prime_til = 1.+omega_prime_infty*F_prime_CA/(2.*j*omega);
    valtype K = (gamma*P/phi)/(gamma-(gamma-1.)/alpha_prime_til);
    valtype c_eq = sqrt(K/rho_eq);
//    valtype k_eq = omega/c_eq;
    valtype k_eq = k_a;


    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);
    valtype k_ax = k_a*cos(beta);
    valtype k_ay = k_a*sin(beta);
    valtype k_eqx = k_eq*cos(beta);
    valtype k_eqy = k_eq*sin(beta);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    valtype C = -pow(k_ax,2)*k_eq*sin(k_ax*l_a)*sin(k_eqx*l_eq)*sigma*d*exp(-j*k_eqy*y)+j*omega*pow(k_ax,2)*rho_eq*sin(k_ax*l_a)*cos(k_eqx*l_eq)*exp(-j*k_ay*y)+j*omega*k_ax*k_eqx*rho_a*sin(k_eqx*l_eq)*cos(k_ax*l_a)*exp(-j*k_eqy*y);
    //     double rc=1.21*340;

//    res(0) = {-k*sin(k*(x*cos(b) + y*sin(b)))*cos(b),  k*cos(b)*cos(k*(x*cos(b) + y*sin(b)))};
//    res(1) = {-k*sin(b)*sin(k*(x*cos(b) + y*sin(b))), k*sin(b)*cos(k*(x*cos(b) + y*sin(b)))};

    if (cdgx<0.){
        res(0) = (-k_ax*(exp(-j*k_ay*y)*exp(-j*k_eqy*y)*j*v0*sigma*d*omega*k_ax*k_eqx*rho_a*sin(k_eqx*l_eq)+exp(-j*k_ay*y)*exp(-j*k_eqy*y)*v0*pow(omega,2)*rho_a*rho_eq*k_ax*cos(k_eqx*l_eq))*sin(k_ax*x)/C+k_ax*exp(-j*k_ay*y)*exp(-j*k_eqy*y)*v0*pow(omega,2)*pow(rho_a, 2)*k_eqx*sin(k_eqx*l_eq)*cos(k_ax*x)/C)/(-j*rho_a*omega);
        res(1) = ((-j*k_ay-j*k_eqy)*(exp(-j*k_ay*y)*exp(-j*k_eqy*y)*j*v0*sigma*d*omega*k_ax*k_eqx*rho_a*sin(k_eqx*l_eq)+(-j*k_ay-j*k_eqy)*exp(-j*k_ay*y)*exp(-j*k_eqy*y)*v0*pow(omega,2)*rho_a*rho_eq*k_ax*cos(k_eqx*l_eq))*cos(k_ax*x)/C+(-j*k_ay-j*k_eqy)*exp(-j*k_ay*y)*exp(-j*k_eqy*y)*v0*pow(omega,2)*pow(rho_a, 2)*k_eqx*sin(k_eqx*l_eq)*sin(k_ax*x)/C)/(-j*rho_a*omega);
    }
    else
    {
        res(0) = (-k_eqx*exp(-j*k_eqy*y)*exp(-j*k_eqy*y)*v0*pow(omega,2)*rho_a*rho_eq*k_ax*cos(k_eqx*l_eq)*sin(k_eqx*x)/C+k_eqx*exp(-j*k_eqy*y)*exp(-j*k_eqy*y)*v0*pow(omega,2)*rho_a*rho_eq*k_ax*sin(k_eqx*l_eq)*cos(k_eqx*x)/C)/(-j*rho_eq*omega);
        res(1) = ((-j*k_eqy-j*k_eqy)*exp(-j*k_eqy*y)*exp(-j*k_eqy*y)*v0*pow(omega,2)*rho_a*rho_eq*k_ax*cos(k_eqx*l_eq)*cos(k_eqx*x)/C+(-j*k_eqy-j*k_eqy)*exp(-j*k_eqy*y)*exp(-j*k_eqy*y)*v0*pow(omega,2)*rho_a*rho_eq*k_ax*sin(k_eqx*l_eq)*sin(k_eqx*x)/C)/(-j*rho_eq*omega);
    }

    return;
}


xEvalPointSourceCplx::xEvalPointSourceCplx(double k_, xtensor::xPoint s_) :
k(k_), O(0.,0.,0.), s(s_) {}

void xEvalPointSourceCplx::operator()(const xGeomElem*  geo_appro, const
xGeomElem* geo_integ, std::complex<double>& res) const {

     xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);
     xtensor::xVector<>  SX(s,xyz);


     res = 0.05*j * hankelH2(0., k * ((x-s(0))*(x-s(0))+(y-s(1))*(y-s(1))));
    // res = exp(-2e-5*k*((x-s(0))*(x-s(0))+(y-s(1))*(y-s(1))));

     return;
}
