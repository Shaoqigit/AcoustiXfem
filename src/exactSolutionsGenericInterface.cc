#include "exactSolutionsGenericInterface.h"
#include "xMesh.h"
#include "Eigen/Dense"
#include "AnalyticalMaterial.h"
#include <boost/math/special_functions/bessel.hpp> 

using namespace std;
using namespace xfem;
using namespace Eigen;
// using namespace sp_bessel;
using namespace boost::math;


xEvalInterfaceCoefficient::xEvalInterfaceCoefficient(double omega_, double beta_, valtype A_, valtype B_, valtype C_, valtype D_) : omega(omega_), beta(beta_), A(A_), B(B_), C(C_), D(D_)
{
    double phi = 0.97, sigma = 57e3, alpha = 1.54, lambda_prime = 73.8e-6, lambda = 24.6e-6;
    rho_a = 1.213;
    // doube beta = 45.;
    Analytical_MatEF mat(omega, beta, phi, sigma, alpha, lambda_prime, lambda);
    rho_eq = mat.rho_eq_til;
    k_eq = mat.k_eq_til;
    k_a = mat.k_a;
    // rho_eq = rho_a;
    // k_eq = k_a;


    valtype k_ax = k_a*cos(beta);
    valtype k_ay = k_a*sin(beta);
    valtype k_eqx = sqrt(pow(k_eq,2)-pow(k_a,2)*pow(sin(beta),2));
    valtype k_eqy = k_ay;

    Matrix2cd M;
    Vector2cd b;

        // pressure 
    // M(0,0) = B*k_ax/(omega*rho_a) - A;
    // M(0,1) = 1.;
    //     // velocity
    // M(1,0) = D*k_ax/(omega*rho_a) - C;
    // M(1,1) = k_eqx/(omega*rho_eq);

    // b(0) = A + B*k_ax/(omega*rho_a);
    // b(1) = C + D*k_ax/(omega*rho_a);

    // Vector2cd c = M.partialPivLu().solve(b);

    // R = c(0);
    // T = c(1);

    M(0,0) = 1;
    M(0,1) = -A-B*k_eqx/(omega*rho_eq);
        // velocity
    M(1,0) = -k_ax/(omega*rho_a);
    M(1,1) = -C - D*k_eqx/(omega*rho_eq);

    b(0) = -1.;
    b(1) = -k_ax/(omega*rho_a);

    Vector2cd c = M.partialPivLu().solve(b);

    R = c(0);
    T = c(1);

    cout<<"reflection coeff:"<<R<<endl;
    cout<<"Transmission coeff:"<<T<<endl;

    

}




xSolObliquePlaneGenericInterface::xSolObliquePlaneGenericInterface(double omega_, double beta_, valtype A_, valtype B_, valtype C_, valtype D_) : omega(omega_), beta(beta_), A(A_), B(B_), C(C_), D(D_) {
    rho_a = 1.213;

    xEvalInterfaceCoefficient coeff(omega, beta, A, B, C, D);
    k_a = coeff.k_a;
    k_eq = coeff.k_eq;
    // k_eq = k_a;
    rho_eq = coeff.rho_eq;
    // rho_eq = rho_a;
    R = coeff.R;
    T = coeff.T;

    cout<<"wave number of EF: "<<k_eq<<endl;
    cout<<"wave number of air: "<<k_a<<endl;
}

void xSolObliquePlaneGenericInterface::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double>& res) const {

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


    if (cdgx<0.){
        res = exp(-j*k_ax*x-j*k_ay*y)+R*exp(j*k_ax*x-j*k_ay*y);
    }
    else
    {
        res = T*exp(-j*k_eqx*x-j*k_eqy*y);
    }

    return;
}


xGradSolObliquePlaneGenericInterface::xGradSolObliquePlaneGenericInterface(double omega_, double beta_, valtype A_, valtype B_, valtype C_, valtype D_) : omega(omega_), beta(beta_), A(A_), B(B_), C(C_), D(D_) {
    rho_a = 1.213;

    xEvalInterfaceCoefficient coeff(omega, beta, A, B, C, D);
    k_a = coeff.k_a;
    k_eq = coeff.k_eq;
    // k_eq = k_a;
    rho_eq = coeff.rho_eq;
    // rho_eq = rho_a;
    R = coeff.R;
    T = coeff.T;
    // rho_eq = rho_a;
}

void xGradSolObliquePlaneGenericInterface::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const {


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


xEvalTubeCoeff::xEvalTubeCoeff(double omega_,  double v0_, valtype M11_, valtype M12_, valtype M21_, valtype M22_) : omega(omega_), v0(v0_), M11(M11_), M12(M12_), M21(M21_), M22(M22_)
{
    double phi = 0.98, sigma = 45e3, alpha = 1.0, lambda_prime = 250e-6, lambda = 110e-6;
    rho_a = 1.213;
    double beta = 0.;
    Analytical_MatEF mat(omega, beta, phi, sigma, alpha, lambda_prime, lambda);
    rho_eq = mat.rho_eq_til;
    k_eq = mat.k_eq_til;
    k_a = mat.k_a;
    // rho_eq = rho_a;
    // k_eq = k_a;

    double l1 = 0.099;
    double l2 = 0.1;

    
    valtype k_ax = k_a*cos(beta);
    valtype k_ay = k_a*sin(beta);
    valtype k_eqx = sqrt(pow(k_eq,2)-pow(k_a,2)*pow(sin(beta),2));
    valtype k_eqy = k_ay;


    Matrix4cd M;
    Vector4cd b;
// normal velocity is zero
    M(0,0) = -j*k_ax*exp(-j*k_ax*(-l1));
    M(0,1) = j*k_ax*exp(j*k_ax*(-l1));
    M(0,2) = 0.;
    M(0,3) = 0;
// wal BC
    M(1,0) = 0.;
    M(1,1) = 0.;
    M(1,2) = -j*k_eqx*exp(-j*k_eqx*(l2));
    M(1,3) = j*k_eqx*exp(j*k_eqx*(l2));
// interface p
    M(2,0) = 1.;
    M(2,1) = 1.;
    M(2,2) = -M11+j*k_eqx*M12/(-j*omega*rho_eq);
    M(2,3) = -M11-j*k_eqx*M12/(-j*omega*rho_eq);
// interface p
    M(3,0) = -j*k_ax/(-j*omega*rho_a);
    M(3,1) = j*k_ax/(-j*omega*rho_a);
    M(3,2) = -M21+j*k_eqx*M22/(-j*omega*rho_eq);
    M(3,3) = -M21-j*k_eqx*M22/(-j*omega*rho_eq);

    cout<<"Coeff 1: "<<exp(-j*k_ax*(-l1))<<endl;
    cout<<"Coeff 2: "<<-M21+j*k_eqx*M22<<endl;


    b(0) = v0;
    b(1) = 0.;
    b(2) = 0.;
    b(3) = 0.;

    Vector4cd c = M.partialPivLu().solve(b);

    A = c(0);
    B = c(1);
    C = c(2);
    D = c(3);
    

}



xPressureTubeNormal::xPressureTubeNormal(double omega_,  double v0_, valtype M11_, valtype M12_, valtype M21_, valtype M22_) : omega(omega_), v0(v0_), M11(M11_), M12(M12_), M21(M21_), M22(M22_){
    rho_a = 1.213;

    xEvalTubeCoeff coeff(omega, v0, M11, M12, M21, M22);
    k_a = coeff.k_a;
    k_eq = coeff.k_eq;
    // k_eq = k_a;
    rho_eq = coeff.rho_eq;
    // rho_eq = rho_a;
    A = coeff.A;
    B = coeff.B;
    C = coeff.C;
    D = coeff.D;

    cout<<"Coeff A: "<<A<<endl;
    cout<<"Coeff B: "<<B<<endl;
    cout<<"Coeff C: "<<C<<endl;
    cout<<"Coeff D: "<<D<<endl;
}

void xPressureTubeNormal::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double>& res) const {

    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);
    double beta = 0.;
    valtype k_ax = k_a*cos(beta);
    valtype k_ay = k_a*sin(beta);
    valtype k_eqx = sqrt(pow(k_eq,2)-pow(k_a,2)*pow(sin(beta),2));
    valtype k_eqy = k_ay;
    // cout<<"ky: "<<k_eqy<<endl;

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    //Check if side was prescribed
    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation < 0) cdgx = -1.e-8;
    else if(orientation > 0) cdgx = 1.e-8;


    if (cdgx<0.){
        res = A*exp(-j*k_ax*x-j*k_ay*y)+B*exp(j*k_ax*x-j*k_ay*y);
    }
    else
    {
        res = C*exp(-j*k_eqx*x-j*k_eqy*y)+D*exp(j*k_eqx*x-j*k_eqy*y);
    }

    return;
}


xEvalCoeffObliqueRigidWall::xEvalCoeffObliqueRigidWall(double omega_,  double beta_, valtype M11_, valtype M12_, valtype M21_, valtype M22_) : omega(omega_), beta(beta_), M11(M11_), M12(M12_), M21(M21_), M22(M22_)
{
    double phi = 0.98, sigma = 45e3, alpha = 1.0, lambda_prime = 250e-6, lambda = 110e-6;
    rho_a = 1.213;
    double beta = 0.;
    Analytical_MatEF mat(omega, beta, phi, sigma, alpha, lambda_prime, lambda);
    rho_eq = mat.rho_eq_til;
    k_eq = mat.k_eq_til;
    k_a = mat.k_a;
    // rho_eq = rho_a;
    // k_eq = k_a;
    double l = 0.1;

    
    valtype k_ax = k_a*cos(beta);
    valtype k_ay = k_a*sin(beta);
    valtype k_eqx = sqrt(pow(k_eq,2)-pow(k_a,2)*pow(sin(beta),2));
    valtype k_eqy = k_ay;


    Matrix3cd M;
    Vector3cd b;

// wal BC
    M(0,0) = 0.;
    M(0,1) = 0.;
    M(0,2) = -j*k_eqx*exp(-j*k_eqx*(l));
    M(0,3) = j*k_eqx*exp(j*k_eqx*(l));
// interface p
    M(1,0) = 1.;
    M(1,1) = -M11+j*k_eqx*M12/(-j*omega*rho_eq);
    M(1,2) = -M11-j*k_eqx*M12/(-j*omega*rho_eq);
// interface v
    M(2,0) = j*k_ax/(-j*omega*rho_a);
    M(2,1) = -M21+j*k_eqx*M22/(-j*omega*rho_eq);
    M(2,2) = -M21-j*k_eqx*M22/(-j*omega*rho_eq);

    cout<<"Coeff 1: "<<exp(-j*k_ax*(-l))<<endl;
    cout<<"Coeff 2: "<<-M21+j*k_eqx*M22<<endl;


    b(0) = 0.;
    b(1) = 1.;
    b(2) = j*k_ax/(-j*omega*rho_a);

    Vector3cd c = M.partialPivLu().solve(b);

    A = c(0);
    B = c(1);
    C = c(2);
    

}



xEvalPressureObliqueRigidWall::xEvalPressureObliqueRigidWall(double omega_,  double beta_, valtype M11_, valtype M12_, valtype M21_, valtype M22_) : omega(omega_), beta(beta_), M11(M11_), M12(M12_), M21(M21_), M22(M22_){
    rho_a = 1.213;

    xEvalCoeffObliqueRigidWall coeff(omega, beta, M11, M12, M21, M22);
    k_a = coeff.k_a;
    k_eq = coeff.k_eq;
    // k_eq = k_a;
    rho_eq = coeff.rho_eq;
    // rho_eq = rho_a;
    A = coeff.A;
    B = coeff.B;
    C = coeff.C;


    cout<<"Reflection A: "<<A<<endl;

}

void xEvalPressureObliqueRigidWall::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double>& res) const {

    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);
    double beta = 0.;
    valtype k_ax = k_a*cos(beta);
    valtype k_ay = k_a*sin(beta);
    valtype k_eqx = sqrt(pow(k_eq,2)-pow(k_a,2)*pow(sin(beta),2));
    valtype k_eqy = k_ay;
    // cout<<"ky: "<<k_eqy<<endl;

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    //Check if side was prescribed
    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation < 0) cdgx = -1.e-8;
    else if(orientation > 0) cdgx = 1.e-8;


    if (cdgx<0.){
        res = exp(-j*k_ax*x-j*k_ay*y)+A*exp(j*k_ax*x-j*k_ay*y);
    }
    else
    {
        res = B*exp(-j*k_eqx*x-j*k_eqy*y)+C*exp(j*k_eqx*x-j*k_eqy*y);
    }

    return;
}