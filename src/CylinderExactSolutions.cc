#include "CylinderExactSolutions.h"
#include "xMesh.h"
#include <iomanip>
#include "complex_bessel.h"
#include "AnalyticalMaterial.h"

// #define CHECK_WITH_EIGEN 1
// #ifdef CHECK_WITH_EIGEN
#include "Eigen/Dense"
// #endif

using namespace std;
using namespace xfem;
using namespace Eigen;
using namespace sp_bessel;




xEvalCylinderCoefficient::xEvalCylinderCoefficient(double omega_, double radius_, double sigmaF_, double d_) : omega(omega_), radius(radius_), sigmaF(sigmaF_), d(d_)
{
     // reference parameters Allard page 270 foam 3/4
    double phi_1 = 0.97, sigma_1 = 57e3, alpha_1 = 1.54, lambda_prime_1 = 73.8e-6, lambda_1 = 24.6e-6, rho_11 = 46, nu_1 = 0.3, E_1 = 214e3, eta_1 = 0.115;
    beta = 0.;
    Analytical_MatBiot mat(omega, beta, phi_1, sigma_1, alpha_1, lambda_prime_1, lambda_1, rho_11, nu_1, E_1, eta_1);
    k_0 = mat.k_a;
    K_eq_til = mat.K_eq_til;
    mu_1 = mat.mu_1;
    mu_2 = mat.mu_2;
    mu_3 = mat.mu_3;

    N = mat.N;
    A_hat = mat.A_hat;

    rho_a = 1.213;
    A_inc = 1.e-7;


    k_1 = mat.delta_1;
    k_2 = mat.delta_2;
    k_3 = mat.delta_3;

    Analytical_MatLimp mat_F(omega, beta, 0.04, 775e3, 1.15, 230e-6, 230e-6, 809, 0.3, 260e6, 0.5);
    valtype Z = mat_F.rho_limp_til * mat_F.c_eq_til;
    valtype k_a = mat_F.k_a;
    valtype k_ay = k_a*sin(beta);
    valtype k_eqx = sqrt(pow(mat_F.rho_limp_til,2)-pow(k_ay,2));
    
    Matrix2cd MT;
    // MT(0,0) = cos(k_eqx * d); 
    // MT(0,1) = -1. * pow(omega,2) * mat_F.rho_limp_til/k_eqx * sin(k_eqx * d); 
    // MT(1,0) = k_eqx/(pow(omega,2) * mat_F.rho_limp_til) * sin(k_eqx * d);
    // MT(1,1) = cos(k_eqx * d); 
    MT(0,0) = cos(k_eqx * d); 
    MT(0,1) = 1. * pow(omega,2) * mat_F.rho_limp_til/k_eqx * sin(k_eqx * d); 
    MT(1,0) = -k_eqx/(pow(omega,2) * mat_F.rho_limp_til) * sin(k_eqx * d);
    MT(1,1) = cos(k_eqx * d); 
    // cout<<"Coeff M_11: "<<std::setprecision(5)<<MT(0,0)<<endl;

    // cout<<"Coeff M_12: "<<std::setprecision(5)<<MT(0,1)<<endl;
    // cout<<"Coeff M_21: "<<std::setprecision(5)<<MT(1,0)<<endl;
    // cout<<"Coeff M_22: "<<std::setprecision(5)<<MT(1,1)<<endl;
    // double a = 1.;  //radius
    double a  = radius;
    inf_mode = max_mode;
    for (int mm=0; mm<inf_mode; mm++) {
        double m = 1.0*mm;
        MatrixXcd M(4,4);

        VectorXcd b(4);
        // u^t
        M(0,0) = m*hankelH2(mm,k_0*a)/a-k_0*hankelH2(mm+1,k_0*a);
        M(0,1) = (mat.mu_1*(-m*besselJ(mm, k_1*a))+k_1*mat.mu_1*a*besselJ(mm+1,k_1*a))/a;
        M(0,2) = (mat.mu_2*(-m*besselJ(mm, k_2*a))+k_2*mat.mu_2*a*besselJ(mm+1,k_2*a))/a;
        M(0,3) = -mat.mu_3*m*besselJ(mm,k_3*a)/a;
        // pressure
        M(1,0) = rho_a*pow(omega,2)*hankelH2(mm,k_0*a);
        M(1,1) = -mat.K_eq_til*mat.mu_1*pow(k_1,2)*besselJ(mm,k_1*a)-j*omega*sigmaF*d*(mat.mu_1*(-m*besselJ(mm, k_1*a))+k_1*mat.mu_1*a*besselJ(mm+1,k_1*a))/a;
        M(1,2) = -mat.K_eq_til*mat.mu_2*pow(k_2,2)*besselJ(mm,k_2*a)-j*omega*sigmaF*d*(mat.mu_2*(-m*besselJ(mm, k_2*a))+k_2*mat.mu_2*a*besselJ(mm+1,k_2*a))/a;
        M(1,3) = 0. +j*omega*sigmaF*d*mat.mu_3*m*besselJ(mm,k_3*a)/a;
        // pressure
        // M(1,0) = rho_a*pow(omega,2)*hankelH2(mm,k_0*a)+sigmaF*d*j*omega*k_0*hankelH2p(mm,k_0*a);
        // M(1,1) = -mat.K_eq_til*mat.mu_1*pow(k_1,2)*besselJ(mm,k_1*a);
        // M(1,2) = -mat.K_eq_til*mat.mu_2*pow(k_2,2)*besselJ(mm,k_2*a);
        // M(1,3) = 0. ;
        // sigma rr
        M(2,0) = 0.;
        M(2,1) = 2.*mat.N/pow(a,2)*(m*(m-1)*besselJ(mm,k_1*a)-pow(k_1,2)*pow(a,2)*besselJ(mm,k_1*a)+k_1*a*besselJ(mm+1,k_1*a))-mat.A_hat*pow(k_1,2)*besselJ(mm,k_1*a);
        M(2,2) = 2.*mat.N/pow(a,2)*(m*(m-1)*besselJ(mm,k_2*a)-pow(k_2,2)*pow(a,2)*besselJ(mm,k_2*a)+k_2*a*besselJ(mm+1,k_2*a))-mat.A_hat*pow(k_2,2)*besselJ(mm,k_2*a);
        M(2,3) = 2.*mat.N*(m*(m-1)*besselJ(mm,k_3*a)-k_3*a*m*besselJ(mm+1,k_3*a))/pow(a,2);
 
        // sigma rt
        M(3,0) = 0.;
        M(3,1) = (2.*mat.N*k_1*a*m*besselJ(mm+1,k_1*a)-mat.N*m*(2.*m-2.)*besselJ(mm,k_1*a))/pow(a,2);
        M(3,2) = (2.*mat.N*k_2*a*m*besselJ(mm+1,k_2*a)-mat.N*m*(2.*m-2.)*besselJ(mm,k_2*a))/pow(a,2);
        M(3,3) = mat.N*((-2.*pow(m,2)+2.*mm+pow(k_3,2)*pow(a,2))*besselJ(mm,k_3*a)-2.*k_3*besselJ(mm+1,k_3*a)*a)/pow(a,2);

// =============================================Generalised Limp interface====================================
        // U^t
        // M(0,0) = m*hankelH2(mm,k_0*a)/a-k_0*hankelH2(mm+1,k_0*a);
        // M(0,1) = MT(1,1)*(mat.mu_1*(-m*besselJ(mm, k_1*a))+k_1*mat.mu_1*a*besselJ(mm+1,k_1*a))/a-MT(1,0)*mat.K_eq_til*mat.mu_1*pow(k_1,2)*besselJ(mm,k_1*a);
        // M(0,2) = MT(1,1)*(mat.mu_2*(-m*besselJ(mm, k_2*a))+k_2*mat.mu_2*a*besselJ(mm+1,k_2*a))/a-MT(1,0)*mat.K_eq_til*mat.mu_2*pow(k_2,2)*besselJ(mm,k_2*a);
        // M(0,3) = -MT(1,1)*mat.mu_3*m*besselJ(mm,k_3*a)/a;
        // // pressure
        // M(1,0) = rho_a*pow(omega,2)*hankelH2(mm,k_0*a);
        // M(1,1) = -MT(0,0)*mat.K_eq_til*mat.mu_1*pow(k_1,2)*besselJ(mm,k_1*a)+MT(0,1)*(mat.mu_1*(-m*besselJ(mm, k_1*a))+k_1*mat.mu_1*a*besselJ(mm+1,k_1*a))/a;
        // M(1,2) = -MT(0,0)*mat.K_eq_til*mat.mu_2*pow(k_2,2)*besselJ(mm,k_2*a)+MT(0,1)*(mat.mu_2*(-m*besselJ(mm, k_2*a))+k_2*mat.mu_2*a*besselJ(mm+1,k_2*a))/a;
        // M(1,3) = 0. -MT(0,1)*mat.mu_3*m*besselJ(mm,k_3*a)/a;
        // // pressure
        // // M(1,0) = rho_a*pow(omega,2)*hankelH2(mm,k_0*a)+sigmaF*d*j*omega*k_0*hankelH2p(mm,k_0*a);
        // // M(1,1) = -mat.K_eq_til*mat.mu_1*pow(k_1,2)*besselJ(mm,k_1*a);
        // // M(1,2) = -mat.K_eq_til*mat.mu_2*pow(k_2,2)*besselJ(mm,k_2*a);
        // // M(1,3) = 0. ;
        // // sigma rr
        // M(2,0) = 0.;
        // M(2,1) = 2.*mat.N/pow(a,2)*(m*(m-1)*besselJ(mm,k_1*a)-pow(k_1,2)*pow(a,2)*besselJ(mm,k_1*a)+k_1*a*besselJ(mm+1,k_1*a))-mat.A_hat*pow(k_1,2)*besselJ(mm,k_1*a);
        // M(2,2) = 2.*mat.N/pow(a,2)*(m*(m-1)*besselJ(mm,k_2*a)-pow(k_2,2)*pow(a,2)*besselJ(mm,k_2*a)+k_2*a*besselJ(mm+1,k_2*a))-mat.A_hat*pow(k_2,2)*besselJ(mm,k_2*a);
        // M(2,3) = 2.*mat.N*(m*(m-1)*besselJ(mm,k_3*a)-k_3*a*m*besselJ(mm+1,k_3*a))/pow(a,2);
 
        // // sigma rt
        // M(3,0) = 0.;
        // M(3,1) = (2.*mat.N*k_1*a*m*besselJ(mm+1,k_1*a)-mat.N*m*(2.*m-2.)*besselJ(mm,k_1*a))/pow(a,2);
        // M(3,2) = (2.*mat.N*k_2*a*m*besselJ(mm+1,k_2*a)-mat.N*m*(2.*m-2.)*besselJ(mm,k_2*a))/pow(a,2);
        // M(3,3) = mat.N*((-2.*pow(m,2)+2.*mm+pow(k_3,2)*pow(a,2))*besselJ(mm,k_3*a)-2.*k_3*besselJ(mm+1,k_3*a)*a)/pow(a,2);


        MatrixXcd M0 = M;
        MatrixXcd M1 = M;
        MatrixXcd M2 = M;
        MatrixXcd M3 = M;


        b(0) = A_inc*k_0*2.*pow((-j),mm)*besselJ(mm+1,k_0*a)-A_inc*2./a*pow((-j),mm)*m*besselJ(mm,k_0*a);
        b(1) = -A_inc*2.*rho_a*pow(omega,2)*pow((-j),mm)*besselJ(mm,k_0*a);
        b(2) = 0.;
        b(3) = 0.;

        // Vector4cd c = M.fullPivHouseholderQr().solve(b);

        M0(0,0) = b(0);
        M0(1,0) = b(1);
        M0(2,0) = b(2);
        M0(3,0) = b(3);


        M1(0,1) = b(0);
        M1(1,1) = b(1);
        M1(2,1) = b(2);
        M1(3,1) = b(3);


        M2(0,2) = b(0);
        M2(1,2) = b(1);
        M2(2,2) = b(2);
        M2(3,2) = b(3);

        M3(0,3) = b(0);
        M3(1,3) = b(1);
        M3(2,3) = b(2);
        M3(3,3) = b(3);

        valtype det_M = M.determinant();
        valtype det_M0 = M0.determinant();
        valtype det_M1 = M1.determinant();
        valtype det_M2 = M2.determinant();
        valtype det_M3 = M3.determinant();

        Vector4cd c;
        c(0) = det_M0/det_M;
        c(1) = det_M1/det_M;
        c(2) = det_M2/det_M;
        c(3) = det_M3/det_M;


        if (mm==0){
            A[0][mm] = 0.5*c(0);
            A[1][mm] = 0.5*c(1);
            A[2][mm] = 0.5*c(2);
            A[3][mm] = 0.5*c(3);
        }
        else {
            A[0][mm] = c(0);
            A[1][mm] = c(1);
            A[2][mm] = c(2);
            A[3][mm] = c(3);
        }

        }
    // cout<<A[0][0]<<endl;
    // cout<<A[0][1]<<endl;
    // cout<<A[0][2]<<endl;
    // cout<<A[0][3]<<endl;
    // cout<<A[0][4]<<endl;
    // cout<<A[0][5]<<endl;
    // cout<<A[0][6]<<endl;
    // cout<<mat.mu_1*k_1*besselJp(2, k_1*a, 1)<<endl;
    // cout<<mat.mu_1*k_1*(2./(k_1*a)*besselJ(2, k_1*a)-besselJ(3, k_1*a))<<endl;
    // cout<<(mat.mu_1*(2.*besselJ(2, k_1*a))-k_1*mat.mu_1*a*besselJ(3,k_1*a))/a<<endl;
    // A[0][0] = (-4.95764487745741e-07-2.1054497856572054e-08*j); A[1][0]=(1.0054975684543387e-06-7.80622525104055e-08*j); A[]
    // A[0][1] = (-1.7945815078777731e-06+4.7461270759793677e-07*j); A[1][1] = (2.2481684716879166e-07-1.8689373540316826e-06j)
    

    }



xEvalExactCylinderUs::xEvalExactCylinderUs(double omega_, double radius_, double sigmaF_, double d_) : omega(omega_), radius(radius_), sigmaF(sigmaF_), d(d_) {
    xEvalCylinderCoefficient Coeff(omega, radius, sigmaF, d);
    inf_mode  = Coeff.inf_mode;
    for (int i=0; i<4; i++) {
        for (int j=0; j<inf_mode; j++) {
            A[i][j] = Coeff.A[i][j];
        }
    }
    // radius  = Coeff.radius;
    k_1 = Coeff.k_1;
    k_2 = Coeff.k_2;
    k_3 = Coeff.k_3;
    inf_mode = Coeff.inf_mode;


}

void xEvalExactCylinderUs::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const{

    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
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

        valtype u_r = 0.;
        valtype u_th = 0.;
        for (int m=0; m<inf_mode; m++) {
            double mm = 1.0*m;
        u_r += A[1][m]*cos(m*th)*(mm/r*besselJ(m, k_1*r)-k_1*besselJ(m+1, k_1*r))
                +A[2][m]*cos(m*th)*(mm/r*besselJ(m, k_2*r)-k_2*besselJ(m+1, k_2*r))
                +mm/r*A[3][m]*cos(m*th)*besselJ(m, k_3*r);
        u_th += -mm/r*A[1][m]*besselJ(m,k_1*r)*sin(m*th)
                -mm/r*A[2][m]*besselJ(m,k_2*r)*sin(m*th)
                -A[3][m]*sin(m*th)*(mm/r*besselJ(m,k_3*r)-k_3*besselJ(m+1, k_3*r)) ;
        }
        res(0) = u_r*cos(th)-u_th*sin(th);
        res(1) = u_r*sin(th)+u_th*cos(th);
   
    }
    else
    {
        res(0) = 0;
        res(1) = 0;


    }
    
    return;
}


xEvalExactCylinderTraction::xEvalExactCylinderTraction(double omega_, double radius_, double sigmaF_, double d_) : omega(omega_), radius(radius_), sigmaF(sigmaF_), d(d_) {
    xEvalCylinderCoefficient Coeff(omega, radius, sigmaF, d);
    inf_mode  = Coeff.inf_mode;
    for (int i=0; i<4; i++) {
        for (int j=0; j<inf_mode; j++) {
            A[i][j] = Coeff.A[i][j];
        }
    }
    // radius  = Coeff.radius;
    k_1 = Coeff.k_1;
    k_2 = Coeff.k_2;
    k_3 = Coeff.k_3;
    N = Coeff.N;
    A_hat = Coeff.A_hat;
    inf_mode = Coeff.inf_mode;
    // radius = Coeff.radius;
}

void xEvalExactCylinderTraction::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<std::complex<double> >& res) const{


    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
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
 
    if (cdgr < radius) {
        valtype e_rr = 0.;
        valtype e_rt = 0.;
        valtype e_tt = 0.;
        res(0,0) = 0. ;
        res(0,1) = 0. ;
        res(1,0) = 0. ;
        res(1,1) = 0. ;
        for (int m=0; m<inf_mode; m++) {
            double mm =1.0*m;
            e_rr += A[1][m]*cos(m*th)*(mm/r*k_1*besselJp(m,k_1*r,1)-mm/pow(r,2)*besselJ(m, k_1*r)-pow(k_1,2)*besselJp(m+1,k_1*r,1))
                    +A[2][m]*cos(m*th)*(mm/r*k_2*besselJp(m, k_2*r,1)-mm/pow(r,2)*besselJ(m, k_2*r)-pow(k_2,2)*besselJp(m+1, k_2*r,1))
                    -mm/pow(r,2)*A[3][m]*cos(m*th)*besselJ(m, k_3*r)+mm/r*A[3][m]*cos(m*th)*k_3*besselJp(m,k_3*r);
            e_tt += 1/r*(-mm*mm/r*A[1][m]*besselJ(m,k_1*r)*cos(m*th)
                    -mm*mm/r*A[2][m]*besselJ(m,k_2*r)*cos(m*th)
                    -A[3][m]*mm*cos(m*th)*(mm/r*besselJ(m,k_3*r)-k_3*besselJ(m+1, k_3*r))
                    +A[1][m]*cos(m*th)*(mm/r*besselJ(m, k_1*r)-k_1*besselJ(m+1, k_1*r))
                    +A[2][m]*cos(m*th)*(mm/r*besselJ(m, k_2*r)-k_2*besselJ(m+1, k_2*r))
                    +mm/r*A[3][m]*cos(m*th)*besselJ(m, k_3*r));
            e_rt += 0.5*(1/r*(-mm*A[1][m]*sin(m*th)*(mm/r*besselJ(m, k_1*r)-k_1*besselJ(m+1, k_1*r))
                    -mm*A[2][m]*sin(m*th)*(mm/r*besselJ(m, k_2*r)-k_2*besselJ(m+1, k_2*r))+
                    -mm*mm/r*A[3][m]*sin(m*th)*besselJ(m, k_3*r))
                    -mm/r*A[1][m]*k_1*besselJp(m,k_1*r,1)*sin(m*th)+mm/pow(r,2)*A[1][m]*besselJ(m,k_1*r)*sin(m*th)
                    -mm/r*A[2][m]*k_2*besselJp(m,k_2*r,1)*sin(m*th)+mm/pow(r,2)*A[2][m]*besselJ(m,k_2*r)*sin(m*th)
                    -A[3][m]*sin(m*th)*(mm*k_3/r*besselJp(m,k_3*r,1)-mm/pow(r,2)*besselJ(m,k_3*r)-pow(k_3,2)*besselJp(m+1, k_3*r,1))
                    -(-mm/pow(r,2)*A[1][m]*besselJ(m,k_1*r)*sin(m*th)
                    -mm/pow(r,2)*A[2][m]*besselJ(m,k_2*r)*sin(m*th)
                    -A[3][m]/r*sin(m*th)*(mm/r*besselJ(m,k_3*r)-k_3*besselJ(m+1, k_3*r))));


            
        }
            res(0,0) += 2.0*N*e_rr+A_hat*(e_rr+e_tt);
            res(1,1) += 2.0*N*e_tt+A_hat*(e_rr+e_tt);
            res(0,1) += 2.0*N*e_rt;
            res(1,0) += res(0,1);
    }
    else {
            res(0,0) = 0. ;
            res(0,1) = 0. ;
            res(1,0) = 0. ;
            res(1,1) = 0. ; 
        }
    

    
    return;
}


xEvalExactCylinderUtotal::xEvalExactCylinderUtotal(double omega_, double radius_, double sigmaF_, double d_) : omega(omega_), radius(radius_), sigmaF(sigmaF_), d(d_) {
    xEvalCylinderCoefficient Coeff(omega, radius, sigmaF, d);
    inf_mode  = Coeff.inf_mode;
    for (int i=0; i<4; i++) {
        for (int j=0; j<inf_mode; j++) {
            A[i][j] = Coeff.A[i][j];
        }
    }

    k_1 = Coeff.k_1;
    k_2 = Coeff.k_2;
    k_3 = Coeff.k_3;
    mu_1 = Coeff.mu_1;
    mu_2 = Coeff.mu_2;
    mu_3 = Coeff.mu_3;

    A_inc = Coeff.A_inc;
    rho_a  = Coeff.rho_a;
    k_0 = Coeff.k_0;
    // radius = Coeff.radius;
}

void xEvalExactCylinderUtotal::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const{

    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
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

    if (cdgr<radius)
    {
        res(0) = 0.;
        res(1) = 0.;
        for (int m=0; m<inf_mode; m++) {
            double mm = 1.0*m;

        res(0) += mu_1*A[1][m]*cos(m*th)*(mm/r*besselJ(m, k_1*r)-k_1*besselJ(m+1, k_1*r))
                +mu_2*A[2][m]*cos(m*th)*(mm/r*besselJ(m, k_2*r)-k_2*besselJ(m+1, k_2*r))
                +mu_3*mm/r*A[3][m]*cos(m*th)*besselJ(m, k_3*r);

        res(1) += -mu_1*mm/r*A[1][m]*besselJ(m,k_1*r)*sin(m*th)
                -mu_2*mm/r*A[2][m]*besselJ(m,k_2*r)*sin(m*th)
                -mu_3*A[3][m]*sin(m*th)*(mm/r*besselJ(m,k_3*r)-k_3*besselJ(m+1, k_3*r)) ;
        }
    }
    else
    {
        res(0) = 0;
        res(1) = 0;
        valtype res11 = 0.;
        valtype res12 = 0.;
        valtype res21 = 0.;
        valtype res22 = 0.;

        for (int mm=0; mm<inf_mode; mm++)
        {
        double m = 1.0*mm;
 
        if (mm==0) {
            res11 += 0;
            res21 += 0;
        }
        else {
            res11 += pow((-j),mm)*(m/r*besselJ(mm,k_0*r)-k_0*besselJ(mm+1,k_0*r))*cos(mm*th);
            res21 += -m/r*pow((-j),mm)*besselJ(mm,k_0*r)*sin(mm*th);
        }

        res12 += A[0][mm]*(m/r*hankelH2(mm, k_0*r)-k_0*hankelH2(mm+1,k_0*r))*cos(mm*th);
        res22 += -m/r*A[0][mm]*hankelH2(mm, k_0*r)*sin(mm*th);
        // cout<<res<<endl;
        }
        res(0) +=  -A_inc*k_0*besselJ(1,k_0*r)+A_inc*2.*res11+res12;

        res(1) += A_inc*2.*res21+res22;
        // res(0) = 0;
        // res(1) = 0;
        
    }
    
    return;
}

xEvalExactCylinderPressure::xEvalExactCylinderPressure(double omega_, double radius_, double sigmaF_, double d_) : omega(omega_), radius(radius_), sigmaF(sigmaF_), d(d_) {
    xEvalCylinderCoefficient Coeff(omega, radius, sigmaF, d);
    
    inf_mode  = Coeff.inf_mode;
    for (int i=0; i<4; i++) {
        for (int j=0; j<inf_mode; j++) {
            A[i][j] = Coeff.A[i][j];
        }
    }

    k_1 = Coeff.k_1;
    k_2 = Coeff.k_2;
    k_3 = Coeff.k_3;
    K_eq_til = Coeff.K_eq_til;
    mu_1 = Coeff.mu_1;
    mu_2 = Coeff.mu_2;
    mu_3 = Coeff.mu_3;

    A_inc = Coeff.A_inc;
    rho_a  = Coeff.rho_a;
    k_0 = Coeff.k_0;

    // radius = Coeff.radius;


    // cout<<"derivative of bessel"<<besselJ(0,j,1)<<endl;
  

}

void xEvalExactCylinderPressure::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double>& res) const{
    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    double r = sqrt(x * x + y * y);
    double th = atan2(y, x);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    double cdgy = cdgxyz(1);

    double cdgr = sqrt(cdgx * cdgx + cdgy * cdgy);

    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation < 0) cdgr = radius-1.e-8;
    else if(orientation > 0) cdgr = radius+1.e-8;


    if (cdgr<radius) {
        res = 0.;
        for (int mm=0; mm<inf_mode; mm++)
        {
        double m = 1.0*mm;
        res += K_eq_til*mu_1*pow(k_1,2)*besselJ(mm,k_1*r)*cos(m*th)*A[1][mm]+K_eq_til*mu_2*pow(k_2,2)*besselJ(mm,k_2*r)*cos(m*th)*A[2][mm];
        // res += K_eq_til*mu_1*pow(k_1,2)*besselJ(mm,k_1*r)*cos(m*th)+K_eq_til*mu_2*pow(k_2,2)*besselJ(mm,k_2*r)*cos(m*th);
        // cout<<Coeff.A[1][mm]<<endl;
        // cout<<res<<endl;
        }    
        
    }

    else
    {
        res = 0.;
        valtype res1 = 0.;
        valtype res2 = 0.;

        for (int mm=0; mm<inf_mode; mm++)
        {
        if (mm==0) {
            res1 += 0;
        }
        else {
            res1 += pow((-j),mm)*besselJ(mm,k_0*r)*cos(mm*th);
        }
        
        res2 += A[0][mm]*hankelH2(mm, k_0*r)*cos(mm*th);
        // cout<<res<<endl;
        }
        res += A_inc*rho_a*pow(omega,2)*besselJ(0,k_0*r)+A_inc*rho_a*pow(omega,2)*2.*res1+rho_a*pow(omega,2)*res2;
    }
    
    return;
}


xEvalExactCylinderVelocity::xEvalExactCylinderVelocity(double omega_, double radius_, double sigmaF_, double d_) : omega(omega_), radius(radius_), sigmaF(sigmaF_), d(d_) {
    xEvalCylinderCoefficient Coeff(omega, radius, sigmaF, d);
    
    inf_mode  = Coeff.inf_mode;
    for (int i=0; i<4; i++) {
        for (int j=0; j<inf_mode; j++) {
            A[i][j] = Coeff.A[i][j];
        }
    }
    k_1 = Coeff.k_1;
    k_2 = Coeff.k_2;
    k_3 = Coeff.k_3;
    K_eq_til = Coeff.K_eq_til;
    mu_1 = Coeff.mu_1;
    mu_2 = Coeff.mu_2;
    mu_3 = Coeff.mu_3;

    A_inc = Coeff.A_inc;
    rho_a  = Coeff.rho_a;
    k_0 = Coeff.k_0;
    cout<<"wave number 1"<<k_1<<endl;
    cout<<"wave number 2"<<k_2<<endl;
    cout<<"wave number 3"<<k_3<<endl;
    // radius = Coeff.radius;


}

void xEvalExactCylinderVelocity::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const{
    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    double r = sqrt(x * x + y * y);
    double th = atan2(y, x);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    double cdgy = cdgxyz(1);

    double cdgr = sqrt(cdgx * cdgx + cdgy * cdgy);

    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation < 0) cdgr = radius-1.e-8;
    else if(orientation > 0) cdgr = radius+1.e-8;

    

    if (cdgr<radius) {
        res(0) = 0.;
        res(1) = 0.;
        // for (int mm=0; mm<inf_mode; mm++)
        // {
        // double m = 1.0*mm;
        // res(0) += K_eq_til*mu_1*pow(k_1,2)*besselJp(mm,k_1*r,1)*cos(m*th)*A[1][mm]+K_eq_til*mu_2*pow(k_2,2)*besselJp(mm,k_2*r,1)*cos(m*th)*A[2][mm];
        // res(1) += -m*K_eq_til*mu_1*pow(k_1,2)*besselJ(mm,k_1*r)*sin(m*th)*A[1][mm]-m*K_eq_til*mu_2*pow(k_2,2)*besselJ(mm,k_2*r)*sin(m*th)*A[2][mm];
        // cout<<Coeff.A[1][mm]<<endl;
        // cout<<res<<endl;
        // }    
        
    }

    else
    {
   
        valtype res11 = 0.;
        valtype res12 = 0.;
        valtype res21 = 0.;
        valtype res22 = 0.;
        valtype vrr = 0.;
        valtype vrt = 0.;

        for (int mm=0; mm<inf_mode; mm++)
        {
        double m = 1.0*mm;
 
        if (mm==0) {
            res11 += 0;
            res21 += 0;
        }
        else {
            res11 += pow((-j),mm)*k_0*besselJp(m,k_0*r,1)*cos(mm*th);
            
            res21 += -m/r*pow((-j),mm)*besselJ(mm,k_0*r)*sin(mm*th);
        }

        res12 += k_0*A[0][mm]*hankelH2p(m, k_0*r,1)*cos(mm*th);
        res22 += -m/r*A[0][mm]*hankelH2(mm, k_0*r)*sin(mm*th);
        // cout<<res<<endl;
        }
        vrr += k_0*rho_a*pow(omega,2)*A_inc*besselJp(0,k_0*r,1)+A_inc*rho_a*pow(omega,2)*2.*res11+rho_a*pow(omega,2)*res12;

        vrt += A_inc*rho_a*pow(omega,2)*2.*res21+rho_a*pow(omega,2)*res22;
        res(0) = vrr*cos(th)-vrt*sin(th);
        res(1) = vrr*sin(th)+vrt*cos(th);
    }
    
    return;
}


xEvalExactCylinderVelocityS::xEvalExactCylinderVelocityS(double omega_, double radius_, double sigmaF_, double d_) : omega(omega_), radius(radius_), sigmaF(sigmaF_), d(d_) {
    xEvalCylinderCoefficient Coeff(omega, radius, sigmaF, d);
    
    inf_mode  = Coeff.inf_mode;
    for (int i=0; i<4; i++) {
        for (int j=0; j<inf_mode; j++) {
            A[i][j] = Coeff.A[i][j];
        }
    }
    k_1 = Coeff.k_1;
    k_2 = Coeff.k_2;
    k_3 = Coeff.k_3;
    K_eq_til = Coeff.K_eq_til;
    mu_1 = Coeff.mu_1;
    mu_2 = Coeff.mu_2;
    mu_3 = Coeff.mu_3;

    A_inc = Coeff.A_inc;
    rho_a  = Coeff.rho_a;
    k_0 = Coeff.k_0;
    // radius = Coeff.radius;
    cout<<"wave number 1"<<k_1<<endl;
    cout<<"wave number 2"<<k_2<<endl;
    cout<<"wave number 3"<<k_3<<endl;
    // cout<<"Coefficient are:"<<A[0][59]<<endl;
    // cout<<"Coefficient are:"<<A[1][59]<<endl;
    // cout<<"Coefficient are:"<<A[2][59]<<endl;
    // cout<<"Coefficient are:"<<A[3][59]<<endl;


}

void xEvalExactCylinderVelocityS::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double> & res) const{
    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    double r = sqrt(x * x + y * y);
    double th = atan2(y, x);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    double cdgy = cdgxyz(1);

    double cdgr = sqrt(cdgx * cdgx + cdgy * cdgy);

    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation < 0) cdgr = radius-1.e-8;
    else if(orientation > 0) cdgr = radius+1.e-8;

    

    if (cdgr<radius) {
        res = 0.;
        // for (int mm=0; mm<inf_mode; mm++)
        // {
        // double m = 1.0*mm;
        // res(0) += K_eq_til*mu_1*pow(k_1,2)*besselJp(mm,k_1*r,1)*cos(m*th)*A[1][mm]+K_eq_til*mu_2*pow(k_2,2)*besselJp(mm,k_2*r,1)*cos(m*th)*A[2][mm];
        // res(1) += -m*K_eq_til*mu_1*pow(k_1,2)*besselJ(mm,k_1*r)*sin(m*th)*A[1][mm]-m*K_eq_til*mu_2*pow(k_2,2)*besselJ(mm,k_2*r)*sin(m*th)*A[2][mm];
        // cout<<Coeff.A[1][mm]<<endl;
        // cout<<res<<endl;
        // }    
        
    }

    else
    {
        // res = 0.;
        valtype res11 = 0.;
        valtype res12 = 0.;

        for (int mm=0; mm<inf_mode; mm++)
        {
        double m = 1.0*mm;
 
        if (mm==0) {
            res11 += 0;

        }
        else {
            res11 += pow((-j),mm)*k_0*besselJp(m,k_0*r,1)*cos(mm*th);
            // res21 += -m/r*pow((-j),mm)*besselJ(mm,k_0*r)*sin(mm*th);
        }

        res12 += k_0*A[0][mm]*hankelH2p(m, k_0*r,1)*cos(mm*th);
        // res22 += -m/r*A[0][mm]*hankelH2(mm, k_0*r)*sin(mm*th);
        // cout<<res<<endl;
        }
        res = k_0*rho_a*pow(omega,2)*A_inc*besselJp(0,k_0*r,1)+A_inc*rho_a*pow(omega,2)*2.*res11+rho_a*pow(omega,2)*res12;

        // res(1) += A_inc*rho_a*pow(omega,2)*2.*res21+rho_a*pow(omega,2)*res22;
        
    }
    
    return;
}


