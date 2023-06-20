#include "BiotexactSolutions.h"
#include "AnalyticalMaterial.h"
#include "xMesh.h"
#include <iomanip>
#include <cmath>
#include "complex_bessel.h"

#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"


using namespace std;
using namespace xfem;
using namespace Eigen;



// double hankel(double v, double x)
// {
//     double k;
//     const std::complex<double> j{0.,1.};
//     k = cyl_bessel_j(v,x)-j*cyl_neumann(v,x);
//     return k;
// }

xEvalCoefficient::xEvalCoefficient(double omega_, double beta_, double sigmaF_, double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_)
{
    const double P = 101325., gamma = 1.4, rho_a = 1.213, mu = 0.1839e-4, Pr = 0.710;
    const double phi = 0.97, sigma = 57e3, alpha = 1.54, lambda_prime = 73.8e-6, lambda = 24.6e-6;
    
    double K_a = gamma * P;
    double c_a = sqrt(K_a/rho_a);
    k_a = omega/c_a;
    double Z_a = rho_a * c_a;
    double nu = mu / rho_a;
    double nu_prime = nu / Pr;

    //// fluid equivalent part######################
    valtype omega_0 = sigma * phi / (rho_a * alpha);
    valtype omega_inf = pow(sigma, 2)*pow(phi,2)*pow(lambda,2)/(4.*mu*rho_a*pow(alpha,2));
    valtype rho_eq_til  = (rho_a*alpha/phi)*(1.+((omega_0)/(j*omega))*sqrt(1.+j*omega/omega_inf));

    valtype omega_prime_infty = (16.*nu_prime)/(pow(lambda_prime, 2));
    valtype F_prime_CA = sqrt(1.+j*omega/omega_prime_infty);
    valtype alpha_prime_til = 1.+omega_prime_infty*F_prime_CA/(2.*j*omega);
    valtype K_eq_til = (gamma*P/phi)/(gamma-(gamma-1.)/alpha_prime_til);
    valtype c_eq = sqrt(K_eq_til/rho_eq_til);
    valtype k_eq_til = omega / c_eq;
    // cout<<"wave number is: "<<k_eq_til<<endl;

    ////// solid part////////
    const double rho_1 = 46., nu_p = 0.3, E = 214e3, eta = 0.115;
    valtype loss = 1. + j*eta;
    valtype rho_12 = -phi*rho_a*(alpha-1.);
    valtype rho_11 = rho_1-rho_12;
    valtype rho_2 = phi*rho_a;
    valtype rho_22 = rho_2-rho_12;

    valtype rho_22_til = pow(phi,2)*rho_eq_til;
    valtype rho_12_til = rho_2-rho_22_til;
    valtype rho_11_til = rho_1-rho_12_til;
    valtype rho_til = rho_11_til-(pow(rho_12_til, 2)/rho_22_til);
    valtype gamma_til = phi*(rho_12_til/rho_22_til-(1.-phi)/phi);
    valtype rho_s_til =  rho_til+pow(gamma_til,2)*rho_eq_til;

    N = E*loss/(2.*(1.+nu_p));
    A_hat = nu_p*E*loss/((1.+nu_p)*(1.-2.*nu_p));
    P_hat = A_hat + 2.*N;

    // biot 1956 elastic coefficients
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

    mu_1 = gamma_til*pow(delta_eq,2)/(pow(delta_1,2)-pow(delta_eq,2));
    mu_2 = gamma_til*pow(delta_eq,2)/(pow(delta_2,2)-pow(delta_eq,2));
    mu_3 = -gamma_til;

    double ky = k_a*sin(beta);
    double k_ax = k_a*cos(beta);
    valtype k_1x = sqrt(pow(delta_1,2)-pow(ky,2));
    valtype k_2x = sqrt(pow(delta_2,2)-pow(ky,2));
    valtype k_3x = sqrt(pow(delta_3,2)-pow(ky,2));

    Matrix4cd M;
    Matrix4cd M1;
    Matrix4cd M2;
    Matrix4cd M3;
    Matrix4cd M4;
    Vector4cd b;
// =====================================perfect interface air-Biot coupling=================================
    // M << k_ax/(rho_a*pow(omega,2)), k_1x/(K_eq_til*pow(delta_1,2)), k_2x/(K_eq_til*pow(delta_2,2)), ky*mu_3,
    //     1., -1. ,-1., 0.,
    //     0., 2.*N*pow(k_1x,2)/(mu_1*K_eq_til*pow(delta_1,2))+A_hat/(mu_1*K_eq_til), 2.*N*pow(k_2x,2)/(mu_2*K_eq_til*pow(delta_2,2))+A_hat/(mu_2*K_eq_til), 2.*N*ky*k_3x,
    //     0., -2.*N*ky*k_1x/(mu_1*K_eq_til*pow(delta_1,2)), -2.*N*ky*k_2x/(mu_2*K_eq_til*pow(delta_2,2)), N*pow(k_3x,2)-N*pow(ky,2);

// ================================Generalized interface coefficients==================================
// film is modelled by fluid with limp propriété

    // Analytical_MatLimp mat_F(omega, beta, 0.04, 775e3, 1.15, 230e-6, 230e-6, 809, 0.3, 260e6, 0.5);
    // valtype Z = mat_F.rho_limp_til * mat_F.c_eq_til;
    // valtype k_a = mat_F.k_a;
    // valtype k_ay = k_a*sin(beta);
    // valtype k_eqx = sqrt(pow(mat_F.rho_limp_til,2)-pow(k_ay,2));
    
    // Matrix2cd MT;
    // MT(0,0) = cos(k_eqx * d); 
    // MT(0,1) = -1. * pow(omega,2) * mat_F.rho_limp_til/k_eqx * sin(k_eqx * d); 
    // MT(1,0) = k_eqx/(pow(omega,2) * mat_F.rho_limp_til) * sin(k_eqx * d);
    // MT(1,1) = cos(k_eqx * d); 

    // // MT(0,0) = 1.; 
    // // MT(0,1) = j*omega*sigmaF*d; 
    // // MT(1,0) = 0.;
    // // MT(1,1) = 1.; 


    // M << j*k_ax/(rho_a*pow(omega,2)), -MT(1,0)+MT(1,1)*j*k_1x/(K_eq_til*pow(delta_1,2)), -MT(1,0)+MT(1,1)*j*k_2x/(K_eq_til*pow(delta_2,2)), MT(1,1)*j*ky*mu_3,
    //     1., -MT(0,0)+MT(0,1)*j*k_1x*mu_1/(mu_1*K_eq_til*pow(delta_1,2)),-MT(0,0)+MT(0,1)*j*k_2x*mu_2/(mu_2*K_eq_til*pow(delta_2,2)), 0.+MT(0,1)*j*ky*mu_3,
    //     0., 2.*N*pow(k_1x,2)/(mu_1*K_eq_til*pow(delta_1,2))+A_hat/(mu_1*K_eq_til), 2.*N*pow(k_2x,2)/(mu_2*K_eq_til*pow(delta_2,2))+A_hat/(mu_2*K_eq_til), 2.*N*ky*k_3x,
    //     0., -2.*N*ky*k_1x/(mu_1*K_eq_til*pow(delta_1,2)), -2.*N*ky*k_2x/(mu_2*K_eq_til*pow(delta_2,2)), N*pow(k_3x,2)-N*pow(ky,2);

    // b << j*k_ax/(rho_a*pow(omega,2)), -1, 0., 0.;
    // Vector4cd c = M.partialPivLu().solve(b);

// ========================================================Biot modelled film=========================================================
#if 1
    Analytical_MatBiot mat_F(omega, beta, 0.04, 775e3, 1.15, 230e-6, 230e-6, 809, 0.3, 260e6, 0.5);
    valtype k_a = mat_F.k_a;
    valtype k_ay = k_a*sin(beta);
    valtype k_eqx = sqrt(pow(mat_F.k_eq_til,2)-pow(k_ay,2));

    MatrixXcd TM_(6,6);
    TM_ <<  0., 0., 0., j*ky*mat_F.A_hat/mat_F.P_hat, j*ky*mat_F.gamma_til, -(pow(mat_F.A_hat,2)-pow(mat_F.P_hat,2))/mat_F.P_hat*pow(ky,2)-mat_F.rho_til*pow(omega,2),
        0., 0., 0., 1./mat_F.P_hat, 0., j*ky*mat_F.A_hat/mat_F.P_hat,
        0., 0., 0., 0., -1./mat_F.K_eq_til+pow(ky,2)/(mat_F.rho_eq_til*pow(omega,2)), -j*ky*mat_F.gamma_til,
        j*ky, -mat_F.rho_s_til*pow(omega,2), -mat_F.rho_eq_til*mat_F.gamma_til*pow(omega,2), 0., 0., 0.,
        0., mat_F.rho_eq_til*mat_F.gamma_til*pow(omega,2), mat_F.rho_eq_til*pow(omega,2), 0., 0., 0.,
        1./mat_F.N, j*ky, 0., 0., 0., 0.;

    MatrixXcd I(6,6);
    I = MatrixXcd::Identity(6,6);
    MatrixXcd TM(6,6);
    // TM = (-d*TM_).exp();
    TM = I-d*TM_;

    // cout<<"TM(3,1)"<<TM(3,1)<<endl;

    Matrix4cd SV;
    SV(0,0) = 1.;
    SV(0,1) = -TM(4,1)*(-j*k_1x)/(K_eq_til*mu_1*pow(delta_1,2))-TM(4,2)*(-j*k_1x/(K_eq_til*pow(delta_1,2)))-TM(4,4);
    SV(0,2) = -TM(4,1)*(-j*k_2x/(K_eq_til*mu_2*pow(delta_2,2)))-TM(4,2)*(-j*k_2x/(K_eq_til*pow(delta_2,2)))-TM(4,4);
    SV(0,3) = -TM(4,1)*(-j*ky)-TM(4,2)*(-j*mu_3*ky);
    // # U^t_x
    SV(1,0) = j*k_ax/(rho_a*pow(omega,2));
    SV(1,1) = -TM(2,2)*(-j*k_1x/(K_eq_til*pow(delta_1,2)))-TM(2,4)-TM(2,5)*(-j*ky)/(K_eq_til*mu_1*pow(delta_1,2));
    SV(1,2) = -TM(2,2)*(-j*k_2x/(K_eq_til*pow(delta_2,2)))-TM(2,4)-TM(2,5)*(-j*ky)/(K_eq_til*mu_2*pow(delta_2,2));
    SV(1,3) = -TM(2,2)*(-j*mu_3*ky)-TM(2,5)*(j*k_3x);
    // # sigma_xx = 0
    SV(2,0) = 0.;
    SV(2,1) = -TM(3,0)*(-2.*N*k_1x*ky/(K_eq_til*mu_1*pow(delta_1,2)))-TM(3,1)*(-j*k_1x)/(K_eq_til*mu_1*pow(delta_1,2))-TM(3,2)*(-j*k_1x/(K_eq_til*pow(delta_1,2)))\
                -TM(3,3)*(((2.*N+A_hat)*-pow(k_1x,2)/(K_eq_til*mu_1*pow(delta_1,2)))+A_hat*-pow(ky,2)/(K_eq_til*mu_1*pow(delta_1,2)));
    SV(2,2) = -TM(3,0)*(-2.*N*k_2x*ky/(K_eq_til*mu_2*pow(delta_2,2)))-TM(3,1)*(-j*k_2x/(K_eq_til*mu_2*pow(delta_2,2)))-TM(3,2)*(-j*k_2x/(K_eq_til*pow(delta_2,2)))\
        -TM(3,3)*((2.*N+A_hat)*-pow(k_2x,2)/(K_eq_til*mu_2*pow(delta_2,2))+A_hat*-pow(ky,2)/(K_eq_til*mu_2*pow(delta_2,2)));
    SV(2,3) = -TM(3,0)*(N*(pow(k_3x,2)-pow(ky,2)))-TM(3,1)*(-j*ky)-TM(3,2)*(-j*mu_3*ky)\
        -TM(3,3)*((2.*N+A_hat)*-k_3x*ky+A_hat*k_3x*ky);
    // # sigma_yx = 0
    SV(3,0) = 0.;
    SV(3,1) = -TM(0,0)*(N*-2.*k_1x*ky/(K_eq_til*mu_1*pow(delta_1,2)))-TM(0,3)*(((2.*N+A_hat)*-pow(k_1x,2)/(K_eq_til*mu_1*pow(delta_1,2)))+A_hat*-pow(ky,2)/(K_eq_til*mu_1*pow(delta_1,2)))\
        -TM(0,4)-TM(0,5)*(-j*ky)/(K_eq_til*mu_1*pow(delta_1,2));
    SV(3,2) = -TM(0,0)*(N*-2.*k_2x*ky/(K_eq_til*mu_2*pow(delta_2,2)))-TM(0,3)*((2.*N+A_hat)*-pow(k_2x,2)/(K_eq_til*mu_2*pow(delta_2,2))+A_hat*-pow(ky,2)/(K_eq_til*mu_2*pow(delta_2,2)))\
        -TM(0,4)-TM(0,5)*(-j*ky)/(K_eq_til*mu_2*pow(delta_2,2));
    SV(3,3) = -TM(0,0)*(N*(pow(k_3x,2)-pow(ky,2)))-TM(0,3)*((2.*N+A_hat)*-k_3x*ky+A_hat*k_3x*ky)-TM(0,5)*(j*k_3x);

    b << -1., j*k_ax/(rho_a*pow(omega,2)), 0., 0.;
    Vector4cd c = SV.partialPivLu().solve(b);
#endif
    R = c(0);
    T1 = c(1);
    T2 = c(2);
    T3 = c(3);

//  ========================================= PRESSURE JUMP COEFFICIENT======================
// The linear solver of eigen is not accurate to calculate the case air-Biot coupling with pressure jump cofficient
    // beta = 45, sigmaF = 775e3, d=0.45e-3  f=20000
    // R = 0.4736916892837464-0.06303193325800326*j;
    // T1 = -0.029672872289803188-0.08721120331262953*j;
    // T2 = 0.8080613937807499-0.059091888644636215*j;
    // T3 = 3.3741011139829993e-12-1.0840356997717598e-12*j;

        // beta = 45, sigmaF = 775e3, d=0.q5e-3 f=1256
    // R = (0.3122292673094504-0.08930058599790822*j);
    // T1 = (-0.017595589988612565-0.02773306107603979*j);
    // T2 = (0.9209502077093961-0.11465606747313398*j);
    // T3 = (3.1398268579164395e-10-1.8836990688845787e-10*j);

        // beta = 30, sigmaF = 775e3, d=0.45e-3  f=12560
    // R = (0.3235007219542232-j*0.031218635645836317);
    // T1 = (0.0009396093076229649-0.0020736115846061967*j);
    // T2 = (0.830000749291403-0.05187536985197039*j);
    // T3 = (1.02413173010619e-12-3.528170974331122e-13*j);

        // beta = 30, sigmaF = 775e3, d=0.45e-3  f=31400
    // R = (0.1182301660077634+0.0018358028827327863*j);
    // T1 = (0.0005380072840555999-0.000694145874995288*j);
    // T2 = (0.7470227601773236+0.0033016648271148863*j);
    // T3 = (1.4768396154891418e-13-5.0715309987938206e-14*j);

    // Generalized interface beta = 30, f=31400, d=0.45e-3 non-woven film
    // R = (0.37909116328790593+0.2060150614727343*j);
    // T1 = (0.00028269073177113587-0.0009589912162235677*j);
    // T2 = (0.7080458757137789-0.2972827189745256*j);
    // T3 = (9.744284375508612e-14-7.079468579311765e-14*j);

    // Generalized interface where the film is modelled by Biot waved
    // R = (0.3656407593630336+0.20951681827406377*j);
    // T1 = (-0.0004311789691831611+0.001055413510124234*j);
    // T2 = (0.7232787268627336-0.3022186245170058*j);
    // T3 = (6.613817867080286e-13-4.2560709499794784e-13*j);

    // Generalized interface where the film is modelled by Biot waved
    // f = 50 Hz
    // R = (0.4661314630961826-0.10436733379688035*j);
    // T1 = (0.23289537285674863-0.6921219634151338*j);
    // T2 = (0.968983207343647+0.4566822066993247*j);
    // T3 = (1.3139084312924968e-08+7.029733103639829e-10*j);

    //*******************************************//
    A_1 = T1/(K_eq_til*mu_1*pow(delta_1, 2));
    A_2 = T2/(K_eq_til*mu_2*pow(delta_2, 2));
    // cout<<"Reflection : "<<std::setprecision(20)<<R<<endl;

    // cout<<"Transmission 1: "<<std::setprecision(20)<<T1<<endl;
    // cout<<"Transmission 2: "<<std::setprecision(20)<<T2<<endl;
    // cout<<"Transmission 3: "<<std::setprecision(20)<<T3<<endl;

   
}


xEvalExactBiotUs::xEvalExactBiotUs(double omega_, double beta_, double sigmaF_, double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_) {
    xEvalCoefficient Coeff(omega, beta, sigmaF, d);
    A_1 = Coeff.A_1;
    A_2 = Coeff.A_2;
    T3 = Coeff.T3;
    k_a = Coeff.k_a;
    delta_1 = Coeff.delta_1;
    delta_2 = Coeff.delta_2;
    delta_3 = Coeff.delta_3;
}

void xEvalExactBiotUs::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const{

    const double rho_a = 1.213;

    double ky = k_a*sin(beta);
    double k_ax = k_a*cos(beta);
    valtype k_1x = sqrt(pow(delta_1,2)-pow(ky,2));
    valtype k_2x = sqrt(pow(delta_2,2)-pow(ky,2));
    valtype k_3x = sqrt(pow(delta_3,2)-pow(ky,2));

    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    // u^s  to use #ifdef endif
    // res(0) = -j*k_1x*exp(-j*k_1x*x-j*ky*y) -j*k_2x*exp(-j*k_2x*x-j*ky*y) -j*ky*exp(-j*k_3x*x-j*ky*y);
    // res(1) = -j*ky*exp(-j*k_1x*x-j*ky*y) -j*ky*exp(-j*k_2x*x-j*ky*y) +j*k_3x*exp(-j*k_3x*x-j*ky*y);
    
    // res(0) = -j*k_1x*exp(-j*k_1x*x-j*ky*y) -j*k_2x*exp(-j*k_2x*x-j*ky*y);
    // res(1) = -j*ky*exp(-j*k_1x*x-j*ky*y) -j*ky*exp(-j*k_2x*x-j*ky*y);
    
    // res(0) = -j*k_1x*exp(-j*k_1x*x-j*ky*y);
    // res(1) = -j*ky*exp(-j*k_1x*x-j*ky*y);
    
    // res(0) = -j*k_2x*exp(-j*k_2x*x-j*ky*y);
    // res(1) = -j*ky*exp(-j*k_2x*x-j*ky*y);
    
    // res(0) = -1e-5*j*ky*exp(-j*k_3x*x-j*ky*y);
    // res(1) = 1e-5*j*k_3x*exp(-j*k_3x*x-j*ky*y);

    // 2D computation numerically
    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);
    if(orientation < 0) cdgx = -1.e-8;
    else if(orientation > 0) cdgx = 1.e-8;
    if (cdgx<0.){
        res(0) = 0.;
        res(1) = 0.;
    }
    else
    {
        res(0) = (-j*k_1x*A_1*exp(-j*k_1x*x-j*ky*y)-j*k_2x*A_2*exp(-j*k_2x*x-j*ky*y)-j*ky*T3*exp(-j*k_3x*x-j*ky*y));
        res(1) = (-j*ky*A_1*exp(-j*k_1x*x-j*ky*y)-j*ky*A_2*exp(-j*k_2x*x-j*ky*y) +j*k_3x*T3*exp(-j*k_3x*x-j*ky*y));
    }
    
    return;
}


xEvalExactBiotTraction::xEvalExactBiotTraction(double omega_, double beta_, double sigmaF_, double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_) {
    xEvalCoefficient Coeff(omega, beta, sigmaF, d);
    A_1 = Coeff.A_1;
    A_2 = Coeff.A_2;
    T3 = Coeff.T3;
    k_a = Coeff.k_a;
    delta_1 = Coeff.delta_1;
    delta_2 = Coeff.delta_2, 
    delta_3 = Coeff.delta_3;
    N = Coeff.N;
    A_hat = Coeff.A_hat;
    P_hat = Coeff.P_hat;
}

void xEvalExactBiotTraction::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<std::complex<double> >& res) const{

    const double rho_a = 1.213;
    // valtype e_xx = -pow(k_1x,2)*exp(-j*k_1x*x-j*ky*y)-pow(k_2x,2)*exp(-j*k_2x*x-j*ky*y);
    // valtype e_yy = -pow(ky,2)*exp(-j*k_1x*x-j*ky*y)-pow(ky,2)*exp(-j*k_2x*x-j*ky*y);
    // valtype e_xy = 0.5*(-k_1x*ky*exp(-j*k_1x*x-j*ky*y)-k_2x*ky*exp(-j*k_2x*x-j*ky*y)
    //                 -k_1x*ky*exp(-j*k_1x*x-j*ky*y)-k_2x*ky*exp(-j*k_2x*x-j*ky*y));
    double ky = k_a*sin(beta);
    double k_ax = k_a*cos(beta);
    valtype k_1x = sqrt(pow(delta_1,2)-pow(ky,2));
    valtype k_2x = sqrt(pow(delta_2,2)-pow(ky,2));
    valtype k_3x = sqrt(pow(delta_3,2)-pow(ky,2));

    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    // 2D computation numerically
    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation > 0) cdgx = -1.e-8;
    else if(orientation < 0) cdgx = 1.e-8;
 
    valtype e_xx = (-A_1*pow(k_1x,2)*exp(-j*k_1x*x-j*ky*y)-A_2*pow(k_2x,2)*exp(-j*k_2x*x-j*ky*y)-T3*k_3x*ky*exp(-j*k_3x*x-j*ky*y));
    valtype e_yy = (-A_1*pow(ky,2)*exp(-j*k_1x*x-j*ky*y)-A_2*pow(ky,2)*exp(-j*k_2x*x-j*ky*y)+T3*k_3x*ky*exp(-j*k_3x*x-j*ky*y));
    valtype e_xy = (0.5*(-2.*k_1x*ky*A_1*exp(-j*k_1x*x-j*ky*y)-2.*k_2x*ky*A_2*exp(-j*k_2x*x-j*ky*y)+(pow(k_3x,2)-pow(ky,2))*T3*exp(-j*k_3x*x-j*ky*y)));
    
    // wave 1
    // valtype e_xx = -pow(k_1x,2)*exp(-j*k_1x*x-j*ky*y);
    // valtype e_yy = -pow(ky,2)*exp(-j*k_1x*x-j*ky*y);
    // valtype e_xy = 0.5*(-k_1x*ky*exp(-j*k_1x*x-j*ky*y)-k_1x*ky*exp(-j*k_1x*x-j*ky*y));
    // wave 2
    // valtype e_xx = -pow(k_2x,2)*exp(-j*k_2x*x-j*ky*y);
    // valtype e_yy = -pow(ky,2)*exp(-j*k_2x*x-j*ky*y);
    // valtype e_xy = 0.5*(-k_2x*ky*exp(-j*k_2x*x-j*ky*y)-k_2x*ky*exp(-j*k_2x*x-j*ky*y));
    // wave 3
    // can be normalised by dividing the K_eq_til
    // valtype e_xx = -1e-5*k_3x*ky*exp(-j*k_3x*x-j*ky*y);
    // valtype e_yy = 1e-5*k_3x*ky*exp(-j*k_3x*x-j*ky*y);
    // valtype e_xy = 0.5*(-1e-5*pow(ky,2)*exp(-j*k_3x*x-j*ky*y)+1e-5*pow(k_3x,2)*exp(-j*k_3x*x-j*ky*y));
    
    // u^s
    // sigma_xx
    if (cdgx<0.){
        res(0,0) = 0;
        res(1,1) = 0;
        res(0,1) = 0;
        res(1,0) = res(0,1);
    }
    else
    {
    res(0,0) = 2.0*N*e_xx+A_hat*(e_xx+e_yy);
    res(1,1) = 2.0*N*e_yy+A_hat*(e_xx+e_yy);
    res(0,1) = 2.0*N*e_xy;
    res(1,0) = res(0,1);
    }
    
    return;
}


xEvalExactBiotUtotal::xEvalExactBiotUtotal(double omega_, double beta_, double sigmaF_ , double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_) {
    xEvalCoefficient Coeff(omega, beta, sigmaF, d);
    A_1 = Coeff.A_1;
    A_2 = Coeff.A_2;
    T3 = Coeff.T3;
    k_a = Coeff.k_a;
    delta_1 = Coeff.delta_1;
    delta_2 = Coeff.delta_2, 
    delta_3 = Coeff.delta_3;
    mu_1 =  Coeff.mu_1;
    mu_2 = Coeff.mu_2;
    mu_3 = Coeff.mu_3;
}

void xEvalExactBiotUtotal::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const{

    const double rho_a = 1.213;

    double ky = k_a*sin(beta);
    double k_ax = k_a*cos(beta);
    valtype k_1x = sqrt(pow(delta_1,2)-pow(ky,2));
    valtype k_2x = sqrt(pow(delta_2,2)-pow(ky,2));
    valtype k_3x = sqrt(pow(delta_3,2)-pow(ky,2));

    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    // u^t
    // res(0) = -j*mu_1*k_1x*exp(-j*k_1x*x-j*ky*y) -j*mu_2*k_2x*exp(-j*k_2x*x-j*ky*y) -j*mu_3*ky*exp(-j*k_3x*x-j*ky*y);
    // res(1) = -j*mu_1*ky*exp(-j*k_1x*x-j*ky*y) -j*mu_2*ky*exp(-j*k_2x*x-j*ky*y) +j*mu_3*k_3x*exp(-j*k_3x*x-j*ky*y);
    
    // res(0) = -j*mu_1*k_1x*exp(-j*k_1x*x-j*ky*y)-j*mu_2*k_2x*exp(-j*k_2x*x-j*ky*y);
    // res(1) = -j*mu_1*ky*exp(-j*k_1x*x-j*ky*y)-j*mu_2*ky*exp(-j*k_2x*x-j*ky*y);
    
    // res(0) = -j*mu_1*k_1x*exp(-j*k_1x*x-j*ky*y);
    // res(1) = -j*mu_1*ky*exp(-j*k_1x*x-j*ky*y);
    
    // res(0) = -j*mu_2*k_2x*exp(-j*k_2x*x-j*ky*y);
    // res(1) = -j*mu_2*ky*exp(-j*k_2x*x-j*ky*y);
    
    // res(0) = -1e-5*j*mu_3*ky*exp(-j*k_3x*x-j*ky*y);
    // res(1) = 1e-5*j*mu_3*k_3x*exp(-j*k_3x*x-j*ky*y);

// 2D computation numerically   
    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation < 0) cdgx = -1.e-8;
    else if(orientation > 0) cdgx = 1.e-8;

    if (cdgx<0.){
        res(0) = 0.;
        res(1) = 0.;
    }
    else
    {
        res(0) = (-A_1*j*mu_1*k_1x*exp(-j*k_1x*x-j*ky*y) -A_2*j*mu_2*k_2x*exp(-j*k_2x*x-j*ky*y) -T3*j*mu_3*ky*exp(-j*k_3x*x-j*ky*y));
        res(1) = (-A_1*j*mu_1*ky*exp(-j*k_1x*x-j*ky*y) -A_2*j*mu_2*ky*exp(-j*k_2x*x-j*ky*y) +T3*j*mu_3*k_3x*exp(-j*k_3x*x-j*ky*y));
    }
    
    return;
}


xEvalExactBiotPressure::xEvalExactBiotPressure(double omega_, double beta_, double sigmaF_, double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_) {
    xEvalCoefficient Coeff(omega, beta, sigmaF, d);
    R = Coeff.R;
    // R = 0.38556426694682383-j*0.034164487;
    // T1 = 0.0009450689124152338-j*0.0021046;
    // T2 = 0.8106409333712645-j*0.06397479570258742;
    T1 = Coeff.T1;
    T2 = Coeff.T2;
    k_a = Coeff.k_a;
    delta_1 = Coeff.delta_1;
    delta_2 = Coeff.delta_2, 
    delta_3 = Coeff.delta_3;
    double absp = 1-pow(abs(R),2);
    cout<<"Absoprtion coeff:"<<std::setprecision(10)<<absp<<endl;
    cout<<"Refection coeff:"<<R<<endl;
    cout<<"wave number solid:"<<delta_1<<endl;
    cout<<"wave number fluid:"<<delta_2<<endl;
    cout<<"wave number shear:"<<delta_3<<endl;
}

void xEvalExactBiotPressure::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double>& res) const{
  
    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    double ky = k_a*sin(beta);
    double k_ax = k_a*cos(beta);
    valtype k_1x = sqrt(pow(delta_1,2)-pow(ky,2));
    valtype k_2x = sqrt(pow(delta_2,2)-pow(ky,2));
    valtype k_3x = sqrt(pow(delta_3,2)-pow(ky,2));

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);

    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation > 0) cdgx = -1.e-10;
    else if(orientation < 0) cdgx = 1.e-10;
    // p
    // res = K_eq_til*mu_1*pow(delta_1,2)*exp(-j*k_1x*x-j*ky*y)+K_eq_til*pow(delta_2,2)*mu_2*exp(-j*k_2x*x-j*ky*y);


    if (cdgx<0.) {
        res = (exp(-j*k_ax*x-j*ky*y)+R*exp(j*k_ax*x-j*ky*y));
    }
    else
    {
        res = (T1*exp(-j*k_1x*x-j*ky*y) + T2*exp(-j*k_2x*x-j*ky*y));
    }
    
    return;
}


xEvalExactBiotVelocity::xEvalExactBiotVelocity(double omega_, double beta_, double sigmaF_, double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_) {
    xEvalCoefficient Coeff(omega, beta, sigmaF, d);
    R = Coeff.R;
    T1 = Coeff.T1;
    T2 = Coeff.T2;
    T3 = Coeff.T3;
    // R = 0.38556426694682383-j*0.034164487;
    // T1 = 0.0009450689124152338-j*0.0021046;
    // T2 = 0.8106409333712645-j*0.06397479570258742;
    k_a = Coeff.k_a;
    delta_1 = Coeff.delta_1;
    delta_2 = Coeff.delta_2, 
    delta_3 = Coeff.delta_3;
   

}

void xEvalExactBiotVelocity::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const{

   const double rho_a = 1.213;
    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    double ky = k_a*sin(beta);
    double k_ax = k_a*cos(beta);
    valtype k_1x = sqrt(pow(delta_1,2)-pow(ky,2));
    valtype k_2x = sqrt(pow(delta_2,2)-pow(ky,2));
    valtype k_3x = sqrt(pow(delta_3,2)-pow(ky,2));

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);

    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation > 0) cdgx = -1.e-8;
    else if(orientation < 0) cdgx = 1.e-8;
    // p
    // res = K_eq_til*mu_1*pow(delta_1,2)*exp(-j*k_1x*x-j*ky*y)+K_eq_til*pow(delta_2,2)*mu_2*exp(-j*k_2x*x-j*ky*y);
    // res = K_eq_til*mu_1*pow(delta_1,2)*exp(-j*k_1x*x-j*ky*y);
    // res = K_eq_til*mu_2*pow(delta_2,2)*exp(-j*k_2x*x-j*ky*y);
    // res = K_eq_til*mu_1*pow(delta_1,2)*exp(-j*k_1x*x-j*ky*y);
    // res = K_eq_til*mu_3*k_3x*ky*exp(-j*k_3x*x-j*ky*y)-K_eq_til*mu_3*k_3x*ky*exp(-j*k_3x*x-j*ky*y);

    if (cdgx<0.) {
        // res(0) = 1./(rho_a*pow(omega,2))*(-j*k_ax*exp(-j*k_ax*x-j*ky*y)+j*k_ax*R*exp(j*k_ax*x-j*ky*y));
        // res(1) = 1./(rho_a*pow(omega,2))*(-j*ky*exp(-j*k_ax*x-j*ky*y)-j*ky*R*exp(j*k_ax*x-j*ky*y));
        // To avoid the bad identification of the material in the part of code 
        // to the imposition of the natural velocity BCs, apply the material property mannuelly
        
        res(0) = (-j*k_ax*exp(-j*k_ax*x-j*ky*y)+j*k_ax*R*exp(j*k_ax*x-j*ky*y));
        res(1) = (-j*ky*exp(-j*k_ax*x-j*ky*y)-j*ky*R*exp(j*k_ax*x-j*ky*y));
    }
    else
    {
        res(0) = 0;
        res(1) = 0;
    }
    return;
}


