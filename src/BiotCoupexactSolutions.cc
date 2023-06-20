#include "BiotCoupexactSolutions.h"
#include "xMesh.h"
#include <iomanip>

#define CHECK_WITH_EIGEN 1
#ifdef CHECK_WITH_EIGEN
#include "Eigen/Dense"
#endif

using namespace std;
using namespace xfem;
using namespace Eigen;


PEM_Biot::PEM_Biot(double omega_, double beta_, double phi_, double sigma_, double alpha_, double lambda_prime_, double lambda_, 
double rho_1_, double nu_, double E_, double eta_) : omega(omega_), beta(beta_),
phi(phi_), sigma(sigma_), alpha(alpha_), lambda_prime(lambda_prime_),
lambda(lambda_), rho_1(rho_1_), nu(nu_), E(E_), eta(eta_)

{
  const double P = 101325., gamma = 1.4, rho_a = 1.213, mu = 0.1839e-4, Pr = 0.710;

  double K_a = gamma * P;
  double c_a = sqrt(K_a/rho_a);
  k_a = omega/c_a;
  double Z_a = rho_a * c_a;
  double nu_p = mu / rho_a;
  double nu_prime = nu_p / Pr;

  valtype omega_0 = sigma * phi / (rho_a * alpha);
  valtype omega_inf = pow(sigma, 2)*pow(phi,2)*pow(lambda,2)/(4.*mu*rho_a*pow(alpha,2));
  rho_eq_til  = ((rho_a*alpha/phi)*(1.+((omega_0)/(j*omega))*sqrt(1.+j*omega/omega_inf)));

  valtype omega_prime_infty = (16.*nu_prime)/(pow(lambda_prime, 2));
  valtype F_prime_CA = sqrt(1.+j*omega/omega_prime_infty);
  valtype alpha_prime_til = 1.+omega_prime_infty*F_prime_CA/(2.*j*omega);
  K_eq_til = ((gamma*P/phi)/(gamma-(gamma-1.)/alpha_prime_til));
  valtype c_eq = sqrt(K_eq_til/rho_eq_til);
  valtype k_eq_til = omega / c_eq;
//   cout<<"Wave number"<<k_eq_til<<endl;

  valtype loss = 1. + j*eta;
//   valtype loss = 1.;
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

  N = E*loss/(2.*(1.+nu));
  A_hat = nu*E*loss/((1.+nu)*(1.-2.*nu));
  P_hat = A_hat + 2.*N;

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
    // cout<<"wave number is: "<<delta_1<<delta_2<<delta_3<<endl;

  mu_1 = gamma_til*pow(delta_eq,2)/(pow(delta_1,2)-pow(delta_eq,2));
  mu_2 = gamma_til*pow(delta_eq,2)/(pow(delta_2,2)-pow(delta_eq,2));
  mu_3 = -gamma_til;


}


xEvalCoupCoefficient::xEvalCoupCoefficient(double omega_, double beta_, double sigmaF_, double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_)
{
     // reference parameters Allard page 270 foam 3/4
    double phi_1 = 0.98, sigma_1 = 3.75e3, alpha_1 = 1.17, lambda_prime_1 = 742e-6, lambda_1 = 110e-6, rho_11 = 22.1, nu_1 = 0.39, E_1 = 70e3, eta_1 = 0.265;
    double phi_2 = 0.98, sigma_2 = 15.5e3, alpha_2 = 1.01, lambda_prime_2 = 200e-6, lambda_2 = 100e-6, rho_12 = 11, nu_2 = 0.35, E_2 = 200e3, eta_2 = 0.1;

    // double phi_1 = 0.98, sigma_1 = 1.55e3, alpha_1 = 1.01, lambda_prime_1 = 200e-6, lambda_1 = 100e-6, rho_11 = 11, nu_1 = 0.35, E_1 = 200e3, eta_1 = 0.1;
    // double phi_2 = 0.98, sigma_2 = 1.55e3, alpha_2 = 1.01, lambda_prime_2 = 200e-6, lambda_2 = 100e-6, rho_12 = 11, nu_2 = 0.35, E_2 = 150e3, eta_2 = 0.1;
    PEM_Biot mat1(omega, beta, phi_1, sigma_1, alpha_1, lambda_prime_1, lambda_1, rho_11, nu_1, E_1, eta_1);
    PEM_Biot mat2(omega, beta, phi_2, sigma_2, alpha_2, lambda_prime_2, lambda_2, rho_12, nu_2, E_2, eta_2);

    ky = (mat1.delta_1*sin(beta)).real();
    kx_11 = sqrt(pow(mat1.delta_1,2)-pow(ky,2));
    kx_12 = sqrt(pow(mat1.delta_2,2)-pow(ky,2));
    kx_13 = sqrt(pow(mat1.delta_3,2)-pow(ky,2));

    kx_21 = sqrt(pow(mat2.delta_1,2)-pow(ky,2));
    kx_22 = sqrt(pow(mat2.delta_2,2)-pow(ky,2));
    kx_23 = sqrt(pow(mat2.delta_3,2)-pow(ky,2));

    MatrixXcd M(6,6);
    VectorXcd b(6);
    // pressure is continue
    M(0,0) = 1.;
    M(0,1) = 1.;
    M(0,2) = 0.;
    M(0,3) = -1.-1*omega*sigmaF*d*mat2.mu_1*kx_21/(mat2.K_eq_til*mat2.mu_1*pow(mat2.delta_1,2));
    M(0,4) = -1.-1*omega*sigmaF*d*mat2.mu_2*kx_22/(mat2.K_eq_til*mat2.mu_2*pow(mat2.delta_2,2));
    M(0,5) = 0.-1.*omega*sigmaF*d*mat2.mu_3*ky;
    // u^s_x
    M(1,0) = kx_11/(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
    M(1,1) = kx_12/(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    M(1,2) = -ky;
    M(1,3) = kx_21/(mat2.K_eq_til*mat2.mu_1*pow(mat2.delta_1,2));
    M(1,4) = kx_22/(mat2.K_eq_til*mat2.mu_2*pow(mat2.delta_2,2));
    M(1,5) = ky;
    // u^s_y
    M(2,0) = -ky/(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
    M(2,1) = -ky/(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    M(2,2) = -kx_13;
    M(2,3) = ky/(mat2.K_eq_til*mat2.mu_1*pow(mat2.delta_1,2));
    M(2,4) = ky/(mat2.K_eq_til*mat2.mu_2*pow(mat2.delta_2,2));
    M(2,5) = -kx_23;
    // u^t_x
    M(3,0) = mat1.mu_1*kx_11/(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
    M(3,1) = mat1.mu_2*kx_12/(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    M(3,2) = -mat1.mu_3*ky;
    M(3,3) = mat2.mu_1*kx_21/(mat2.K_eq_til*mat2.mu_1*pow(mat2.delta_1,2));
    M(3,4) = mat2.mu_2*kx_22/(mat2.K_eq_til*mat2.mu_2*pow(mat2.delta_2,2));
    M(3,5) = mat2.mu_3*ky;
    // sigma_xx
    M(4,0) = (2.*mat1.N*pow(kx_11,2)+mat1.A_hat*pow(mat1.delta_1,2))/(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
    M(4,1) = (2.*mat1.N*pow(kx_12,2)+mat1.A_hat*pow(mat1.delta_2,2))/(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    M(4,2) = -2.*mat1.N*kx_13*ky;
    M(4,3) = -(2.*mat2.N*pow(kx_21,2)+mat2.A_hat*pow(mat2.delta_1,2))/(mat2.K_eq_til*mat2.mu_1*pow(mat2.delta_1,2));
    M(4,4) = -(2.*mat2.N*pow(kx_22,2)+mat2.A_hat*pow(mat2.delta_2,2))/(mat2.K_eq_til*mat2.mu_2*pow(mat2.delta_2,2));
    M(4,5) = -2.*mat2.N*kx_23*ky;
    // sigma_xy
    M(5,0) = 2.*mat1.N*kx_11*ky/(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
    M(5,1) = 2.*mat1.N*kx_12*ky/(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    M(5,2) = mat1.N*(pow(kx_13,2)-pow(ky,2));
    M(5,3) = 2.*mat2.N*kx_21*ky/(mat2.K_eq_til*mat2.mu_1*pow(mat2.delta_1,2));
    M(5,4) = 2.*mat2.N*kx_22*ky/(mat2.K_eq_til*mat2.mu_2*pow(mat2.delta_2,2));
    M(5,5) = -mat2.N*(pow(kx_23,2)-pow(ky,2));

// first compressionel wave
    b(0) = -1.;
    b(1) = kx_11/(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
    b(2) = ky/(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
    b(3) = (mat1.mu_1*kx_11)/(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
    b(4) = -(2.*mat1.N*pow(kx_11,2)+mat1.A_hat*pow(mat1.delta_1,2))/(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
    b(5) = 2.*mat1.N*kx_11*ky/(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
// second compresseionel wave
    // b(0) = -1.;
    // b(1) = kx_12/(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    // b(2) = ky/(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    // b(3) = (mat1.mu_2*kx_12)/(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    // b(4) = -(2.*mat1.N*pow(kx_12,2)+mat1.A_hat*pow(mat1.delta_2,2))/(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    // b(5) = 2.*mat1.N*kx_12*ky/(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    VectorXcd c = M.partialPivLu().solve(b);
    R1 = c(0);
    R2 = c(1);
    R3 = c(2);
    T1 = c(3);
    T2 = c(4);
    T3 = c(5);
// 0 1256
    // R1 = 0.08628202625737383-0.15217210650364268*j;
    // R2 = 0.5572790539722784+0.17876796052964372*j;
    // R3 = 0.*j;
    // T1 = -0.3795697022179179-1.1638033597995514*j;
    // T2 = 2.0632696455606956+1.2960996984136206*j;
    // T3 = -0.*j;
// // 45 1256
    // R1 = -0.3964746993705698-0.09355957327578461*j;
    // R2 = 0.714335671843966+0.8003087616160831*j;
    // R3 = 3.015674397130168e-08+8.572857999076112e-09*j;
    // T1 = 0.25496881327705895-0.5495977518571338*j;
    // T2 = 1.5297195652703555+1.3721969310997082*j;
    // T3 = -1.7417621151033124e-08-8.421803150728344e-08*j;
// 0 12560
    // R1 = 0.02820616351614082-0.05965961019496064*j;
    // R2 = 0.2523577881430488-0.22946175085065354*j;
    // R3 = 0.*j;
    // T1 = 2.1943989451927224-3.1960163409640274*j;
    // T2 = -0.9238525807849587+2.7495593942101926*j;
    // T3 = -0.*j;

// first compressionel wave
    A0 = 1./(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
    A11 = R1/(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
    A12 = R2/(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    A21 = T1/(mat2.K_eq_til*mat2.mu_1*pow(mat2.delta_1,2));
    A22 = T2/(mat2.K_eq_til*mat2.mu_2*pow(mat2.delta_2,2));
// second compresseionel wave
    // A0 = 1./(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    // A11 = c(0)/(mat1.K_eq_til*mat1.mu_1*pow(mat1.delta_1,2));
    // A12 = c(1)/(mat1.K_eq_til*mat1.mu_2*pow(mat1.delta_2,2));
    // A21 = c(3)/(mat2.K_eq_til*mat2.mu_1*pow(mat2.delta_1,2));
    // A22 = c(4)/(mat2.K_eq_til*mat2.mu_2*pow(mat2.delta_2,2));

    mat1_N = mat1.N;
    mat1_A_hat = mat1.A_hat;
    mat1_delta_1 = mat1.delta_1;
    mat1_delta_2 = mat1.delta_2;
    mat1_delta_3 = mat1.delta_3;
    mat1_mu_1 = mat1.mu_1;
    mat1_mu_2 = mat1.mu_2;
    mat1_mu_3 = mat1.mu_3;

    // cout<<"wave number 11: "<<mat1.delta_1<<endl;
    // cout<<"wave number 12: "<<mat1.delta_2<<endl;
    // cout<<"wave number 13: "<<mat1.delta_3<<endl;

    // cout<<"ratio velocity 1: "<<mat1.mu_1<<endl;
    // cout<<"ratio velocity 2: "<<mat1.mu_2<<endl;
    // cout<<"ratio velocity 3: "<<mat1.mu_3<<endl;

    mat2_N = mat2.N;
    mat2_A_hat = mat2.A_hat;
    mat2_delta_1 = mat2.delta_1;
    mat2_delta_2 = mat2.delta_2;
    mat2_delta_3 = mat2.delta_3;
    mat2_mu_1 = mat2.mu_1;
    mat2_mu_2 = mat2.mu_2;
    mat2_mu_3 = mat2.mu_3;    
    // cout<<"wave number 21: "<<mat2.delta_1<<endl;
    // cout<<"wave number 22: "<<mat2.delta_2<<endl;
    // cout<<"wave number 23: "<<mat2.delta_3<<endl;
    mat1_rho_eq_til = mat1.rho_eq_til;
    mat1_K_eq_til = mat1.K_eq_til;
    mat2_rho_eq_til = mat2.rho_eq_til;
    mat2_K_eq_til = mat2.K_eq_til;


}


xEvalExactCoupBiotUs::xEvalExactCoupBiotUs(double omega_, double beta_, double sigmaF_, double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_){
    xEvalCoupCoefficient Coeff(omega, beta, sigmaF, d);
    A0 = Coeff.A0;
    A11 = Coeff.A11;
    A12 = Coeff.A12;
    A21 = Coeff.A21;
    A22 = Coeff.A22;
    R3 = Coeff.R3;
    T3 = Coeff.T3;
    ky = Coeff.ky;
    kx_11 = Coeff.kx_11;
    kx_12 = Coeff.kx_12;
    kx_13 = Coeff.kx_13;
    kx_21 = Coeff.kx_21;
    kx_22 = Coeff.kx_22;
    kx_23 = Coeff.kx_23;

}

void xEvalExactCoupBiotUs::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const{

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
    if(orientation < 0) cdgx = -1.e-8;
    else if(orientation > 0) cdgx = 1.e-8;
    if (cdgx<0.){
        res(0) = 10.e-2*(-j*kx_11*A0*exp(-j*kx_11*x-j*ky*y)+j*kx_11*A11*exp(j*kx_11*x-j*ky*y)
                +j*kx_12*A12*exp(j*kx_12*x-j*ky*y)-j*ky*R3*exp(j*kx_13*x-j*ky*y));
        res(1) = 10.e-2*(-j*ky*A0*exp(-j*kx_11*x-j*ky*y)-j*ky*A11*exp(j*kx_11*x-j*ky*y)
                -j*ky*A12*exp(j*kx_12*x-j*ky*y)-j*kx_13*R3*exp(j*kx_13*x-j*ky*y));

        // res(0) = 10.e3*(-j*kx_12*A0*exp(-j*kx_12*x-j*ky*y)+j*kx_11*A11*exp(j*kx_11*x-j*ky*y)
        //         +j*kx_12*A12*exp(j*kx_12*x-j*ky*y)-j*ky*R3*exp(j*kx_13*x-j*ky*y));
        // res(1) = 10.e3*(-j*ky*A0*exp(-j*kx_12*x-j*ky*y)-j*ky*A11*exp(j*kx_11*x-j*ky*y)
        //         -j*ky*A12*exp(j*kx_12*x-j*ky*y)-j*kx_13*R3*exp(j*kx_13*x-j*ky*y));
    }
    else
    {
        res(0) = 10.e-2*(-j*kx_21*A21*exp(-j*kx_21*x-j*ky*y)-j*kx_22*A22*exp(-j*kx_22*x-j*ky*y)-j*ky*T3*exp(-j*kx_23*x-j*ky*y));
        res(1) = 10.e-2*(-j*ky*A21*exp(-j*kx_21*x-j*ky*y)-j*ky*A22*exp(-j*kx_22*x-j*ky*y) +j*kx_23*T3*exp(-j*kx_23*x-j*ky*y));


    }
    
    return;
}


xEvalExactCoupBiotTraction::xEvalExactCoupBiotTraction(double omega_, double beta_, double sigmaF_, double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_){
    xEvalCoupCoefficient Coeff(omega, beta, sigmaF, d);
    A0 = Coeff.A0;
    A11 = Coeff.A11;
    A12 = Coeff.A12;
    R3 = Coeff.R3;
    A21 = Coeff.A21;
    A22 = Coeff.A22;
    T3 = Coeff.T3;
    ky = Coeff.ky;
    kx_11 = Coeff.kx_11;
    kx_12 = Coeff.kx_12;
    kx_13 = Coeff.kx_13;
    kx_21 = Coeff.kx_21;
    kx_22 = Coeff.kx_22;
    kx_23 = Coeff.kx_23;
    mat1_N = Coeff.mat1_N;
    mat1_A_hat =Coeff.mat1_A_hat;
    mat2_N = Coeff.mat2_N;
    mat2_A_hat = Coeff.mat2_A_hat;
}

void xEvalExactCoupBiotTraction::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xTensor2<std::complex<double> >& res) const{


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

    if(orientation < 0) cdgx = -1.e-8;
    else if(orientation > 0) cdgx = 1.e-8;
 
    valtype e_xx1 = 10.e-2*(-A0*pow(kx_11,2)*exp(-j*kx_11*x-j*ky*y)-A11*pow(kx_11,2)*exp(j*kx_11*x-j*ky*y)-A12*pow(kx_12,2)*exp(j*kx_12*x-j*ky*y)+R3*kx_13*ky*exp(j*kx_13*x-j*ky*y));
    valtype e_yy1 = 10.e-2*(-A0*pow(ky,2)*exp(-j*kx_11*x-j*ky*y)-A11*pow(ky,2)*exp(j*kx_11*x-j*ky*y)-A12*pow(ky,2)*exp(j*kx_12*x-j*ky*y)-R3*kx_13*ky*exp(j*kx_13*x-j*ky*y));
    valtype e_xy1 = 10.e-2*0.5*(-2.*kx_11*ky*A0*exp(-j*kx_11*x-j*ky*y)+2.*kx_11*ky*A11*exp(j*kx_11*x-j*ky*y)+2.*kx_12*ky*A12*exp(j*kx_12*x-j*ky*y)+(pow(kx_13,2)-pow(ky,2))*R3*exp(j*kx_13*x-j*ky*y));

    // valtype e_xx1 = 10.e3*(-A0*pow(kx_12,2)*exp(-j*kx_12*x-j*ky*y)-A11*pow(kx_11,2)*exp(j*kx_11*x-j*ky*y)-A12*pow(kx_12,2)*exp(j*kx_12*x-j*ky*y)+R3*kx_13*ky*exp(j*kx_13*x-j*ky*y));
    // valtype e_yy1 = 10.e3*(-A0*pow(ky,2)*exp(-j*kx_12*x-j*ky*y)-A11*pow(ky,2)*exp(j*kx_11*x-j*ky*y)-A12*pow(ky,2)*exp(j*kx_12*x-j*ky*y)-R3*kx_13*ky*exp(j*kx_13*x-j*ky*y));
    // valtype e_xy1 = 10.e3*0.5*(-2.*kx_12*ky*A0*exp(-j*kx_12*x-j*ky*y)+2.*kx_11*ky*A11*exp(j*kx_11*x-j*ky*y)+2.*kx_12*ky*A12*exp(j*kx_12*x-j*ky*y)+(pow(kx_13,2)-pow(ky,2))*R3*exp(j*kx_13*x-j*ky*y));
    
    valtype e_xx2 = 10.e-2*(-A21*pow(kx_21,2)*exp(-j*kx_21*x-j*ky*y)-A22*pow(kx_22,2)*exp(-j*kx_22*x-j*ky*y)-T3*kx_23*ky*exp(-j*kx_23*x-j*ky*y));
    valtype e_yy2 = 10.e-2*(-A21*pow(ky,2)*exp(-j*kx_21*x-j*ky*y)-A22*pow(ky,2)*exp(-j*kx_22*x-j*ky*y)+T3*kx_23*ky*exp(-j*kx_23*x-j*ky*y));
    valtype e_xy2 = 10.e-2*0.5*(-2.*kx_21*ky*A21*exp(-j*kx_21*x-j*ky*y)-2.*kx_22*ky*A22*exp(-j*kx_22*x-j*ky*y)+(pow(kx_23,2)-pow(ky,2))*T3*exp(-j*kx_23*x-j*ky*y));
    
    
    // u^s
    // sigma_xx
    if (cdgx<0.){
    res(0,0) = 2.0*mat1_N*e_xx1+mat1_A_hat*(e_xx1+e_yy1);
    res(1,1) = 2.0*mat1_N*e_yy1+mat1_A_hat*(e_xx1+e_yy1);
    res(0,1) = 2.0*mat1_N*e_xy1;
    res(1,0) = res(0,1);
    }
    else
    {
    res(0,0) = 2.0*mat2_N*e_xx2+mat2_A_hat*(e_xx2+e_yy2);
    res(1,1) = 2.0*mat2_N*e_yy2+mat2_A_hat*(e_xx2+e_yy2);
    res(0,1) = 2.0*mat2_N*e_xy2;
    res(1,0) = res(0,1);
    }
    
    return;
}


xEvalExactCoupBiotUtotal::xEvalExactCoupBiotUtotal(double omega_, double beta_, double sigmaF_, double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_){
    xEvalCoupCoefficient Coeff(omega, beta, sigmaF, d);
    A0 = Coeff.A0;
    A11 = Coeff.A11;
    A12 = Coeff.A12;
    A21 = Coeff.A21;
    A22 = Coeff.A22;
    R3 = Coeff.R3;
    T3 = Coeff.T3;

    ky = Coeff.ky;
    kx_11 = Coeff.kx_11;
    kx_12 = Coeff.kx_12;
    kx_13 = Coeff.kx_13;
    kx_21 = Coeff.kx_21;
    kx_22 = Coeff.kx_22;
    kx_23 = Coeff.kx_23;
    mat1_mu_1 =  Coeff.mat1_mu_1;
    mat1_mu_2 = Coeff.mat1_mu_2;
    mat1_mu_3 = Coeff.mat1_mu_3;
    mat2_mu_1 =  Coeff.mat2_mu_1;
    mat2_mu_2 = Coeff.mat2_mu_2;
    mat2_mu_3 = Coeff.mat2_mu_3;
}

void xEvalExactCoupBiotUtotal::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, xtensor::xVector<std::complex<double> >& res) const{

    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);
    // u^t

// 2D computation numerically   
    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation < 0) cdgx = -1.e-8;
    else if(orientation > 0) cdgx = 1.e-8;

    if (cdgx<0.){
        res(0) = 10.e-2*(-mat1_mu_1*j*kx_11*A0*exp(-j*kx_11*x-j*ky*y)+mat1_mu_1*j*kx_11*A11*exp(j*kx_11*x-j*ky*y)
                +mat1_mu_2*j*kx_12*A12*exp(j*kx_12*x-j*ky*y)-mat1_mu_3*j*ky*R3*exp(j*kx_13*x-j*ky*y));
        res(1) = 10.e-2*(-mat1_mu_1*j*ky*A0*exp(-j*kx_11*x-j*ky*y)-mat1_mu_1*j*ky*A11*exp(j*kx_11*x-j*ky*y)
                -mat1_mu_2*j*ky*A12*exp(j*kx_12*x-j*ky*y)-mat1_mu_3*j*kx_13*R3*exp(j*kx_13*x-j*ky*y));

        // res(0) = 10.e3*(-mat1_mu_2*j*kx_12*A0*exp(-j*kx_12*x-j*ky*y)+mat1_mu_1*j*kx_11*A11*exp(j*kx_11*x-j*ky*y)
        //         +mat1_mu_2*j*kx_12*A12*exp(j*kx_12*x-j*ky*y)-mat1_mu_3*j*ky*R3*exp(j*kx_13*x-j*ky*y));
        // res(1) = 10.e3*(-mat1_mu_2*j*ky*A0*exp(-j*kx_12*x-j*ky*y)-mat1_mu_1*j*ky*A11*exp(j*kx_11*x-j*ky*y)
        //         -mat1_mu_2*j*ky*A12*exp(j*kx_12*x-j*ky*y)-mat1_mu_3*j*kx_13*R3*exp(j*kx_13*x-j*ky*y));
    }
    else
    {
        res(0) = 10.e-2*(-mat2_mu_1*j*kx_21*A21*exp(-j*kx_21*x-j*ky*y)-mat2_mu_2*j*kx_22*A22*exp(-j*kx_22*x-j*ky*y)-mat2_mu_3*j*ky*T3*exp(-j*kx_23*x-j*ky*y));
        res(1) = 10.e-2*(-mat2_mu_1*j*ky*A21*exp(-j*kx_21*x-j*ky*y)-mat2_mu_2*j*ky*A22*exp(-j*kx_22*x-j*ky*y) + mat2_mu_3*j*kx_23*T3*exp(-j*kx_23*x-j*ky*y));
    }
    
    return;
}


xEvalExactCoupBiotPressure::xEvalExactCoupBiotPressure(double omega_, double beta_, double sigmaF_, double d_) : omega(omega_), beta(beta_), sigmaF(sigmaF_), d(d_){
    xEvalCoupCoefficient Coeff(omega, beta, sigmaF, d);
    R1 = Coeff.R1;
    R2 = Coeff.R2;
    R3 = Coeff.R3;
    T1 = Coeff.T1;
    T2 = Coeff.T2;
    T3 = Coeff.T3;

    cout<<"Reflection 1: "<<R1<<endl;
    cout<<"Reflection 2: "<<R2<<endl;
    cout<<"Reflection 3: "<<R3<<endl;
    cout<<"Transmission 1: "<<T1<<endl;
    cout<<"Transmission 2: "<<T2<<endl;
    cout<<"Transmission 3: "<<T3<<endl;
    ky = Coeff.ky;
    kx_11 = Coeff.kx_11;
    kx_12 = Coeff.kx_12;
    kx_13 = Coeff.kx_13;
    kx_21 = Coeff.kx_21;
    kx_22 = Coeff.kx_22;
    kx_23 = Coeff.kx_23;

}

void xEvalExactCoupBiotPressure::operator()(const xGeomElem*  geo_appro, const xGeomElem* geo_integ, std::complex<double>& res) const{
  
    // Trellis_Util::mPoint xyz = geo_appro->getXYZ();
    xtensor::xPoint xyz = geo_appro->getXYZ();
    double x = xyz(0);
    double y = xyz(1);

    auto cdgxyz =geo_integ->getCDGxyz();
    double cdgx = cdgxyz(0);

    unsigned int side_tag;
    side_tag= AOMD::AOMD_Util::Instance()->lookupMeshDataId("side_tag");
    int orientation=(geo_integ->getEntity())->getAttachedInt(side_tag);

    if(orientation < 0) cdgx = -1.e-8;
    else if(orientation > 0) cdgx = 1.e-8;
    // p
    // res = K_eq_til*mu_1*pow(delta_1,2)*exp(-j*k_1x*x-j*ky*y)+K_eq_til*pow(delta_2,2)*mu_2*exp(-j*k_2x*x-j*ky*y);


    if (cdgx<0.) {
        res = 10.e-2*(exp(-j*kx_11*x-j*ky*y)+R1*exp(j*kx_11*x-j*ky*y)+R2*exp(j*kx_12*x-j*ky*y));
        // res = 10.e3*(exp(-j*kx_12*x-j*ky*y)+R1*exp(j*kx_11*x-j*ky*y)+R2*exp(j*kx_12*x-j*ky*y));
    }
    else
    {
        res = 10.e-2*(T1*exp(-j*kx_21*x-j*ky*y) + T2*exp(-j*kx_22*x-j*ky*y));
    }
    
    return;
}



