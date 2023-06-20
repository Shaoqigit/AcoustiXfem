#include "xMaterial.h"
#include "xParseData.h"

using namespace xfem;

namespace xfem {


class xAcousticLiquid : virtual public xMaterial {

public:
  xAcousticLiquid();
//Definition of the pure virtual functions
  void checkProperties() override;

  void  sensitivityTo(const std::string& phys_token, std::complex<double> &sensitivity) override;

protected:
  // Useful tensors derived from the above
  // set up in the constructor
  double fluid_density;
  double sound_celerity;
  double Z_a;
};


class xAcousticAir : virtual public xMaterial {

public:
using valtype = std::complex<double>;
  xAcousticAir();
//Definition of the pure virtual functions
  void checkProperties() override;

  void  sensitivityTo(const std::string& phys_token, std::complex<double> &sensitivity) override;
  void  sensitivityTo(const std::string& phys_token, xtensor::xTensor4<std::complex<double>> &sensitivity) override;

protected:
  // Useful tensors derived from the above
  // set up in the constructor
  double rho_a;
  double c_a;
  double K_a, gamma, P, mu, nu, nu_prime, Z_a;
  double fluid_density;
  double sound_celerity;
};


class xAcousticEquivalentFluid : virtual public xMaterial {

public:
  xAcousticEquivalentFluid();
//Definition of the pure virtual functions
  void checkProperties() override;
  void  sensitivityTo(const std::string& phys_token, std::complex<double> &sensitivity) override;
//  void  sensitivityTo(const std::string& phys_token, double &sensitivity) override;
//  inline static double delta(int i, int j) {return ((i == j)?1.0:0.0);}
//  static void SetThermicConductivityIsotropic(double k, xtensor::xTensor2<>& conductivity);

  static void setOmega(double freq){
      omega = 2.*M_PI*freq;
//      this->checkProperties();
  }

protected:
  const std::complex<double> j{0.,1.};

  double fluid_density, sound_celerity, d;
  std::complex<double> rho_eq, c_eq;
  std::complex<double> F_prime_CA, alpha_prime_til, K;
  std::complex<double> P, Pr, rho_a, mu, nu, nu_prime, gamma,  phi, sigma, alpha, lambda, lambda_prime;
  std::complex<double> omega_0, omega_inf, omega_prime_infty;

  static double omega;

};


class xAcousticLimpFluid : virtual public xMaterial {

public:
  xAcousticLimpFluid();
//Definition of the pure virtual functions
  void checkProperties() override;
  void  sensitivityTo(const std::string& phys_token, std::complex<double> &sensitivity) override;
//  void  sensitivityTo(const std::string& phys_token, double &sensitivity) override;
//  inline static double delta(int i, int j) {return ((i == j)?1.0:0.0);}
//  static void SetThermicConductivityIsotropic(double k, xtensor::xTensor2<>& conductivity);

  static void setOmega(double freq){
      omega = 2.*M_PI*freq;
//      this->checkProperties();
  }

protected:
  using valtype = std::complex<double>;
  const std::complex<double> j{0.,1.};

  double fluid_density, sound_celerity, d;
  std::complex<double> rho_eq_til, c_eq, rho_limp_til;
  std::complex<double> F_prime_CA, alpha_prime_til, K;
  std::complex<double> P, Pr, rho_a, mu, nu, nu_prime, gamma,  phi, sigma, alpha, lambda, lambda_prime;
  std::complex<double> omega_0, omega_inf, omega_prime_infty;

  double rho_1;  // Mass of solid per unit volume of aggregate
  double nu_p; // poisson ratio
  double E;  // Young's modulus
  double eta;  // loss facotr
  double loss_type;  // type of lossed for the frame (anelastic or structural)
  valtype lambda_, rho_12, rho_11, rho_2, rho_22, gamma_til;
  valtype rho_11_til, rho_12_til, rho_22_til, rho_til, rho_s_til, loss;
  double N, A_hat, P_hat;

  static double omega;

};


class xAcousticPEMBiot : virtual public xMaterial {

public:
    xAcousticPEMBiot();
    void checkProperties() override;
    void  sensitivityTo(const std::string& phys_token, std::complex<double> &sensitivity) override;
    void  sensitivityTo(const std::string& phys_token, xtensor::xTensor4<std::complex<double>> &sensitivity) override;
    // void  sensitivityTo(const std::string& phys_token, xtensor::xTensor4Isotropic &sensitivity)
  static inline double delta(int i, int j) {return ((i == j)?1.0:0.0);}
  static void SetElasticStiffnessIsotropic   (double E, double nu_p, std::complex<double> loss, xtensor::xTensor4< std::complex<double> >&  tensor);

  static void setOmega(double freq){
    omega = 2.*M_PI*freq;
//      this->checkProperties();
  }

protected:
    using valtype = std::complex<double>;
    const valtype j{0.,1.};

    valtype rho_eq_til, c_eq;
    valtype F_prime_CA, alpha_prime_til, K_eq_til;
    valtype P, Pr, rho_a, mu, nu, nu_prime, gamma,  phi, sigma, alpha, lambda, lambda_prime;
    valtype omega_0, omega_inf, omega_prime_infty;

    double rho_1;  // Mass of solid per unit volume of aggregate
    double nu_p; // poisson ratio
    double E;  // Young's modulus
    double eta;  // loss facotr
    double loss_type;  // type of lossed for the frame (anelastic or structural)
    valtype lambda_, rho_12, rho_11, rho_2, rho_22, gamma_til;
    valtype rho_11_til, rho_12_til, rho_22_til, rho_til, rho_s_til, loss;
    double N, A_hat, P_hat;
    xtensor::xTensor4<valtype> elastic_stiffness;

    static double omega;

};


}
