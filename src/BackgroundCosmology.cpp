#include <iostream>
#include <cmath>
#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================

BackgroundCosmology::BackgroundCosmology(
    double h,
    double OmegaB,
    double OmegaCDM,
    double OmegaK,
    double Neff,
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff),
  TCMB(TCMB)
{

  OmegaR = 2*pow(M_PI,2)/30 * pow(Constants.k_b*TCMB,4)/ (pow(Constants.hbar,3)*pow(Constants.c,5))*8*M_PI*Constants.G/(3*pow(Constants.H0_over_h,2));
  OmegaLambda = 1 - (OmegaK+OmegaB+OmegaCDM+OmegaR+OmegaNu);
  OmegaNu = 0;
  H0 = Constants.H0_over_h * h;
  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaLambda, H0, ...
  //=============================================================================
  //...
  //...
  //...
  //...
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  const double xmin = -20;
  const double xmax = 2.0;
  const int    npts = 1000;
  Vector x_array = Utils::linspace(xmin, xmax, npts);
  //nt x = 0;

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //double detadx = Constants.c/Hp

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    //...
    //...

    detadx[0] = Constants.c/Hp_of_x(x);

    return GSL_SUCCESS;
  };
  double yini = Constants.c/Hp_of_x(0);
  Vector y_ic{yini};

  ODESolver ode;
  ode.solve(detadx, x_array, y_ic);

  auto y_array = ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, y_array, "$\eta(x)$");

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline
  //=============================================================================
  // ...
  // ...
  // ...
  // ...

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{


  double Omega_value1 = ((OmegaB+OmegaCDM)*exp(-3*x) + (OmegaR)*exp(-4*x) + OmegaLambda);
  //double Omega_value2 = -3 * (OmegaB+OmegaCDM)*exp(-3*x) - 4 * (OmegaR+OmegaNu)*exp(-4*x);



  double H = H0 * sqrt(Omega_value1);

  /*double dHdx = H0*0.5 * 1 / sqrt(Omega_value1) * Omega_value2;

  double ddHddx = H0 * 0.5 * (-0.5 * 1/pow(Omega_value1,-3/2))*Omega_value2 + H0 * 0.5 * 1/sqrt(Omega_value1) * (9 * (OmegaB0+OmegaCDM)*exp(-3*x)) * 16 * OmegaR * exp(-4*x);
  */

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return H;
}


double BackgroundCosmology::dHdx_of_x(double x) const{


  double Omega_value1 = ((OmegaB+OmegaCDM)*exp(-3*x) + (OmegaR)*exp(-4*x) + OmegaLambda);
  double Omega_value2 = (-3 * (OmegaB+OmegaCDM)*exp(-3*x) - 4 * (OmegaR)*exp(-4*x));



  double dHdx = 0.5 * H0 * sqrt(Omega_value1) * Omega_value2;

  /*double dHdx = H0*0.5 * 1 / sqrt(Omega_value1) * Omega_value2;

  double ddHddx = H0 * 0.5 * (-0.5 * 1/pow(Omega_value1,-3/2))*Omega_value2 + H0 * 0.5 * 1/sqrt(Omega_value1) * (9 * (OmegaB0+OmegaCDM)*exp(-3*x)) * 16 * OmegaR * exp(-4*x);
  */

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return dHdx;
}

double BackgroundCosmology::ddHddx_of_x(double x) const{


  double Omega_value1 = ((OmegaB+OmegaCDM)*exp(-3*x) + (OmegaR)*exp(-4*x) + OmegaLambda);
  double Omega_value2 = (-3 * (OmegaB+OmegaCDM)*exp(-3*x) - 4 * (OmegaR)*exp(-4*x));
  double Omega_value3 = (9*(Omega_value3+OmegaCDM)*exp(-3*x) + 16*OmegaR*exp(-4*x));



  double ddHddx = 0.5 * H0 * (0.5 * pow(Omega_value1,(-3/2))) * Omega_value2 + 0.5*H0*1/sqrt(Omega_value1) * Omega_value3;

  /*double dHdx = H0*0.5 * 1 / sqrt(Omega_value1) * Omega_value2;

  double ddHddx = H0 * 0.5 * (-0.5 * 1/pow(Omega_value1,-3/2))*Omega_value2 + H0 * 0.5 * 1/sqrt(Omega_value1) * (9 * (OmegaB0+OmegaCDM)*exp(-3*x)) * 16 * OmegaR * exp(-4*x);
  */

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return ddHddx;
}

double BackgroundCosmology::Hp_of_x(double x) const{

  double Omega_value_Hp1 = ((OmegaB+OmegaCDM)*exp(-x) + OmegaR*exp(-4*x) + OmegaLambda);

  double Hp = H0 * sqrt(Omega_value_Hp1);

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  double Omega_value_Hp1 = ((OmegaB+OmegaCDM)*exp(-x) + OmegaR*exp(-4*x) + OmegaLambda);
  double Omega_value_Hp2 = (-(OmegaB+OmegaCDM)*exp(-x) -2 * OmegaR*exp(-2*x) + OmegaLambda*exp(-2*x));



  double dHpdx = 0.5 * H0 * 1/sqrt(Omega_value_Hp1) * Omega_value_Hp2;



  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{


  double Omega_value_Hp1 = ((OmegaB+OmegaCDM)*exp(-x) + OmegaR*exp(-4*x) + OmegaLambda);
  double Omega_value_Hp2 = (-(OmegaB+OmegaCDM)*exp(-x) -2 * OmegaR*exp(-2*x) + OmegaLambda*exp(-2*x));
  double Omega_value_Hp3 = ( -(OmegaB+OmegaCDM)*exp(-x) + 4*OmegaR + 4*OmegaLambda*exp(-2*x)  );



  double ddHpddx = 0.5*H0 * (-0.5*pow(Omega_value_Hp1,(-3/2))) * Omega_value_Hp2 + H0*0.5/sqrt(Omega_value_Hp1) * Omega_value_Hp3;


  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return ddHpddx;
}

double BackgroundCosmology::get_OmegaB(double x) const{
  if(x == 0.0) return OmegaB;


  return OmegaB / (pow(exp(x),3)*(pow(H_of_x(x),2))/pow(H0,2));
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  //return OmegaB;
}

double BackgroundCosmology::get_OmegaR(double x) const{
  if(x == 0.0) return OmegaR;


  return OmegaB/ (exp(x)*(pow(H_of_x(x),2))/pow(H0,2));
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

}

double BackgroundCosmology::get_OmegaNu(double x) const{
  if(x == 0.0) return OmegaNu;

  return OmegaNu / (pow(exp(x),4)*(pow(H_of_x(x),2))/pow(H0,2));

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

}

double BackgroundCosmology::get_OmegaCDM(double x) const{
  if(x == 0.0) return OmegaCDM;

  return OmegaCDM / (pow(exp(x),3)*(pow(H_of_x(x),2))/pow(H0,2));
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

}

double BackgroundCosmology::get_OmegaLambda(double x) const{
  if(x == 0.0) return OmegaLambda;


  return OmegaLambda / (pow(exp(x),2)*(pow(H_of_x(x),2))/pow(H0,2));
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  //return OmegaLambda;
}

double BackgroundCosmology::get_OmegaK(double x) const{
  if(x == 0.0) return OmegaK;

  return OmegaK / (pow(exp(x),2)*(pow(H_of_x(x),2))/pow(H0,2));

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{
  return H0;
}

double BackgroundCosmology::get_h() const{
  return h;
}

double BackgroundCosmology::get_Neff() const{
  return Neff;
}

double BackgroundCosmology::get_TCMB(double x) const{
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x);
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{
  std::cout << OmegaLambda << "\n";
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -20.0;
  const double x_max =  2.0;
  const int    n_pts =  1000;

  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)/Constants.Mpc        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
