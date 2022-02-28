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
  OmegaB(OmegaB), //BAryonic
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff),
  TCMB(TCMB)
{

  OmegaR = 2*pow(M_PI,2)/30 * pow(Constants.k_b*TCMB,4)/ (pow(Constants.hbar,3)*pow(Constants.c,5))*8*M_PI*Constants.G/(3*pow(Constants.H0_over_h*h,2)); //Radiation
  OmegaLambda = 1 - (OmegaK+OmegaB+OmegaCDM+OmegaR+OmegaNu); //Dark energy
  OmegaNu = 0;
  H0 = Constants.H0_over_h * h;
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
  const double xmax = 3.0;
  const int    npts = 1000;
  Vector x_array = Utils::linspace(xmin, xmax, npts);
  //nt x = 0;

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    detadx[0] = Constants.c/Hp_of_x(x);

    return GSL_SUCCESS;
  };
  double yini = Constants.c/Hp_of_x(-20);
  Vector y_ic{yini};

  ODESolver ode;
  ode.solve(detadx, x_array, y_ic);

  auto y_array = ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, y_array, "$\eta(x)$");


  Utils::EndTiming("Eta");
}

void BackgroundCosmology::solve_t(){
  Utils::StartTiming("t");
  const double xmin = -20;
  const double xmax = 3;
  const int    npts = 1000;
  Vector x_array = Utils::linspace(xmin, xmax, npts);


  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){



    dtdx[0] = 1/H_of_x(x);

    return GSL_SUCCESS;
  };
  double yini = 1.0/(2*H_of_x(x_array[0]));
  Vector y_ic{yini};

  ODESolver ode;
  ode.solve(dtdx, x_array, y_ic);

  auto t_array = ode.get_data_by_component(0);

  t_of_x_spline.create(x_array, t_array, "$t(x)$");


  Utils::EndTiming("t");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{


  double Omega_value1 = ((OmegaB+OmegaCDM)*exp(-3*x) + (OmegaR+OmegaNu)*exp(-4*x) + OmegaLambda + OmegaK*exp(-2*x));



  double H = H0 * sqrt(Omega_value1);


  return H;
}


double BackgroundCosmology::dHdx_of_x(double x) const{


  double Omega_value1 = ((OmegaB+OmegaCDM)*exp(-3*x) + (OmegaR)*exp(-4*x) + OmegaLambda);
  double Omega_value2 = (-3 * (OmegaB+OmegaCDM)*exp(-3*x) - 4 * (OmegaR)*exp(-4*x));



  double dHdx = 0.5 * H0 * sqrt(Omega_value1) * Omega_value2;


  return dHdx;
}

double BackgroundCosmology::ddHddx_of_x(double x) const{


  double Omega_value1 = ((OmegaB+OmegaCDM)*exp(-3*x) + (OmegaR)*exp(-4*x) + OmegaLambda);
  double Omega_value2 = (-3 * (OmegaB+OmegaCDM)*exp(-3*x) - 4 * (OmegaR)*exp(-4*x));
  double Omega_value3 = (9*(Omega_value3+OmegaCDM)*exp(-3*x) + 16*OmegaR*exp(-4*x));



  double ddHddx = 0.5 * H0 * (0.5 * pow(Omega_value1,(-3/2))) * Omega_value2 + 0.5*H0*1/sqrt(Omega_value1) * Omega_value3;

  return ddHddx;
}

double BackgroundCosmology::Hp_of_x(double x) const{

  double Omega_value_Hp1 = ((OmegaB+OmegaCDM)*exp(-x) + OmegaR*exp(-2*x) + OmegaLambda*exp(2*x) + OmegaK);

  double Hp = H0 * sqrt(Omega_value_Hp1);

  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  double Omega_value_Hp1 = ((OmegaB+OmegaCDM)*exp(-x) + OmegaR*exp(-2*x) + OmegaLambda*exp(2*x) + OmegaK);
  double Omega_value_Hp2 = (- (OmegaB+OmegaCDM) *exp(-x) -2 * OmegaR*exp(-2*x) + 2*OmegaLambda*exp(2*x));



  double dHpdx = 0.5 * H0 * 1/sqrt(Omega_value_Hp1) * Omega_value_Hp2;

  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{


  double Omega_value_Hp1 = ((OmegaB+OmegaCDM)*exp(-x) + OmegaR*exp(-2*x) + OmegaLambda*exp(2*x) + OmegaK);
  double Omega_value_Hp2 = (- (OmegaB+OmegaCDM) *exp(-x) -2 * OmegaR*exp(-2*x) + 2*OmegaLambda*exp(2*x));
  double Omega_value_Hp3 = ( (OmegaB+OmegaCDM)*exp(-x) + 4*OmegaR*exp(-2*x) + 4*OmegaLambda*exp(2*x)  );



  double ddHpddx = 0.5 * H0 / sqrt(Omega_value_Hp1) * (Omega_value_Hp3-0.5*Omega_value_Hp2*Omega_value_Hp2/Omega_value_Hp1);


  return ddHpddx;
}

double BackgroundCosmology::get_OmegaB(double x) const{ //Finding OmegaB
  if(x == 0.0) return OmegaB;


  return OmegaB * pow(H0,2) / (exp(3*x) * (pow(H_of_x(x),2)));

}

double BackgroundCosmology::get_OmegaR(double x) const{ //Finding OmegaR
  if(x == 0.0) return OmegaR;


  return OmegaR / (exp(4*x) * (pow(H_of_x(x),2))/pow(H0,2));
}

double BackgroundCosmology::get_OmegaNu(double x) const{ //Finding OmegaNu
  if(x == 0.0) return OmegaNu;

  return OmegaNu / (pow(exp(x),4)*(pow(H_of_x(x),2))/pow(H0,2));


}

double BackgroundCosmology::get_OmegaCDM(double x) const{ //Finding OmegaCDM
  if(x == 0.0) return OmegaCDM;

  return OmegaCDM / (exp(3*x)*(pow(H_of_x(x),2))/pow(H0,2));

}

double BackgroundCosmology::get_OmegaLambda(double x) const{ //Finding OmegaLambda
  if(x == 0.0) return OmegaLambda;


  return OmegaLambda / ((pow(H_of_x(x),2))/pow(H0,2));

}

double BackgroundCosmology::get_OmegaK(double x) const{ //Finding OmegaK
  if(x == 0.0) return OmegaK;

  return OmegaK / (exp(2*x)*(pow(H_of_x(x),2))/pow(H0,2));


}

double BackgroundCosmology::eta_of_x(double x) const{ //Finding eta
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{ //Finding t
  return t_of_x_spline(x);
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



double BackgroundCosmology::get_z(double x) const{ //Redshift
  return 1/exp(x) - 1;
}
double BackgroundCosmology::get_chi(double x) const{ //Finding chi
  return (eta_of_x_spline(Constants.xhi0) - eta_of_x(x));
}

double BackgroundCosmology::get_dA(double x) const{ //Finding dA
  return exp(x) * get_chi(x);
}

double BackgroundCosmology::get_dL(double x) const{ //Finding dL
  return get_dA(x) / exp(2*x);
}

double BackgroundCosmology::t_of_0(double x) const{
  return t_of_x(0);
}

double BackgroundCosmology::time_today(double x) const{
  return eta_of_x(0)/Constants.c;
}

double BackgroundCosmology::eta_Hp_c(double x) const{
  return eta_of_x(x)*Hp_of_x(x)/Constants.c;
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
  const double x_min = log(1e-8);
  const double x_max =  2.0;
  const int    n_pts =  1000;

  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  double a_rm = OmegaR/(OmegaB+OmegaCDM); //Scalefactor a in Radiotion matter equality
  double x_rm = log(a_rm); //x in Radiotion matter equality
  double z_rm = 1/a_rm - 1; //Redshift z in Radiotion matter equality
  double t_rm = t_of_x((x_rm))/(pow(10,9)*365*24*60*60); //t in Radiotion matter equality

  double a_MDE_ = (pow(( (get_OmegaB(0)+get_OmegaCDM(0) ) / get_OmegaLambda(0)  ),(1./3))); //Scalefactor a in matter-dark energy equality
  double x_MDE_ = log(a_MDE_); //x in matter-dark energy equality
  double z_MDE_ = 1/a_MDE_ - 1; //Redshift z in matter-dark energy equality
  double t_MDE_ = t_of_x(log(a_MDE_))/(pow(10,9)*365*24*60*60); //t in matter-dark energy equality

  double a_acc = pow( (OmegaB+OmegaCDM)/ (2*OmegaLambda), (1./3)  ); //Scalefactor a when the universe starts to accelerate
  double x_acc = log(a_acc); //x when the universe starts to accelerate
  double z_acc = 1/a_acc - 1; //Redshift z when the universe starts to accelerate
  double t_acc = t_of_x(log(a_acc))/(pow(10,9)*365*24*60*60); //t when the universe starts to accelerate



  double age_today = t_of_x(0)/(pow(10,9)*365*24*60*60); //Universes age today
  double CT_today = eta_of_x(0)/Constants.c/(pow(10,9)*365*24*60*60); //Conformal time today


  std::cout << "omegaR " << OmegaR << "\n";
  std::cout << "omegaB " << OmegaB << "\n";
  std::cout << "omegaCDM " << OmegaCDM << "\n";
  std::cout << "a_rm: " << a_rm << "\n";
  std::cout << "x_rm: " << x_rm << "\n";
  std::cout << "z_rm: " << z_rm << "\n";
  std::cout << "t_rm: " << t_rm << "\n";

  std::cout << "x_MDE: " << x_MDE_ << "\n";
  std::cout << "z_MDE: " << z_MDE_ << "\n";
  std::cout << "t_MDE: " << t_MDE_ << "\n";

  std::cout << "x_acc: " << x_acc << "\n";
  std::cout << "z_acc: " << z_acc << "\n";
  std::cout << "t_acc: " << t_acc << "\n";

  std::cout << "Age of the Universe today: " << age_today << '\n';
  std::cout << "Conformal time today: " << CT_today << '\n';

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)/Constants.Mpc       << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp << t_of_x(x)        << " ";
    fp << ddHpddx_of_x(x) << " ";
    fp << H_of_x(x) << " ";
    fp << H0 << " ";
    fp << get_z(x) << " ";
    fp << get_dA(x)/(Constants.Mpc) << " ";
    fp << get_dL(x)/(Constants.Mpc) << " ";
    fp << get_chi(x)/(Constants.Mpc) << " ";
    fp << t_of_0(x) << " ";
    fp << time_today(x) << " ";
    fp << eta_Hp_c(x) << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
