#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================

RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo,
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{

}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){

  // Compute and spline Xe, ne
  solve_number_density_electrons();

  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");

  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);
  Vector Xe_saha_arr(npts_rec_arrays);

  // Physical constants and cosmological parameters
  const double G           = Constants.G;
  const double OmegaB      = cosmo->get_OmegaB();
  const double H0          = cosmo->get_H0();
  const double m_H         = Constants.m_H;
  double rho_c0            = 3*pow(H0,2)/(8*M_PI*G);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;
    const double Xe_saha_current = Xe_ne_data.first;

    double Xe_current_peebles;
    double ne_current_peebles;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){

      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      //...
      //...
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;
      Xe_saha_arr[i] = Xe_saha_current;

    } else {




      //==============================================================
      // TODO: Compute X_e from current time til today by solving
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };

      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result
      //=============================================================================
      //...
      //...
      double x_tmp_1 = x_array[i-1];
      double x_tmp_2 = x_array[i];
      double y_ic = Xe_arr[i-1];

      Vector x_tmp_arr{x_tmp_1,x_tmp_2};
      Vector Xe_ic{y_ic};
      peebles_Xe_ode.solve(dXedx, x_tmp_arr, Xe_ic);

      double n_H = OmegaB*rho_c0/(m_H*exp(3*x_array[i]));

      Xe_current_peebles = peebles_Xe_ode.get_data_by_component(0)[1];
      ne_current_peebles =  n_H*Xe_current_peebles;

      Xe_arr[i] = Xe_current_peebles;
      ne_arr[i] = ne_current_peebles;

    }


  }




  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x
  // functions are working
  //=============================================================================
  //...
  //...
  auto log_ne_arr = log(ne_arr);

  Xe_of_x_spline.create(x_array, Xe_arr, "Xe");
  ne_of_x_spline.create(x_array, log_ne_arr, "ne");
  Xe_saha_spline.create(x_array,Xe_saha_arr, "Xe_saha");


  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);

  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;

  // Fetch cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB();
  const double OmegaR      = cosmo->get_OmegaR();
  const double OmegaNu     = cosmo->get_OmegaNu();
  const double OmegaCDM    = cosmo->get_OmegaCDM();
  const double OmegaLambda = cosmo->get_OmegaLambda();
  const double OmegaK      = cosmo->get_OmegaK();
  const double H0          = cosmo->get_H0();
  const double h           = cosmo->get_h();
  const double Neff        = cosmo->get_Neff();
  const double TCMB        = cosmo->get_TCMB();


  const double rho_c0      = 3*pow(H0,2)/(8*M_PI*G);
  const double n_H         = OmegaB*rho_c0/(m_H*pow(a,3));
  const double T_b         = TCMB/a;

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;

  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  //...
  //...
  //Solution from Saha eq.
  double c = 1.0/ n_H * pow(m_e*k_b*T_b / (2.0*M_PI*pow(hbar,2.0)),3.0/2.0) * exp(-epsilon_0/(k_b*T_b));

  if(c > 1e14){
    Xe = 2.0/(sqrt(1.0+4.0/c)+1);
  }
  else{
    Xe = -c/2.0 + sqrt(pow(c, 2.0) + 4.0*c)/2.0;
  }

  ne = Xe*n_H;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;
  const double alpha       = 1/137.0359992;

  // Cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB();
  const double OmegaR      = cosmo->get_OmegaR();
  const double OmegaNu     = cosmo->get_OmegaNu();
  const double OmegaCDM    = cosmo->get_OmegaCDM();
  const double OmegaLambda = cosmo->get_OmegaLambda();
  const double OmegaK      = cosmo->get_OmegaK();
  const double H0          = cosmo->get_H0();
  const double h           = cosmo->get_h();
  const double Neff        = cosmo->get_Neff();
  const double TCMB        = cosmo->get_TCMB();
  const double H           = cosmo->H_of_x(x);

  const double rho_c0      = 3*pow(H0,2)/(8*M_PI*G);
  const double T_b         = TCMB/a;

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  //...
  //...
  double n_H             = (1-Yp) * OmegaB*rho_c0/(m_H*pow(a,3));
  double n_1s            = (1 - X_e)*n_H;
  double lambda_alpha    = H*(pow(3*epsilon_0,3.0)) / (pow(8*M_PI,2.0)*n_1s) / pow(c*hbar,3.0);
  double phi_2           = 0.448 * log(epsilon_0/(k_b*T_b));
  double alpha2          = 64 * M_PI / sqrt(27*M_PI) * pow(alpha/m_e,2) * sqrt(epsilon_0/(k_b*T_b)) * phi_2 * pow(hbar,2)/c;
  double beta            = alpha2 * pow(m_e*k_b*T_b/(2*M_PI),3.0/2.0)*exp(-epsilon_0/(k_b*T_b)) * pow(1/pow(hbar,2.0),3.0/2.0);
  double beta2           = alpha2 * pow(m_e*k_b*T_b/(2*M_PI),3.0/2.0)*exp(-epsilon_0/(4*k_b*T_b)) * pow(1/pow(hbar,2.0),3.0/2.0);
  double Cr              = (lambda_2s1s + lambda_alpha)/(lambda_2s1s + lambda_alpha + beta2);

  dXedx[0]               = Cr / H * (beta * (1-X_e) - n_H * alpha2 * pow(X_e,2));

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...
    double H       = cosmo->H_of_x(x);

    // Set the derivative for photon optical depth
    dtaudx[0] = -Constants.c*exp(log_ne_of_x(x))*Constants.sigma_T/H;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...
  Vector x_tmp_arr = Utils::linspace(x_end, x_start, npts_rec_arrays);
  Vector y_ic{0.0};                                                             // Sjekk at initialbetingelsene her er riktig!!!
  ODESolver tau_ode;
  tau_ode.solve(dtaudx, x_tmp_arr, y_ic);
  auto tau_array_ode = tau_ode.get_data_by_component(0);


  Vector x_array(npts_rec_arrays);
  Vector tau_array(npts_rec_arrays);
  Vector dtaudx_array(npts_rec_arrays);

  double H;

  int j = 0;
  while(j<npts_rec_arrays){
    x_array[j] = x_tmp_arr[npts_rec_arrays-j-1];
    tau_array[j] = tau_array_ode[npts_rec_arrays-j-1];
    H = cosmo->H_of_x(x_array[j]);
    dtaudx_array[j] = -Constants.c*exp(log_ne_of_x(x_array[j]))*Constants.sigma_T/H;
    j++;
  }


  tau_of_x_spline.create(x_array, tau_array, "tau");
  dtaudx_spline.create(x_array, dtaudx_array, "dtaudx");

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...
  Vector g_tilde_arr(npts_rec_arrays);

  int k = 0;
  while(k<npts_rec_arrays){
    H = cosmo->H_of_x(x_array[k]);
    g_tilde_arr[k] = Constants.c*exp(log_ne_of_x(x_array[k]))*Constants.sigma_T/H*exp(-tau_of_x_spline(x_array[k]));
    k++;
  }

  g_tilde_of_x_spline.create(x_array, g_tilde_arr, "g_tilde_of_x");

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...

  return dtaudx_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return dtaudx_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return Xe_of_x_spline(x);
}

double RecombinationHistory::log_ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return ne_of_x_spline(x);
}

double RecombinationHistory::get_x_decoupling() const{
  // Compute time of decoupling
  return Utils::binary_search_for_value(tau_of_x_spline, 1);
}

double RecombinationHistory::get_a_decoupling() const{
  return exp(get_x_decoupling());
}

double RecombinationHistory::get_z_decoupling() const{
  return 1/get_a_decoupling() -1;
}

double RecombinationHistory::get_x_recombination() const{
  return Utils::binary_search_for_value(Xe_of_x_spline, 0.5);
}

double RecombinationHistory::get_a_recombination() const{
  return exp(get_x_recombination());
}

double RecombinationHistory::get_z_recombination() const{
  return 1/get_a_recombination() -1;
}


double RecombinationHistory::x_saha_predicton() const{
  return Utils::binary_search_for_value(Xe_saha_spline, 0.5);
}

double RecombinationHistory::a_saha_predicton() const{
  return exp(x_saha_predicton());
}

double RecombinationHistory::z_saha_predicton() const{
  return 1/a_saha_predicton() - 1;
}

double RecombinationHistory::Freeze_out_abundance() const{
  return Xe_of_x(0.0);
}


//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";

  std::cout << "x_decoupling " << get_x_decoupling()<<"\n";
  std::cout << "z_decoupling " << get_z_decoupling()<<"\n";
  std::cout << "x_recoupling " << get_x_recombination()<<"\n";
  std::cout << "z_recoupling " << get_z_recombination()<<"\n";
  std::cout << "x Saha prediction " << x_saha_predicton()<<"\n";
  std::cout << "z Saha prediction " << z_saha_predicton()<<"\n";
  std::cout << "Freeze-out abundance " << Freeze_out_abundance()<<"\n";

  std::cout << std::endl;
}

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 2000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << log_ne_of_x(x)       << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << get_x_decoupling() << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
