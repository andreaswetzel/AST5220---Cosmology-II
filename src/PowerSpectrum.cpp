#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo,
    RecombinationHistory *rec,
    Perturbations *pert):
  cosmo(cosmo),
  rec(rec),
  pert(pert)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  Vector k_array     = Utils::linspace(k_min, k_max, n_k);
  Vector log_k_array = log(k_array);

  //=========================================================================
  // TODO: Make splines for j_ell.
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array);

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");

  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");

  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());

  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================
  double xmax  = k_max*cosmo->eta_of_x(0);
  int N        = xmax/(2.0*M_PI/16.0);
  Vector x     = Utils::linspace(0, xmax, N);

  Vector j_ell(N);

  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];

    int j = 0;
    while(j < N){
      j_ell[j] = Utils::j_ell(ell, x[j]);
      j++;
    }

    // Make the j_ell_splines[i] spline
    j_ell_splines[i].create(x, j_ell);
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array,
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));
  double F_ell;


  double xmin   = Constants.x_start;
  double xmax   = 0.0;
  double deltax = 0.03;
  int N         = (xmax - xmin)/deltax;
  Vector x      = Utils::linspace(xmin, xmax, N);
  double eta0   = cosmo->eta_of_x(0);


  for(size_t ik = 0; ik < k_array.size(); ik++){

    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral
    // F_ell(k) = Int deltax jell(k(eta-eta0)) * S(x,k) for all the ell values for the
    // given value of k
    //=============================================================================
    double k = k_array[ik];
    for(size_t iell = 0; iell < ells.size(); iell++){
      // Using the trapezoidal rule:
      F_ell = 0.5*pert->get_Source_T(x[0], k)*j_ell_splines[iell](k*(eta0 - cosmo->eta_of_x(x[0])));

      int ix = 0;
      while(ix < (N-1)){
        F_ell = F_ell + pert->get_Source_T(x[ix], k)*j_ell_splines[iell](k*(eta0 - cosmo->eta_of_x(x[ix])));
        ix++;
      }
      result[iell][ik] = (F_ell + 0.5*pert->get_Source_T(x[N-1], k)*j_ell_splines[iell](k*(eta0 - cosmo->eta_of_x(x[N-1]))))*deltax;
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();

  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  for(int i = 0; i < nells; i++){
    thetaT_ell_of_k_spline[i].create(k_array, thetaT_ell_of_k[i]);
  }

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
}

//====================================================
// Compute Cell (could be TT or TE or EE)
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell deltak/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell deltak/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================

  double deltak = 2.0*M_PI/(32.0*cosmo->eta_of_x(0.0));
  int N     = (k_max-k_min)/deltak;                         // Number of k values
  double eta0 = cosmo->eta_of_x(0.0);
  Vector k  = Utils::linspace(k_min, k_max, N);
  Vector cell(nells);
  Vector result(nells);

  // Integration with trapezoidal rule:
  for(int i = 0; i < nells; i++){
    cell[i] = 0.5*primordial_power_spectrum(k[0])*f_ell_spline[i](k[0])*g_ell_spline[i](k[0])/k[0];
    int j = 0;
    while(j < N-1){
      cell[i] = cell[i] + primordial_power_spectrum(k[j])*f_ell_spline[i](k[j])*g_ell_spline[i](k[j])/k[j];
      j++;
    }
    result[i] = 4.0*M_PI*deltak*(cell[i] + primordial_power_spectrum(k[N-1])*f_ell_spline[i](k[N-1])*g_ell_spline[i](k[N-1])/k[N-1]);
  }

  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k) const{
  double pofk = 0.0;

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================
  double OmegaM  = cosmo->get_OmegaCDM() + cosmo->get_OmegaB();
  double Phi     = pert->get_Phi(0, k);
  double H0      = cosmo->get_H0();

  double  DeltaM = pow(Constants.c,2.0)*pow(k,2.0) * Phi / (3.0/2.0 * OmegaM*pow(exp(x),-1.0)*pow(H0,2.0));
  double Delta   = pow(Constants.c*k,2)*Phi*exp(x)/(1.5*OmegaM*pow(H0,2.0));


  // pofk           = 2.0*pow(M_PI,2.0)/pow(k, 3.0)*primordial_power_spectrum(k)*pow(abs(Delta),2.0);


  pofk = pow(abs(DeltaM),2.0) * primordial_power_spectrum(k);

  return pofk;
}

double PowerSpectrum::get_k_eq() const{
  double a_eq = cosmo->get_OmegaR()/(cosmo->get_OmegaB()+cosmo->get_OmegaCDM());
  double x_eq = log(a_eq);
  double k_eq = cosmo->Hp_of_x(x_eq) / Constants.c * Constants.Mpc / cosmo->get_h();
  return k_eq;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================


void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI)
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}



void PowerSpectrum::output_theta_ell(std::string filename) const{
  const int npts = 20000.0;
  std::ofstream fp(filename.c_str());
  auto k_array = exp(Utils::linspace(log(k_min), log(k_max), npts));

  fp << "# k (1/Mpc)  Theta_ell's\n";
  auto print_data = [&](const double k) {

    // 1
    fp << k * Constants.Mpc << " ";
    fp << k * Constants.c/cosmo->get_H0() << " ";

    // 2: Theta_ell
    fp << thetaT_ell_of_k_spline[4](k) << " ";
    fp << thetaT_ell_of_k_spline[19](k) << " ";
    fp << thetaT_ell_of_k_spline[24](k) << " ";
    fp << thetaT_ell_of_k_spline[32](k) << " ";
    fp << thetaT_ell_of_k_spline[42](k) << " ";

    fp << cosmo->eta_of_x(0.0) << " ";
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}


void PowerSpectrum::output_matter(std::string filename) const{
  // Print matter power spectrum:
  std::ofstream fp(filename.c_str());
  auto k_array = Utils::linspace(k_min, k_max, n_k);  // [Mpc]
  auto print_P_data = [&](const double k) {
    fp << k*Constants.Mpc/cosmo->get_h()                                  << " ";
    fp << get_matter_power_spectrum(0.0, k) * 2 * pow(M_PI,2) /pow(Constants.Mpc/cosmo->get_h() * k,3.0)  << " ";
    fp << get_k_eq() << " ";
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_P_data);
}
