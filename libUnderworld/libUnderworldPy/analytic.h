
void solcx(
  double pos[], 
  double _eta_A,
  double _eta_B, /* Input parameters: density, viscosity A, viscosity B */
  double _x_c,
  int _n, /* Input parameters: viscosity jump location, wavenumber in x */
  double vel[],
  double* presssure, 
  double total_stress[],
  double strain_rate[]   );

void solkx( 
  double  pos[],
  double  _sigma, /* density */
  double  _m,
  int     _n, /* wavelength in z, wavenumber in x */
  double  _B, /* viscosity parameter */
  double  vel[],
  double* presssure, 
  double  total_stress[],
  double  strain_rate[],
  double* viscosity );

void solkz( 
  double  pos[],
  double  _sigma, /* density */
  double  _km,
  int     _n, /* wavelength in z, wavenumber in x */
  double  _B, /* viscosity parameter */
  double  vel[],
  double* presssure, 
  double  total_stress[],
  double  strain_rate[] );
