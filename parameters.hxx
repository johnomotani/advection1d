struct Parameters {
  /// Number of spatial grid points
  int Nz = 32;

  /// Size of spatial domain
  double L = 1.0;

  /// internal time-step
  double dt = .01;

  /// time between outputs
  double t_out = 1.0;

  /// number of outputs
  int N_out = 100;
};

const Parameters createParameters();
