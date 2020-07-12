class Parameters; // Forward-declare to avoid including header here

class Model {
  Model(Parameters);

private:
  const int Nz;
  const double L;
  const double dz;
};
