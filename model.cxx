#include "parameters.hxx"

Model::Model(Parameters parameters)
    : Nz(parameters.Nz), L(parameters.L), dz(parameters.L/parameters.Nz) {}
