#include <iostream>

#include "parameters.hxx"

int main() {

  const Parameters parameters;

  std::cout<<"running 1d advection..."<<std::endl;

  std::cout<<std::endl<<"Parameters"<<std::endl<<"----------"<<std::endl;
  std::cout<<"Nz\t"<<parameters.Nz<<std::endl;
  std::cout<<"L\t"<<parameters.L<<std::endl;
  std::cout<<"dt\t"<<parameters.dt<<std::endl;
  std::cout<<"t_out\t"<<parameters.t_out<<std::endl;
  std::cout<<"N_out\t"<<parameters.N_out<<std::endl;

  return 0;
}
