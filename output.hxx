#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <fstream>
#include <iomanip>

#include "array.hxx"

class Output {
public:
  Output()
      : f_file("advection_f.dat", std::fstream::out),
        t_file("advection_t.dat", std::fstream::out) {
    // increase precision of outputs
    f_file << std::setprecision(16);
    t_file << std::setprecision(16);
  }
  ~Output() {
    f_file.close();
    t_file.close();
  }

  void writeStep(const double t, const Array &f) {
    t_file << t << std::endl;

    for (const auto &value : f) {
      f_file << value << " ";
    }
    f_file << std::endl;
  }

private:
  std::fstream f_file;
  std::fstream t_file;
};

#endif // __OUTPUT_H__
