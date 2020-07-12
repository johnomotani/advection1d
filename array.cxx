#include "array.hxx"

Array createArray(const int N) {
  Array result;
  result.resize(N);

  // initialise to 0.0
  std::fill(result.begin(), result.end(), 0.0);

  // ensure result has the minimum amount of memory allocated
  result.shrink_to_fit();

  return result;
}
