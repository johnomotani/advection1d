// Convenience stuff for working with std::vector

#ifndef __ARRAY_H__
#define __ARRAY_H__

#include <fftw3.h>
#include <limits>
#include <vector>

// Copied from template at https://en.cppreference.com/w/cpp/named_req/Allocator
// Modified to use fftw_malloc() and fftw_free()
// Allows us to create a std::vector which has memory allocated by fftw_malloc,
// required for correct alignment by FFTW, which (sometimes) allows better
// performance.
template <class T> struct FFTWAllocator {
  typedef T value_type;
  FFTWAllocator() = default;
  template <class U>
  constexpr FFTWAllocator(const FFTWAllocator<U> &) noexcept {}

  [[nodiscard]] T *allocate(std::size_t n) {
    if (n > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
      throw std::bad_alloc();
    }

    if (auto p = static_cast<T *>(fftw_malloc(n * sizeof(T)))) {
      return p;
    }

    throw std::bad_alloc();
  }

  void deallocate(T *p, std::size_t) noexcept { fftw_free(p); }
};

template <class T, class U>
bool operator==(const FFTWAllocator<T> &, const FFTWAllocator<U> &) {
  return true;
}
template <class T, class U>
bool operator!=(const FFTWAllocator<T> &, const FFTWAllocator<U> &) {
  return false;
}

using Array = std::vector<double, FFTWAllocator<double>>;

Array createArray(const int N);

#endif // __ARRAY_H__
