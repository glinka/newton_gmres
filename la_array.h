#ifndef LA_ARRAY_H
#define LA_ARRAY_H
#include <vector>

namespace la {

  template <typename T, unsigned int N>
    class Array {
    friend class Array<T, N+1>;
  public:
    Array();
    Array(const std::vector<int> dims);
    Array(const std::vector<int> dims, const T init_val);
    ~Array() {
      for(int i = 0; i < dims_[0]; i++) {
	delete A[i];
      }
      delete[] A;
    }
    
    std::vector<int> get_dims() const {
      return dims_;
    }
    
    Array<T, N-1>& operator[](unsigned int i);
    const Array<T, N-1>& operator[](unsigned int i) const;
  private:
    Array(std::vector<int>::const_iterator dims);
    Array(const std::vector<int>::const_iterator dims, const T init_val);

    const std::vector<int> dims_;
    Array<T, N-1>** A;
  };

  template <typename T>
    class Array<T, 1> {
    // why doesn't the friend function work?
    friend class Array<T, 2>;
  public:
    Array();
    Array(const unsigned int size);
    Array(const unsigned int size, const T init_val);
    Array(const std::vector<int> dims);
    Array(const std::vector<int> dims, const T init_val);
    ~Array() {
      delete[] V;
    }

    std::vector<int> get_dims() const {
      return std::vector<int>(1, size_);
    }
    int get_size() const {
      return size_;
    }

    T& operator[](unsigned int i);
    const T& operator[](unsigned int i) const;

    /* friend Array<T, 2>::Array(std::vector<int>::const_iterator dims); */
    /* friend Array<T, 2>::Array(std::vector<int>::const_iterator dims, const T init_val); */
  private:
    // make following constructor public for immediate use
    Array(std::vector<int>::const_iterator dims);
    Array(std::vector<int>::const_iterator dims, const T init_val);

    const int size_;
    T* V;
  };
}

#include "la_array.tpp"

#endif
