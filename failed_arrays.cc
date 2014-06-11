.cc:

  // // doesn't seem good
  // template <typename T> Array<T>::Array(): dims_(std::vector<int>()), n_(-1) {
  // }

  // // template <typename T> Array<T>::Array(const std::vector<int>* dims): dims_(*dims), n_(dims->size()) {
  // template <typename T> Array<T>::Array(const std::vector<int> dims): dims_(dims), n_(dims.size())  {
  //   if(n_ > 2) {
  //     A = new Array<T>*[dims.back()];
  //     std::vector<int> next_dims = dims;
  //     next_dims.pop_back();
  //     for(int i = 0; i < dims.back(); i++) {
  // 	A[i] = new Array<T>(next_dims);
  //     }
  //   }
  //   else {
  //     A = new Array<T>*[dims.back()];
  //     for(int i = 0; i < dims.back(); i++) {
  // 	A[i] = new Vector<T>(dims[0]);
  //     }
  //   }
  // }

  // template <typename T> Array<T>::Array(const std::vector<int> dims, const T init_val): dims_(dims), n_(dims.size()) {
  //   if(n_ > 2) {
  //     A = new Array<T>*[dims.back()];
  //     std::vector<int> next_dims = dims;
  //     next_dims.pop_back();
  //     for(int i = 0; i < dims.back(); i++) {
  // 	A[i] = new Array<T>(next_dims, init_val);
  //     }
  //   }
  //   else {
  //     A = new Array<T>*[dims.back()];
  //     for(int i = 0; i < dims.back(); i++) {
  // 	A[i] = new Vector<T>(dims[0], init_val);
  //     }
  //   }
  // }

  // template <typename T> Vector<T>::Vector(const int size): Array<T>(), size_(size) {
  //   V = new T[size];
  // }

  // template <typename T> Vector<T>::Vector(const int size, const T init_val): Array<T>(), size_(size) {
  //   V = new T[size];
  //   for(int i = 0; i < size; i++) {
  //     V[i] = init_val;
  //   }
  // }

  // template <typename T> Array<T>& Array<T>::operator[](const unsigned int i) {

  //   return *(A[i]);

  // }

  // template <typename T> Array<T>& Vector<T>::operator[](const unsigned int i) {
  //   return 
  // }

  // template <typename T> T& Vector<T>::operator[](const unsigned int i) {
  //   return V[i];
  // }

  // template <typename T> T Array<T>::get_val(const std::vector<int> v) {
  //   std::vector<int> v2 = v;
  //   v2.erase(v2.begin());
  //   return A[v[0]]->get_val(v2);
  // }

  // template <typename T> T Vector<T>::get_val(const std::vector<int> v) {
  //   return V[v[0]];
  // }

  // // template <typename T, typename D> array<T, D>::array(const std::vector<int> dims) {
  // //   A = new array<T, D>

.h:

  // the seedy underbelly of Array

  /* template <typename T, unsigned int N> */
  /*   class array { */
  /* public: */
  /*   explicit array(const std::vector<int> dims); */
  /*   array(const std::vector<int> dims, const int init_val); */
  /*   ~array() {} */
  /*   std::vector<int> get_dims() { */
  /*     return dims_; */
  /*   } */
  /* private: */
  /*   const std::vector<int> dims_; */
  /*   const int n_; */
  /*   array** A; */
  /* }; */

  /* template <typename T> */
  /*   class array<T, 1> { */
  /* public: */
  /*   explicit array(const int dim); */
  /*   array(const int dim, const int init_val); */
  /*   ~array() {} */
  /*   int get_dim() { */
  /*     return size; */
  /*   } */
  /* private: */
  /*   const int size; */
  /*   T* A; */
  /* }; */




  /* template <typename T>   */
  /*   class Array { */
  /* public: */
  /*   Array(); */
  /*   explicit Array(const std::vector<int> dims); */
  /*   /\* explicit Array(const std::vector<int>* dims); *\/ */
  /*   Array(const std::vector<int> dims, const T init_val); */
  /*   virtual Array<T>& operator[](const unsigned int i); */
  /*   ~Array() {} */

  /*   virtual T get_val(const std::vector<int> v); */
  /*   std::vector<int> get_dims() { */
  /*     return dims_; */
  /*   } */
  /* private: */
  /*   const std::vector<int> dims_; */
  /*   const int n_; */
  /*   Array<T>** A; */
  /* }; */

  /* template <typename T> */
  /*   class Vector : public Array<T> { */
  /* public: */
  /*   explicit Vector(const int size); */
  /*   Array<T>& operator[](const unsigned int i); */
  /*   T& operator[](const unsigned int i); */
  /*   Vector(const int size, const T init_val); */
  /*   ~Vector() {} */

  /*   T get_val(const std::vector<int> v); */
  /* private: */
  /*   const int size_; */
  /*   T* V; */
  /* };  */



