#include <iostream>
#include "la_array.h"
#include "la_ops.h"

// int main() {
//   la::Array<double, 2> A(std::vector<int>(2, 2), 2);
//   std::vector<int> dims(2);
//   dims[0] = 1;
//   dims[1] = 2;
//   la::Array<double, 1> B(2, 1);
//   la::Array<double, 1> C = A*B;
//   for(int i = 0 ; i < C.get_dims()[0]; i++) {
//     // for(int j = 0; j < C.get_dims()[1]; j++) {
//       std::cout << C[i] << ", ";
//     // }
//     std::cout << std::endl;
//   }
// }
