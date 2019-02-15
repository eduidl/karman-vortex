#pragma once

#include <cmath>
#include <string>
#include <vector>

namespace karman_vortex {
const int mx = 401;
const int my = 201;

using vector_2d = std::vector<std::vector<double>>;

class Grid {
 public:
  Grid()
      : u_(mx + 2, std::vector<double>(my + 2, 1.0)),
        v_(mx + 2, std::vector<double>(my + 2, 0.0)),
        p_(mx + 2, std::vector<double>(my + 2, 0.0)){};
  ~Grid() = default;

  void BoundaryConditionP();
  void BoundaryConditionUV();
  void PoissonEquation();
  void VelocityEquation();
  void WriteP(const std::string &file_name);

 private:
  vector_2d u_;
  vector_2d v_;
  vector_2d p_;
};
}
