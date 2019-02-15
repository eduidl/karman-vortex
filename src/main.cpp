#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "Karman.h"

int main() {
  namespace chrono = std::chrono;

  const auto nlast = 10000;

  karman_vortex::Grid grid;
  grid.BoundaryConditionP();
  grid.BoundaryConditionUV();

  auto start = chrono::system_clock::now();
  for (int i = 1; i <= nlast; ++i) {
    grid.PoissonEquation();
    grid.VelocityEquation();
    grid.BoundaryConditionUV();

    const auto end = chrono::system_clock::now();
    const auto elapsed_time =
        chrono::duration_cast<chrono::milliseconds>(end - start).count() /
        1000.0;
    std::cout << "step:" << i << "/" << nlast
              << " elapsed_time:" << elapsed_time << "[s]" << std::endl;

    std::ostringstream ss;
    ss << "../data/" << std::setw(6) << std::setfill('0') << i << ".csv";
    grid.WriteP(ss.str());
  }
}
