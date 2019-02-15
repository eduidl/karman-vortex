#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "Karman.h"

namespace karman_vortex {

constexpr double Re = 70.0;  // Reynolds Number
constexpr double cfl = 0.2;  // CFL Number

/*SOR Pamameters*/
constexpr double omegap = 1.0;
constexpr int maxitp = 100;
constexpr double errorp = 0.0001;

/* set x-grid parameters*/
constexpr int i_1 = 96;
constexpr int i_2 = 106;

/* set y-grid parameters*/
constexpr int j_1 = 96;
constexpr int j_2 = 106;

/* set delta x,y,t*/
constexpr double dx = 1.0 / (i_2 - i_1);
constexpr double dy = 1.0 / (j_2 - j_1);
constexpr double dt = cfl * fmin(dx, dy);

inline bool within_square(int x, int y) {
  return i_1 < x && x < i_2 && j_1 < y && y < j_2;
}

void Grid::BoundaryConditionP() {
  // 上下
  for (int j = 1; j <= my; j++) {
    p_[1][j] = 0.0;
    p_[mx][j] = 0.0;
  }

  // 左右
  for (int i = 1; i <= mx; i++) {
    p_[i][1] = 0.0;
    p_[i][my] = 0.0;
  }

  // 箱の隅
  p_[i_1][j_1] = p_[i_1 - 1][j_1 - 1];
  p_[i_1][j_2] = p_[i_1 - 1][j_2 + 1];
  p_[i_2][j_1] = p_[i_2 + 1][j_1 - 1];
  p_[i_2][j_2] = p_[i_2 + 1][j_2 + 1];

  // 箱の上下
  for (int j = j_1 + 1; j < j_2; j++) {
    p_[i_1][j] = p_[i_1 - 1][j];
    p_[i_2][j] = p_[i_2 + 1][j];
  }

  // 箱の左右
  for (int i = i_1 + 1; i < i_2; i++) {
    p_[i][j_1] = p_[i][j_1 - 1];
    p_[i][j_2] = p_[i][j_2 + 1];
  }
}

void Grid::BoundaryConditionUV() {
  // 上下
  for (int j = 1; j <= my; j++) {
    u_[1][j] = 1.0;
    u_[0][j] = 1.0;

    v_[1][j] = 0.0;
    v_[0][j] = 0.0;

    u_[mx][j] = 2.0 * u_[mx - 1][j] - u_[mx - 2][j];
    u_[mx + 1][j] = 2.0 * u_[mx][j] - u_[mx - 1][j];

    v_[mx][j] = 2.0 * v_[mx - 1][j] - v_[mx - 2][j];
    v_[mx + 1][j] = 2.0 * v_[mx][j] - v_[mx - 1][j];
  }

  // 左右
  for (int i = 1; i <= mx; i++) {
    u_[i][1] = 2.0 * u_[i][2] - u_[i][3];
    u_[i][0] = 2.0 * u_[i][1] - u_[i][2];

    v_[i][1] = 2.0 * v_[i][2] - v_[i][3];
    v_[i][0] = 2.0 * v_[i][1] - v_[i][2];

    u_[i][my] = 2.0 * u_[i][my - 1] - u_[i][my - 2];
    u_[i][my + 1] = 2.0 * u_[i][my] - u_[i][my - 1];

    v_[i][my] = 2.0 * v_[i][my - 1] - v_[i][my - 2];
    v_[i][my + 1] = 2.0 * v_[i][my] - v_[i][my - 1];
  }

  // 箱の中
  for (int i = i_1; i <= i_2; i++) {
    for (int j = j_1; j <= j_2; j++) {
      u_[i][j] = 0.0;
      v_[i][j] = 0.0;
    }
  }
}

void Grid::PoissonEquation() {
  vector_2d rhs(mx + 2, std::vector<double>(my + 2, 0.0));
  for (int i = 2; i < mx; i++) {
    for (int j = 2; j < my; j++) {
      if (within_square(i, j)) continue;

      double ux = (u_[i + 1][j] - u_[i - 1][j]) / (2.0 * dx);
      double uy = (u_[i][j + 1] - u_[i][j - 1]) / (2.0 * dy);
      double vx = (v_[i + 1][j] - v_[i - 1][j]) / (2.0 * dx);
      double vy = (v_[i][j + 1] - v_[i][j - 1]) / (2.0 * dy);
      rhs[i][j] = (ux + vy) / dt - (ux * ux + 2.0 * uy * vx + vy * vy);
    }
  }

  for (int iter = 1; iter <= maxitp; iter++) {
    double res = 0.0;
    for (int i = 2; i < mx; i++) {
      for (int j = 2; j < my; j++) {
        if (within_square(i, j)) continue;

        double dp = ((p_[i + 1][j] + p_[i - 1][j]) / (dx * dx) +
                     (p_[i][j + 1] + p_[i][j - 1]) / (dy * dy) - rhs[i][j]) /
                        (2.0 / (dx * dx) + 2.0 / (dy * dy)) -
                    p_[i][j];
        res += dp * dp;
        p_[i][j] += omegap * dp;
      }
    }
    BoundaryConditionP();
    res = sqrt(res / double(mx * my));
    if (res < errorp) break;
  }
}

// Kawamura-Kuwahara スキーム
inline double KKSchemeX(const double u, const vector_2d &f, const int i,
                        const int j) {
  return (u * (-f[i + 2][j] + 8.0 * (f[i + 1][j] - f[i - 1][j]) + f[i - 2][j]) /
              (12.0 * dx) +
          abs(u) * (f[i + 2][j] - 4.0 * (f[i + 1][j] + f[i - 1][j]) +
                    6.0 * f[i][j] + f[i - 2][j]) /
              (4.0 * dx));
}

// Kawamura-Kuwahara スキーム
inline double KKSchemeY(const double v, const vector_2d &f, const int i,
                        const int j) {
  return (v * (-f[i][j + 2] + 8.0 * (f[i][j + 1] - f[i][j - 1]) + f[i][j - 2]) /
              (12.0 * dy) +
          abs(v) * (f[i][j + 2] - 4.0 * (f[i][j + 1] + f[i][j - 1]) +
                    6.0 * f[i][j] + f[i][j - 2]) /
              (4.0 * dy));
}

void Grid::VelocityEquation() {
  vector_2d urhs(mx + 2, std::vector<double>(my + 2, 0.0));
  vector_2d vrhs(mx + 2, std::vector<double>(my + 2, 0.0));

  for (int i = 2; i < mx; i++) {
    for (int j = 2; j < my; j++) {
      if (within_square(i, j)) continue;

      urhs[i][j] =
          -(p_[i + 1][j] - p_[i - 1][j]) / (2.0 * dx) +
          (u_[i + 1][j] - 2.0 * u_[i][j] + u_[i - 1][j]) / (Re * dx * dx) +
          (u_[i][j + 1] - 2.0 * u_[i][j] + u_[i][j - 1]) / (Re * dy * dy) -
          KKSchemeX(u_[i][j], u_, i, j) - KKSchemeY(v_[i][j], u_, i, j);
      vrhs[i][j] =
          -(p_[i][j + 1] - p_[i][j - 1]) / (2.0 * dy) +
          (v_[i + 1][j] - 2.0 * v_[i][j] + v_[i - 1][j]) / (Re * dx * dx) +
          (v_[i][j + 1] - 2.0 * v_[i][j] + v_[i][j - 1]) / (Re * dy * dy) -
          KKSchemeX(u_[i][j], v_, i, j) - KKSchemeY(v_[i][j], v_, i, j);
    }
  }

  for (int i = 2; i < mx; i++) {
    for (int j = 2; j < my; j++) {
      if (within_square(i, j)) continue;

      u_[i][j] += dt * urhs[i][j];
      v_[i][j] += dt * vrhs[i][j];
    }
  }
}

void Grid::WriteP(const std::string &file_name) {
  std::ofstream ofs(file_name);
  if (!ofs) {
    std::cerr << "ファイルを開けません" << std::endl;
    exit(-1);
  }

  for (int i = my; i >= 1; i--) {
    for (int j = 1; j <= mx; j++) {
      ofs << p_[j][i];
      ofs << ((j == mx) ? "\n" : ",");
    }
  }
}
}
