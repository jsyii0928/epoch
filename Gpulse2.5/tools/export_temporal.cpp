/**
 * @file export_temporal.cpp
 * @brief 激光脉冲时域功率分布导出工具
 *
 * 计算通过焦平面 (z=0) 的总瞬时功率随时间的变化。
 * 扫描范围: -50fs 到 +50fs
 *
 * 物理公式: P(t) = Integral [ I(x,y,t) dA ]
 * 其中 I(x,y,t) = c * epsilon0 * |E|^2 (近似真空光强)
 */

#include "Gpulse.h"
#include "constants.h"
#include <H5Cpp.h>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

using namespace std;
using namespace H5;

// 网格参数 (与 export_field 保持一致以确保覆盖光斑)
const int GRID_N = 100;
const double SCAN_RANGE = RANGE; // +/- 4um

// 时间扫描参数
const double T_START = -50e-15; // -50 fs
const double T_END = 50e-15;    // +50 fs
const int T_STEPS = 200;        // 时间步数

int main() {
  cout << "启动时域功率分析工具..." << endl;
  cout << "Polarization Mode: " << (RADIAL ? "Radial (Mosaic)" : "Linear")
       << endl;

  // 1. 初始化
  Gpulse_Init();

  // 准备空间网格
  vector<double> axis_x(GRID_N), axis_y(GRID_N);
  double dx = (2 * SCAN_RANGE) / (GRID_N - 1);
  double dy = (2 * SCAN_RANGE) / (GRID_N - 1);
  double dA = dx * dy; // 面积微元

  for (int i = 0; i < GRID_N; ++i) {
    axis_x[i] = -SCAN_RANGE + i * dx;
    axis_y[i] = -SCAN_RANGE + i * dy;
  }

  // 准备时间步
  vector<double> time_points(T_STEPS);
  vector<double> power_values(T_STEPS);
  double dt = (T_END - T_START) / (T_STEPS - 1);

  for (int i = 0; i < T_STEPS; ++i) {
    time_points[i] = T_START + i * dt;
  }

  cout << "扫描配置:" << endl;
  cout << "  时间: [" << T_START * 1e15 << ", " << T_END * 1e15 << "] fs, "
       << T_STEPS << " 步" << endl;
  cout << "  空间: " << GRID_N << "x" << GRID_N << " 网格" << endl;

  // 2. 主循环 (遍历时间)
  for (int t_idx = 0; t_idx < T_STEPS; ++t_idx) {
    double t = time_points[t_idx];

    Complex Ex, Ey, Ez;
    Complex Bx, By, Bz;

    // 计算 z=0 中心点场 (x=0, y=0)
    GpulseComplex(t, 0.0, 0.0, 0.0, Ex, Ey, Ez, Bx, By, Bz);

    // 计算光强 I = 0.5 * c * eps0 * |E|^2 (Time-averaged Envelope Intensity)
    // 注意: Gpulse 计算的是 Complex Envelope Amplitude.
    // 物理光强 I = 0.5 * c * eps0 * (|Ex|^2 + |Ey|^2 + |Ez|^2)
    // Factor 0.5 comes from time-averaging of sinusoidal field OR convention
    // for peak field. However, standard I = c * eps0 * E_rms^2 = 0.5 * c * eps0
    // * E_peak^2 Gpulse returns peak envelope E_0. 所以: I_envelope = 0.5 * C *
    // EPS0 * (norm(Ex) + norm(Ey) + norm(Ez))

    double E_sq_norm =
        std::norm(Ex) + std::norm(Ey) + std::norm(Ez); // std::norm is |z|^2
    double Intensity = 0.5 * C * EPS0 * E_sq_norm;

    power_values[t_idx] = Intensity;

    // 进度条
    if (t_idx % 10 == 0 || t_idx == T_STEPS - 1) {
      printf("\r进度: %d/%d (t=%.1f fs)", t_idx + 1, T_STEPS, t * 1e15);
      fflush(stdout);
    }
  }
  cout << endl << "计算完成。" << endl;

  // 3. 保存数据到 HDF5
  string filename = "temporal_intensity.h5";
  cout << "保存结果到 " << filename << "..." << endl;

  try {
    H5File file(filename, H5F_ACC_TRUNC);

    // 保存时间数组
    hsize_t dim_t[1] = {static_cast<hsize_t>(T_STEPS)};
    DataSpace space_t(1, dim_t);
    DataSet ds_t = file.createDataSet("time", PredType::NATIVE_DOUBLE, space_t);
    ds_t.write(time_points.data(), PredType::NATIVE_DOUBLE);

    // 保存光强数组
    DataSet ds_p =
        file.createDataSet("intensity", PredType::NATIVE_DOUBLE, space_t);
    ds_p.write(power_values.data(), PredType::NATIVE_DOUBLE);

    // Save Polarization Attribute
    DataSpace attr_space(H5S_SCALAR);
    int rad_flag = RADIAL;
    Attribute attr_p =
        file.createAttribute("radial", PredType::NATIVE_INT, attr_space);
    attr_p.write(PredType::NATIVE_INT, &rad_flag);

  } catch (Exception &e) {
    cerr << "HDF5 Error: " << e.getDetailMsg() << endl;
    return -1;
  }

  cout << "导出成功！" << endl;

  // 4. 调用 Python 脚本进行画图
  cout << "正在调用 Python 脚本进行画图..." << endl;
  int result = system("python3 plot_temporal.py");
  if (result != 0) {
    cerr << "Warning: Plotting script failed with error code " << result
         << endl;
  }

  return 0;
}
