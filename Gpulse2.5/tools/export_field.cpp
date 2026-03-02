/**
 * @file export_field.cpp
 * @brief 高斯激光脉冲场分布导出工具
 *
 * 该程序用于计算并导出高斯激光脉冲在 z=0 平面上的电场和磁场分布。
 * 通过调用 Gpulse 库获取激光场在二维网格上的数值，并将结果保存为 HDF5 格式。
 *
 * 输出文件: field_t0.h5 (包含 Ex, Ey, Ez, Bx, By, Bz 六个场分量)
 *
 * @note 依赖 Gpulse 库进行激光场计算
 * @note 使用 HDF5 C++ 库进行数据存储
 */

#include "Gpulse.h"
#include "constants.h" // 引入 RANGE 常量（空间扫描范围）
#include <H5Cpp.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace H5;

/**
 * @brief 导出网格分辨率
 *
 * 定义了输出网格的分辨率，以 100x100 的网格进行场分布采样。
 * 该值可以根据需要调整，以获得更高或更低的空间分辨率。
 */
const int EXPORT_RES = 100; // 100x100 二维网格

int main(int argc, char **argv) {
  cout << "启动激光场导出工具..." << endl;

  /**
   * @brief 第一步：准备数据容器
   *
   * 为电场和磁场的六个分量分配内存空间。
   * 现在存储复数分量 (Real, Imag) 以支持正确的投影计算。
   * 输出将被保存为 Ex_r, Ex_i 等 dataset
   */
  int N = EXPORT_RES;
  vector<double> data_Ex_r(N * N), data_Ex_i(N * N);
  vector<double> data_Ey_r(N * N), data_Ey_i(N * N);
  vector<double> data_Ez_r(N * N), data_Ez_i(N * N);

  vector<double> data_Bx_r(N * N), data_Bx_i(N * N);
  vector<double> data_By_r(N * N), data_By_i(N * N);
  vector<double> data_Bz_r(N * N), data_Bz_i(N * N);

  // ... (Grid setup same) ...
  double range = RANGE;
  cout << "扫描网格: " << N << "x" << N << " 覆盖范围 +/- " << range << " 米"
       << endl;
  auto linspace = [](double start, double end, int num) {
    vector<double> ret(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i)
      ret[i] = start + i * step;
    return ret;
  };
  vector<double> axis = linspace(-range, range, N);

  // ... (Args parsing same) ...
  double t = 0.0;
  if (argc > 1) {
    try {
      t = std::stod(argv[1]);
    } catch (...) {
    }
  }

  double z = 0.0;
  Gpulse_Init();

  // Loop
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      double x = axis[i];
      double y = axis[j];
      Complex Ex, Ey, Ez, Bx, By, Bz;

      GpulseComplex(t, x, y, z, Ex, Ey, Ez, Bx, By, Bz);

      int idx = i * N + j;
      data_Ex_r[idx] = Ex.real();
      data_Ex_i[idx] = Ex.imag();
      data_Ey_r[idx] = Ey.real();
      data_Ey_i[idx] = Ey.imag();
      data_Ez_r[idx] = Ez.real();
      data_Ez_i[idx] = Ez.imag();

      data_Bx_r[idx] = Bx.real();
      data_Bx_i[idx] = Bx.imag();
      data_By_r[idx] = By.real();
      data_By_i[idx] = By.imag();
      data_Bz_r[idx] = Bz.real();
      data_Bz_i[idx] = Bz.imag();
    }
    if (i % 10 == 0)
      cout << "进度: " << i << "/" << N << endl;
  }

  string filename = "field_snapshot.h5";
  cout << "保存到 " << filename << " (Complex)"
       << "..." << endl;

  try {
    H5File file(filename, H5F_ACC_TRUNC);
    hsize_t dims[2] = {static_cast<hsize_t>(N), static_cast<hsize_t>(N)};
    DataSpace dataspace(2, dims);

    auto save_ds = [&](const string &name, const vector<double> &data) {
      DataSet ds = file.createDataSet(name, PredType::NATIVE_DOUBLE, dataspace);
      ds.write(data.data(), PredType::NATIVE_DOUBLE);
    };

    save_ds("Ex_r", data_Ex_r);
    save_ds("Ex_i", data_Ex_i);
    save_ds("Ey_r", data_Ey_r);
    save_ds("Ey_i", data_Ey_i);
    save_ds("Ez_r", data_Ez_r);
    save_ds("Ez_i", data_Ez_i);

    save_ds("Bx_r", data_Bx_r);
    save_ds("Bx_i", data_Bx_i);
    save_ds("By_r", data_By_r);
    save_ds("By_i", data_By_i);
    save_ds("Bz_r", data_Bz_r);
    save_ds("Bz_i", data_Bz_i);

    // Save Range Attribute
    DataSpace attr_space(H5S_SCALAR);
    Attribute attr_r =
        file.createAttribute("range", PredType::NATIVE_DOUBLE, attr_space);
    attr_r.write(PredType::NATIVE_DOUBLE, &range);

    // Save Time Attribute
    Attribute attr_t =
        file.createAttribute("t", PredType::NATIVE_DOUBLE, attr_space);
    attr_t.write(PredType::NATIVE_DOUBLE, &t);

    // Save Polarization Attribute
    // Use int type for RADIAL flag
    int rad_flag = RADIAL;
    Attribute attr_p =
        file.createAttribute("radial", PredType::NATIVE_INT, attr_space);
    attr_p.write(PredType::NATIVE_INT, &rad_flag);

  } catch (FileIException error) {
    error.printErrorStack();
    return -1;
  } catch (DataSetIException error) {
    error.printErrorStack();
    return -1;
  } catch (DataSpaceIException error) {
    error.printErrorStack();
    return -1;
  }

  cout << "导出完成！" << endl;

  // Automatically call the plotting script
  cout << "Executing plot_field.py..." << endl;
  int ret_code = system("python3 plot_field.py");
  if (ret_code != 0) {
    cerr << "Error: plot_field.py failed to execute." << endl;
  }

  return 0;
}
