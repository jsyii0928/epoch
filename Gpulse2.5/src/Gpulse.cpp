/**
 * @file Gpulse.cpp
 * @brief Gpulse 激光场求解器实现文件
 *
 * 功能描述：
 * 计算径向偏振光经过抛物面镜聚焦后的矢量电磁场分布。
 * 使用基于 Stratton-Chu 衍射积分的物理模型，支持超短脉冲宽光谱特性。
 *
 * 核心算法：
 * 1. 矢量衍射积分：计算完整的电磁矢量分量 (Ex, Ey, Ez, Bx, By, Bz)
 * 2. 宽光谱叠加：模拟飞秒脉冲的时空聚焦特性
 * 3. 内存驻留设计：预计算整个焦场的时空分布并驻留内存，实现快速查询
 * 4. 双线性插值：支持任意空间坐标 (x,y) 的连续查询
 */

#include "Gpulse.h"
#include "constants.h"
#include <cmath>
#include <complex>
#include <cstring> // For memcpy/memcmp
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <omp.h>
#include <vector>

using namespace std;

// ==========================================
// 1. 内部数据结构 (Internal Data Structures)
// ==========================================

/**
 * @class SpectralField
 * @brief 频谱场数据存储类
 *
 * 存储预计算的场数据，每个网格点存储多个频率的场分量。
 * 数据布局：[pixel_x][pixel_y][frequency] 的三维数组结构
 */
class SpectralField {
public:
  /**
   * @struct Cell
   * @brief 场数据单元格
   *
   * 存储单个网格点、单个频率下的所有电磁场分量
   */
  struct Cell {
    Complex Ex, Ey, Ez; // 电场分量
    Complex Bx, By, Bz; // 磁场分量
  };

  vector<Cell> data;      // 场数据数组 (展平的三维结构)
  vector<double> omegas;  // 频率列表 (每个存储点对应的角频率)

  // 网格参数
  int N_pixels_x;         // x 方向网格点数
  int N_pixels_y;         // y 方向网格点数
  int N_freq;             // 频率采样点数

  // 空间范围
  double min_x, max_x;    // x 方向最小/最大坐标
  double min_y, max_y;    // y 方向最小/最大坐标
  double dx, dy;          // x/y 方向网格间距

  /**
   * @brief 调整数据容器大小
   * @param n_pix_x x 方向网格点数
   * @param n_pix_y y 方向网格点数
   * @param n_freq 频率采样点数
   */
  void resize(int n_pix_x, int n_pix_y, int n_freq) {
    N_pixels_x = n_pix_x;
    N_pixels_y = n_pix_y;
    N_freq = n_freq;
    data.resize(N_pixels_x * N_pixels_y * N_freq);
    omegas.resize(N_freq);
  }

  /**
   * @brief 设置指定位置的场值
   * @param ix x 方向网格索引
   * @param iy y 方向网格索引
   * @param freq_idx 频率索引
   * @param ex, ey, ez 电场分量
   * @param bx, by, bz 磁场分量
   * @param weight 权重因子（积分权重）
   *
   * 索引映射：idx = (ix * N_pixels_y + iy) * N_freq + freq_idx
   */
  inline void set_at(int ix, int iy, int freq_idx, Complex ex, Complex ey,
                     Complex ez, Complex bx, Complex by, Complex bz,
                     double weight) {
    int pix_idx = ix * N_pixels_y + iy;
    int idx = pix_idx * N_freq + freq_idx;
    data[idx].Ex = ex * weight;
    data[idx].Ey = ey * weight;
    data[idx].Ez = ez * weight;
    data[idx].Bx = bx * weight;
    data[idx].By = by * weight;
    data[idx].Bz = bz * weight;
  }

  /**
   * @brief 获取指定索引和频率的时域场
   * @param ix x 方向网格索引
   * @param iy y 方向网格索引
   * @param t 模拟时间
   * @param z z 坐标（用于传播相位）
   * @param sum_E_x, sum_E_y, sum_E_z 输出：时域电场
   * @param sum_B_x, sum_B_y, sum_B_z 输出：时域磁场
   *
   * 通过逆傅里叶变换将频域场转换为时域场：
   * E(t,z) = Σ_n E_n(ω_n) * exp(-i*ω_n*(t - z/c))
   */
  inline void get_field_at_index_freq(int ix, int iy, double t, double z,
                                      Complex &sum_E_x, Complex &sum_E_y,
                                      Complex &sum_E_z, Complex &sum_B_x,
                                      Complex &sum_B_y,
                                      Complex &sum_B_z) const {
    int pix_idx = ix * N_pixels_y + iy;
    int base_idx = pix_idx * N_freq;

    // 遍历所有频率分量，叠加得到时域场
    for (int k = 0; k < N_freq; ++k) {
      // 时间相位因子：ω*t
      double theta_t = omegas[k] * t;
      // 传播相位因子：ω*z/c
      double theta_z = omegas[k] * z / C;
      // 总相位：θ = ω*(t - z/c)
      // 这里的相位因子对应 exp(-i*θ)，因为物理场是 e^(-iωt) 的叠加
      double theta = theta_t + theta_z;
      Complex phase_factor(cos(theta), -sin(theta));

      const Cell &cell = data[base_idx + k];

      // 叠加每个频率分量的贡献
      sum_E_x += cell.Ex * phase_factor;
      sum_E_y += cell.Ey * phase_factor;
      sum_E_z += cell.Ez * phase_factor;
      sum_B_x += cell.Bx * phase_factor;
      sum_B_y += cell.By * phase_factor;
      sum_B_z += cell.Bz * phase_factor;
    }
  }
};

// 全局场数据指针（单例模式）
static unique_ptr<SpectralField> g_field;
// 初始化标志（用于 std::call_once 确保只初始化一次）
static std::once_flag g_init_flag;

/**
 * @brief 将角度映射到 (-π, π] 区间
 * @param angle 输入角度（任意值）
 * @return 映射后的角度
 *
 * 用于角度的规范化处理，确保角度在标准区间内
 */
double wrapToPi(double angle) {
  angle = fmod(angle + PI, 2 * PI);
  if (angle < 0)
    angle += 2 * PI;
  return angle - PI;
}

// ==========================================
// 2. 求解器核心逻辑 (Refactored Solver Logic)
// ==========================================

/**
 * @brief 运行主求解器，计算并存储焦场数据
 *
 * 执行步骤：
 * 1. 生成反射镜表面网格（抛物面）
 * 2. 计算超高频谱分布
 * 3. 对每个频率执行 Stratton-Chu 矢量衍射积分
 * 4. 存储结果到内存供快速查询
 */
void run_solver() {
  // 创建场数据容器
  g_field = make_unique<SpectralField>();

  // 从 constants.h 读取网格参数
  int N_x = num_x;          // x 方向网格数
  int N_y = num_y;          // y 方向网格数
  double range_val = RANGE;  // 观察窗口半宽

  // 初始化场容器
  g_field->resize(N_x, N_y, N_FREQ);
  g_field->min_x = -range_val;
  g_field->max_x = range_val;
  g_field->min_y = -range_val;
  g_field->max_y = range_val;
  g_field->dx = (2.0 * range_val) / (N_x - 1);
  g_field->dy = (2.0 * range_val) / (N_y - 1);

  std::cout << "[Gpulse] 正在初始化场解算器 (Initializing Field Solver - "
               "Linearized)..."
            << std::endl;

  // --- 1. 源网格生成 (Source Grid Generation) ---
  // 在抛物面镜表面生成积分网格
  // 复制自 radical_debug.m / radical_debug.cpp 逻辑
  int Nr = (1 << NI) + 1;       // 径向网格点数 (2^NI + 1)
  int Nth = (1 << (NI + 1)) + 1; // 角向网格点数 (2^(NI+1) + 1)

  // 径向坐标向量 r
  vector<double> r_vec(Nr);
  // 径向间距：从内孔边缘到外径边缘
  double dr =
      (D / 2.0 - HOLE / 2.0) / (Nr - 1);
  for (int i = 0; i < Nr; ++i)
    r_vec[i] = HOLE / 2.0 + i * dr;

  // 角度向量 theta
  vector<double> th_vec(Nth);
  double dth = (2 * PI) / (Nth - 1); // -PI 到 PI
  for (int i = 0; i < Nth; ++i)
    th_vec[i] = -PI + i * dth;

  int total_src_points = Nr * Nth; // 总源点数

  // 源几何数组（存储反射镜表面的离散点信息）
  vector<double> src_x(total_src_points);     // x 坐标
  vector<double> src_y(total_src_points);     // y 坐标
  vector<double> src_z(total_src_points);     // z 坐标
  vector<double> src_Nx(total_src_points);    // 法向量 x 分量
  vector<double> src_Ny(total_src_points);    // 法向量 y 分量
  vector<double> src_Nz(total_src_points);    // 法向量 z 分量
  vector<double> src_dS(total_src_points);   // 面积微元（含辛普森权重）
  vector<double> src_r(total_src_points);     // 极坐标半径

  // 偏振映射系数（用于 Mosaic 相位板）
  vector<double> cos_2psi(total_src_points);
  vector<double> sin_2psi(total_src_points);

  // --- 网格构建循环 ---
  int idx = 0;
  for (int i = 0; i < Nth; ++i) {
    double th = th_vec[i];

    // 辛普森积分权重 - theta 方向
    // 权重模式：1, 4, 2, 4, 2, ..., 4, 1
    double w_th = 1.0;
    if (i > 0 && i < Nth - 1)
      w_th = (i % 2 == 1) ? 4.0 : 2.0;

    for (int j = 0; j < Nr; ++j) {
      double r = r_vec[j];

      src_r[idx] = r;

      // --- 坐标计算 ---
      // 抛物面方程：z = r²/(4f) - f
      // 使用极坐标转笛卡尔坐标
      src_x[idx] = r * cos(th);
      src_y[idx] = r * sin(th);
      src_z[idx] = (r * r) / (4 * F0) - F0;

      // --- 几何因子计算 ---
      // 计算抛物面的法向量
      // 从抛物面方程 z = r²/(4f) - f 推导：
      // dz/dx = x/(2f), dz/dy = y/(2f)
      double dz_dx = src_x[idx] / (2 * F0);
      double dz_dy = src_y[idx] / (2 * F0);
      // 法向量归一化因子
      double norm_N = sqrt(dz_dx * dz_dx + dz_dy * dz_dy + 1.0);

      // 指向焦点的单位法向量
      src_Nx[idx] = -dz_dx / norm_N;
      src_Ny[idx] = -dz_dy / norm_N;
      src_Nz[idx] = 1.0 / norm_N;

      // 辛普森积分权重 - r 方向
      double w_r = 1.0;
      if (j > 0 && j < Nr - 1)
        w_r = (j % 2 == 1) ? 4.0 : 2.0;

      // --- 有效面积微元 ---
      // dS = r * dr * dtheta * sqrt(1 + (dz/dr)² + ...) * 辛普森权重/9
      double W_simpson = w_th * w_r;
      src_dS[idx] = (r * dr * dth / 9.0) * W_simpson * norm_N;

      // --- Mosaic 偏振映射逻辑 ---
      // 用于将线偏振转换为近似径向偏振
      // 将镜面分为 4 个象限，每个象限赋予不同的偏振旋转角 psi
      double psi = 0.0;
      if (RADIAL) { // 径向偏振模式
        // 旋转参考系（π/4 偏移使象限边界在 45° 方向）
        double theta_rot = wrapToPi(th + PI / 4.0);
        double aa = PI / 4.0; // 象限旋转角

        // 四象限映射（Mosaic 相位板近似）
        if (theta_rot >= 0 && theta_rot < PI / 2.0) {
          psi = aa;          // 第一象限：π/4
        } else if (theta_rot >= PI / 2.0 && theta_rot <= PI) {
          psi = 2 * aa;     // 第二象限：π/2
        } else if (theta_rot >= -PI && theta_rot < -PI / 2.0) {
          psi = 3 * aa;     // 第三象限：3π/4
        } else {
          psi = 0;          // 第四象限：0
        }
      } else {
        psi = 0; // 线偏振模式：无旋转
      }
      // 预计算偏振系数
      cos_2psi[idx] = cos(2 * psi);
      sin_2psi[idx] = sin(2 * psi);

      idx++;
    }
  }

  // --- 2. 频谱计算 (Spectrum Calculation) ---
  // 生成超高频谱分布并归一化能量
  double omega_0 = 2 * PI * C / LAMBDA0; // 中心角频率
  double delta_lambda = WIDTH;             // 波长带宽 (FWHM)
  // 将波长带宽转换为角频率带宽
  // Δω ≈ (2πc/λ²) * Δλ
  double delta_omega = delta_lambda * omega_0 * omega_0 / (2 * PI * C);
  int super_gauss_order = ORDER;          // 超高斯阶数
  // 计算超高斯宽度参数
  double scaling_factor = pow(log(2), 1.0 / super_gauss_order);
  double width_param = (delta_omega / 2.0) / scaling_factor;

  // 频率采样范围：±Δω
  double w_range = 1 * delta_omega;
  int N_SAMPLES = N_FREQ;
  vector<double> omega_list(N_SAMPLES);
  double d_omega = (2 * w_range) / (N_SAMPLES - 1);
  // 生成频率列表（均匀采样）
  for (int i = 0; i < N_SAMPLES; ++i)
    omega_list[i] = (omega_0 - w_range) + i * d_omega;

  // 将频率列表存储到 g_field（用于时域重构）
  for (int k = 0; k < N_SAMPLES; ++k)
    g_field->omegas[k] = omega_list[k];

  // --- 计算频域振幅 ---
  vector<double> A_vals(N_SAMPLES);  // 超高斯振幅
  vector<double> W_sim(N_SAMPLES, 1.0); // 辛普森积分权重

  for (int i = 0; i < N_SAMPLES; ++i) {
    double w = omega_list[i];
    // 超高斯函数：A(w) = exp(-0.5 * |(w-ω0)/σ|^m)
    A_vals[i] = exp(
        -0.5 * pow(std::abs((w - omega_0) / width_param), super_gauss_order));

    // 辛普森权重：首尾点为1，奇数点为4，偶数点为2，最后除以3
    if (i > 0 && i < N_SAMPLES - 1)
      W_sim[i] = (i % 2 == 1) ? 4.0 : 2.0;
    W_sim[i] /= 3.0;
  }

  // 计算积分后的总场强（未归一化）
  double E_s = 0.0;
  for (int i = 0; i < N_SAMPLES; ++i)
    E_s += A_vals[i] * d_omega * W_sim[i];

  // 计算目标峰值场强
  // 从峰值功率计算峰值光强：I_peak = 2*P/(π*w0²)
  double I0_peak = 2 * P_PEAK / (PI * W0 * W0);
  const double N_MED = 1.0; // 介质折射率（真空为1）
  // 峰值电场：E = sqrt(2*I/(n*c*ε0))
  double E_peak_target = sqrt(2 * I0_peak / (N_MED * C * EPS0));

  // 计算归一化因子并应用
  double Rat = E_peak_target / E_s;
  vector<double> E_n_list(N_SAMPLES);
  for (int i = 0; i < N_SAMPLES; ++i)
    E_n_list[i] = Rat * A_vals[i];

  // --- 3. 主积分循环 (Main Integration Loop) ---

  // --- 生成目标网格（展平用于 OpenMP 并行）---
  // 手动生成线性网格以匹配原代码逻辑
  vector<double> xf(N_x * N_y), yf(N_x * N_y);
  // 生成 x 方向均匀网格
  vector<double> lin_x =
      linspace(-range_val, range_val, N_x);

  // 展平网格坐标存储（便于并行处理）
  for (int i = 0; i < N_x; ++i) {
    double x_val = -range_val + i * (2 * range_val) / (N_x - 1);
    for (int j = 0; j < N_y; ++j) {
      double y_val = -range_val + j * (2 * range_val) / (N_y - 1);
      xf[i * N_y + j] = x_val;
      yf[i * N_y + j] = y_val;
    }
  }

  // 源平面到焦平面的传播时间（用于参考相位）
  double t0 = 2 * F0 / C;

  // --- 遍历所有频率 ---
  for (int n = 0; n < N_SAMPLES; ++n) {
    double wn = omega_list[n];      // 当前角频率
    double lam_n = 2 * PI * C / wn; // 对应波长
    double kn = 2 * PI * N_MED / lam_n; // 波数
    double E0_curr = E_n_list[n];       // 当前频率的电场振幅

    // 传播相位因子（包含源平面相位和时间参考相位）
    Complex time_comp_phase = exp(-I * wn * t0);
    Complex common_phase = exp(-I * PI / 2.0) * time_comp_phase;

    // 调试输出（仅在第 50 个频率）
    if (n == 50) {
#pragma omp critical
      {
        std::cout << "\n[DEBUG Gpulse] Freq 50" << std::endl;
        std::cout << "[DEBUG Gpulse] wn = " << wn << std::endl;
        std::cout << "[DEBUG Gpulse] kn = " << kn << std::endl;
        std::cout << "[DEBUG Gpulse] E0_curr = " << E0_curr << std::endl;
        std::cout << "[DEBUG Gpulse] common = " << common_phase << std::endl;
      }
    }

    // --- 优化：预计算源场项（每个频率只计算一次）---
    // 这避免了在每个像素点重新计算相同的源场项
    // 显著提高了计算效率

    int total_pixels = N_x * N_y;

    // --- 预计算源场项 ---
    // 这些项对所有观察点都是相同的，因此只需计算一次
    vector<Complex> vec_NxB_x(total_src_points);
    vector<Complex> vec_NxB_y(total_src_points);
    vector<Complex> vec_NxB_z(total_src_points);
    vector<Complex> vec_N_dot_E(total_src_points);

#pragma omp parallel for
    for (int s = 0; s < total_src_points; ++s) {
      double r = src_r[s];
      // 高斯光束相位因子：exp(-(r/w0)²)
      Complex gaussian_phase = exp(-pow(r / W0, 2));
      // 总相位因子：包含源平面传播和参考相位
      Complex phase_factor = exp(-I * kn * src_z[s]) * common_phase;
      // 入射场标量振幅
      Complex E_scalar = E0_curr * gaussian_phase * phase_factor;

      // --- 根据偏振模式计算入射场 ---
      if (RADIAL) {
        // Mosaic 近似的径向偏振
        // Ex = E0 * sin(2ψ), Ey = -E0 * cos(2ψ)
        Complex Ex = E_scalar * sin_2psi[s];
        Complex Ey = -E_scalar * cos_2psi[s];
        Complex Ez = 0.0;

        // 磁场：B = E/c，偏振旋转 90°
        Complex B_scalar = E_scalar / C;
        Complex Bx = B_scalar * cos_2psi[s];
        Complex By = B_scalar * sin_2psi[s];
        Complex Bz = 0.0;

        // 计算法向量与磁场的叉积 N × B
        double nx = src_Nx[s];
        double ny = src_Ny[s];
        double nz = src_Nz[s];

        vec_NxB_x[s] = ny * Bz - nz * By;
        vec_NxB_y[s] = nz * Bx - nx * Bz;
        vec_NxB_z[s] = nx * By - ny * Bx;
        // 计算法向量与电场的点积 N · E
        vec_N_dot_E[s] = nx * Ex + ny * Ey + nz * Ez;
      } else {
        // 线偏振模式（y 方向偏振）
        Complex Ex = 0.0;
        Complex Ey = E_scalar;
        Complex Ez = 0.0;

        Complex Bx = -E_scalar / C;
        Complex By = 0.0;
        Complex Bz = 0.0;

        double nx = src_Nx[s];
        double ny = src_Ny[s];
        double nz = src_Nz[s];

        vec_NxB_x[s] = ny * Bz - nz * By;
        vec_NxB_y[s] = nz * Bx - nx * Bz;
        vec_NxB_z[s] = nx * By - ny * Bx;
        vec_N_dot_E[s] = nx * Ex + ny * Ey + nz * Ez;
      }
    }

    // --- 对所有观察点进行积分 ---
#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < total_pixels; ++k) {
      double obs_x = xf[k]; // 观察点 x 坐标
      double obs_y = yf[k]; // 观察点 y 坐标
      double obs_z = 0.0;  // 观察点 z 坐标（焦平面）

      Complex sEx = 0, sEy = 0, sEz = 0; // 电场累加器
      Complex sBx = 0, sBy = 0, sBz = 0; // 磁场累加器

      // --- Stratton-Chu 矢量衍射积分 ---
      // 遍历所有源点，叠加贡献
      for (int s = 0; s < total_src_points; ++s) {
        // 从源点到观察点的距离矢量
        double key_Rx = obs_x - src_x[s];
        double key_Ry = obs_y - src_y[s];
        double key_Rz = obs_z - src_z[s];
        double R_sq = key_Rx * key_Rx + key_Ry * key_Ry + key_Rz * key_Rz;
        double R_norm = sqrt(R_sq);

        // --- 格林函数及其导数 ---
        // 格林函数：G = exp(ikR)/R
        Complex G = exp(I * kn * R_norm) / R_norm;
        // 格林函数梯度：∇G = (ik/R - 1/R²) * G
        Complex dG = (I * kn / R_norm - 1.0 / R_sq) * G;

        Complex nxe = vec_N_dot_E[s]; // 法向电场分量

        // --- Stratton-Chu 电场公式 ---
        // E = (iωμ₀)∮(n×B)G ds - (1/ε₀)∮(n·E)∇G ds
        // 这里使用简化形式，合并了常数因子
        sEx += (I * wn * vec_NxB_x[s] * G + nxe * dG * key_Rx) * src_dS[s];
        sEy += (I * wn * vec_NxB_y[s] * G + nxe * dG * key_Ry) * src_dS[s];
        sEz += (I * wn * vec_NxB_z[s] * G + nxe * dG * key_Rz) * src_dS[s];

        // --- Stratton-Chu 磁场公式 ---
        // B = ∮[(n×B)×∇G] ds
        // 展开叉积：(n×B)×∇G = (n×B)×(dG·R/R)
        sBx += (vec_NxB_y[s] * (dG * key_Rz) - vec_NxB_z[s] * (dG * key_Ry)) *
               src_dS[s];
        sBy += (vec_NxB_z[s] * (dG * key_Rx) - vec_NxB_x[s] * (dG * key_Rz)) *
               src_dS[s];
        sBz += (vec_NxB_x[s] * (dG * key_Ry) - vec_NxB_y[s] * (dG * key_Rx)) *
               src_dS[s];
      }

      // --- 应用积分常数因子和频率权重 ---
      double factor = 1.0 / (2 * PI);           // Stratton-Chu 积分常数
      double freq_weight = W_sim[n] * d_omega; // 当前频率的积分权重

      // 将结果存储到场数据容器
      // 权重 = factor * freq_weight
      int ix = k / N_y;
      int iy = k % N_y;

      g_field->set_at(ix, iy, n, sEx, sEy, sEz, sBx, sBy, sBz,
                      factor * freq_weight);
    }

    // 显示进度
    std::cout << "\r[Gpulse] Calculated " << n + 1 << "/" << N_SAMPLES
              << " freqs" << std::flush;
  }
  std::cout << std::endl
            << "[Gpulse] 初始化完成 (Initialization Complete)." << std::endl;
}

// ==========================================
// 3. 缓存与公共 API (Caching & Public API)
// ==========================================

/**
 * @struct CacheHeader
 * @brief 缓存文件头结构
 *
 * 用于验证缓存文件是否与当前参数匹配
 * 如果参数改变，缓存将失效并重新计算
 */
struct CacheHeader {
  // 网格参数
  int n_x, n_y, n_freq;
  double range;

  // 激光参数
  double lambda0, w0, p_peak;
  double width;
  int order;

  // 聚焦系统参数
  double f0, d, hole;
  int radial;

  /**
   * @brief 构造函数：从当前参数初始化头
   */
  CacheHeader() {
    n_x = num_x;
    n_y = num_y;
    n_freq = N_FREQ;
    range = RANGE;
    lambda0 = LAMBDA0;
    w0 = W0;
    p_peak = P_PEAK;
    width = WIDTH;
    order = ORDER;
    f0 = F0;
    d = D;
    hole = HOLE;
    radial = RADIAL;
  }

  /**
   * @brief 比较运算符：检查两个头是否相同
   * @param another 另一个 CacheHeader
   * @return 是否相同
   */
  bool operator==(const CacheHeader &other) const {
    return n_x == other.n_x && n_y == other.n_y && n_freq == other.n_freq &&
           std::abs(range - other.range) < 1e-12 &&
           std::abs(lambda0 - other.lambda0) < 1e-12 &&
           std::abs(w0 - other.w0) < 1e-12 &&
           std::abs(p_peak - other.p_peak) < 1e-12 &&
           std::abs(width - other.width) < 1e-12 && order == other.order &&
           std::abs(f0 - other.f0) < 1e-12 && std::abs(d - other.d) < 1e-12 &&
           std::abs(hole - other.hole) < 1e-12 && radial == other.radial;
  }
  bool operator!=(const CacheHeader &other) const { return !(*this == other); }
};

// 缓存文件名
const std::string CACHE_FILENAME = ".gpulse_cache.bin";

/**
 * @brief 保存场数据到缓存文件
 *
 * 将预计算的场数据写入二进制文件
 * 文件格式：[CacheHeader][omegas][field_data]
 */
void save_cache() {
  if (!g_field)
    return;
  std::ofstream out(CACHE_FILENAME, std::ios::binary);
  if (!out)
    return;
  // 写入头
  CacheHeader header;
  out.write(reinterpret_cast<const char *>(&header), sizeof(CacheHeader));
  // 写入频率列表
  out.write(reinterpret_cast<const char *>(g_field->omegas.data()),
            g_field->omegas.size() * sizeof(double));
  // 写入场数据
  out.write(reinterpret_cast<const char *>(g_field->data.data()),
            g_field->data.size() * sizeof(SpectralField::Cell));
  std::cout << "[Gpulse] Field data saved to cache: " << CACHE_FILENAME
            << std::endl;
}

/**
 * @brief 从缓存文件加载场数据
 * @return 是否加载成功
 *
 * 检查缓存文件是否存在且参数匹配
 * 如果参数不匹配，返回 false 需要重新计算
 */
bool load_cache() {
  std::ifstream in(CACHE_FILENAME, std::ios::binary);
  if (!in)
    return false;

  // 读取并验证头
  CacheHeader saved;
  in.read(reinterpret_cast<char *>(&saved), sizeof(CacheHeader));
  if (in.gcount() != sizeof(CacheHeader))
    return false;

  // 检查参数是否匹配
  CacheHeader current;
  if (saved != current) {
    std::cout << "[Gpulse] Cache parameters mismatch. Recalculating..."
              << std::endl;
    return false;
  }

  // 创建并填充场容器
  g_field = make_unique<SpectralField>();
  g_field->resize(saved.n_x, saved.n_y, saved.n_freq);
  g_field->min_x = -saved.range;
  g_field->max_x = saved.range;
  g_field->min_y = -saved.range;
  g_field->max_y = saved.range;
  g_field->dx = (2.0 * saved.range) / (saved.n_x - 1);
  g_field->dy = (2.0 * saved.range) / (saved.n_y - 1);

  // 读取频率列表
  in.read(reinterpret_cast<char *>(g_field->omegas.data()),
          g_field->omegas.size() * sizeof(double));
  // 读取场数据
  in.read(reinterpret_cast<char *>(g_field->data.data()),
          g_field->data.size() * sizeof(SpectralField::Cell));

  if (!in) {
    g_field.reset();
    return false;
  }
  std::cout << "[Gpulse] Loaded field data from cache." << std::endl;
  return true;
}

/**
 * @brief 初始化实现（带缓存）
 *
 * 尝试从缓存加载，失败则运行求解器并保存缓存
 */
void Initialize_Impl() {
  if (load_cache())
    return;
  run_solver();
  save_cache();
}

/**
 * @brief 双线性插值函数
 * @param c00, c10, c01, c11 四个角点的复数值
 * @param tx x 方向插值权重 (0~1)
 * @param ty y 方向插值权重 (0~1)
 * @return 插值结果
 *
 * 标准 2D 双线性插值公式
 */
Complex bilinear_interp(Complex c00, Complex c10, Complex c01, Complex c11,
                        double tx, double ty) {
  // 先沿 x 方向插值
  Complex c0 = c00 * (1.0 - tx) + c10 * tx;
  Complex c1 = c01 * (1.0 - tx) + c11 * tx;
  // 再沿 y 方向插值
  return c0 * (1.0 - ty) + c1 * ty;
}

// ==========================================
// 4. 公共 API 实现 (Public API Implementation)
// ==========================================

/**
 * @brief 初始化 Gpulse 求解器
 *
 * 使用 std::call_once 确保只初始化一次（线程安全）
 * 首次调用时自动执行初始化，后续调用被忽略
 */
void Gpulse_Init() { std::call_once(g_init_flag, Initialize_Impl); }

/**
 * @brief 释放 Gpulse 求解器占用的内存
 *
 * 清空预计算的场数据，释放内存
 * 调用后需要重新初始化才能使用
 */
void Gpulse_Free() { g_field.reset(); }

/**
 * @brief 获取指定位置和时间的实数场值
 * @param t 时间 (秒)，t=0 对应脉冲峰值到达焦平面
 * @param x x 坐标 (米)
 * @param y y 坐标 (米)
 * @param z z 坐标 (米)，在焦深范围内使用 2.5D 近似
 * @param Ex, Ey, Ez 输出：电场分量 (V/m)
 * @param Bx, By, Bz 输出：磁场分量 (Tesla)
 *
 * 返回场的实部，适用于需要瞬时场值的场景
 */
void Gpulse(double t, double x, double y, double z, double &Ex, double &Ey,
            double &Ez, double &Bx, double &By, double &Bz) {
  Complex cEx, cEy, cEz;
  Complex cBx, cBy, cBz;

  // 调用复数版本获取场
  GpulseComplex(t, x, y, z, cEx, cEy, cEz, cBx, cBy, cBz);

  // 输出实部
  Ex = cEx.real();
  Ey = cEy.real();
  Ez = cEz.real();
  Bx = cBx.real();
  By = cBy.real();
  Bz = cBz.real();
}

/**
 * @brief 获取指定位置和时间的复数场值
 * @param t 时间 (秒)
 * @param x x 坐标 (米)
 * @param y y 坐标 (米)
 * @param z z 坐标 (米)
 * @param Ex, Ey, Ez 输出：复数电场分量
 * @param Bx, By, Bz 输出：复数磁场分量
 *
 * 返回复数包络场，适用于投影计算等需要相位信息的场景
 *
 * 实现步骤：
 * 1. 检查坐标是否在网格范围内
 * 2. 使用双线性插值获取四个角点的场
 * 3. 叠加频率分量得到时域场
 * 4. 对四个角点进行双线性插值得到最终结果
 */
void GpulseComplex(double t, double x, double y, double z, Complex &Ex,
                   Complex &Ey, Complex &Ez, Complex &Bx, Complex &By,
                   Complex &Bz) {
  // 确保已初始化
  Gpulse_Init();

  // 边界检查：超出网格范围返回零场
  if (x < g_field->min_x || x > g_field->max_x || y < g_field->min_y ||
      y > g_field->max_y) {
    Ex = Ey = Ez = Bx = By = Bz = 0.0;
    return;
  }

  // 将物理坐标转换为网格索引
  double u = (x - g_field->min_x) / g_field->dx;
  double v = (y - g_field->min_y) / g_field->dy;
  int ix = static_cast<int>(floor(u)); // 左下角 x 索引
  int iy = static_cast<int>(floor(v)); // 左下角 y 索引

  // 边界钳制（防止越界）
  if (ix < 0)
    ix = 0;
  if (ix >= g_field->N_pixels_x - 1)
    ix = g_field->N_pixels_x - 2;
  if (iy < 0)
    iy = 0;
  if (iy >= g_field->N_pixels_y - 1)
    iy = g_field->N_pixels_y - 2;

  // 计算插值权重
  double tx = u - ix; // x 方向插值权重
  double ty = v - iy; // y 方向插值权重

  // 存储四个角点的场值
  Complex E_corners[2][2][3]; // [dx_i][dy_i][component]
  Complex B_corners[2][2][3];

  // 获取四个角点的时域场
  for (int dx_i = 0; dx_i <= 1; ++dx_i) {
    for (int dy_i = 0; dy_i <= 1; ++dy_i) {
      int cx = ix + dx_i;
      int cy = iy + dy_i;
      Complex Ex_c = 0, Ey_c = 0, Ez_c = 0;
      Complex Bx_c = 0, By_c = 0, Bz_c = 0;
      // 叠加所有频率分量得到时域场
      g_field->get_field_at_index_freq(cx, cy, t, z, Ex_c, Ey_c, Ez_c, Bx_c,
                                       By_c, Bz_c);
      E_corners[dx_i][dy_i][0] = Ex_c;
      E_corners[dx_i][dy_i][1] = Ey_c;
      E_corners[dx_i][dy_i][2] = Ez_c;
      B_corners[dx_i][dy_i][0] = Bx_c;
      B_corners[dx_i][dy_i][1] = By_c;
      B_corners[dx_i][dy_i][2] = Bz_c;
    }
  }

  // 双线性插值得到最终结果
  Ex = bilinear_interp(E_corners[0][0][0], E_corners[1][0][0],
                       E_corners[0][1][0], E_corners[1][1][0], tx, ty);
  Ey = bilinear_interp(E_corners[0][0][1], E_corners[1][0][1],
                       E_corners[0][1][1], E_corners[1][1][1], tx, ty);
  Ez = bilinear_interp(E_corners[0][0][2], E_corners[1][0][2],
                       E_corners[0][1][2], E_corners[1][1][2], tx, ty);
  Bx = bilinear_interp(B_corners[0][0][0], B_corners[1][0][0],
                       B_corners[0][1][0], B_corners[1][1][0], tx, ty);
  By = bilinear_interp(B_corners[0][0][1], B_corners[1][0][1],
                       B_corners[0][1][1], B_corners[1][1][1], tx, ty);
  Bz = bilinear_interp(B_corners[0][0][2], B_corners[1][0][2],
                       B_corners[0][1][2], B_corners[1][1][2], tx, ty);
}
