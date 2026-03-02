#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <complex>
#include <vector>

/**
 * 物理常数与模拟配置头文件
 * 此文件定义了模拟所需的所有核心常数、几何参数和光谱参数。
 */

// --- 基础物理常数 (Basic Physical Constants) ---
// --- 基础物理常数 (Basic Physical Constants) ---
constexpr double PI = 3.14159265358979323846;

// --- 类型定义 ---
using Complex = std::complex<double>;
extern const Complex I; // 虚数单位 i

// 模拟单位制选择 (Unit System Toggle)
// true: 使用归一化单位 (Normalized Units)
// false: 使用国际单位 (SI Units)
constexpr bool USE_NORMALIZED_UNITS = false;

// --- SI 基准值 (SI Reference Values) ---
constexpr double C_SI = 299792458.0;          // 光速 (m/s)
constexpr double EPS0_SI = 8.854187817e-12;   // 真空介电常数 (F/m)
constexpr double Me_SI = 9.10938356e-31;      // 电子质量 (kg)
constexpr double Qe_SI = 1.60217663e-19;      // 基本电荷 (C)

// --- 物理参数输入的 SI 值 (Input Parameters in SI) ---
constexpr double LAMBDA0_SI = 1.8e-6;   // 中心波长 (1.8 um)
constexpr double W0_SI = 8.9e-3;        // 束腰半径 (8.9 mm)
constexpr double P_PEAK_SI = 0.3e12;    // 峰值功率 (0.3 TW)
constexpr double WIDTH_SI = 700e-9;     // 带宽 (700 nm)
constexpr double F0_SI = 0.00635;       // 焦距 (6.35 mm)
constexpr double D_SI = 0.0254;         // 直径 (25.4 mm)
constexpr double HOLE_SI = 0.005;       // 孔径 (5 mm)
constexpr double RANGE_SI = 18e-6;      // 窗口半宽 (18 um)

// --- 归一化基准量计算 (Reference Scales Calculation) ---
// 基于 lambda0 对应的频率 omega_r 进行归一化
const double omega_r = 2.0 * PI * C_SI / LAMBDA0_SI;

const double T_r = 1.0 / omega_r;                // 时间基准
const double L_r = C_SI / omega_r;               // 长度基准
const double V_r = C_SI;                         // 速度基准
const double E_r = Me_SI * C_SI * omega_r / Qe_SI; // 电场基准
const double B_r = Me_SI * omega_r / Qe_SI;      // 磁场基准

// 导出基准: 强度与功率
// I_r = 0.5 * c * eps0 * E_r^2
const double I_r = 0.5 * C_SI * EPS0_SI * E_r * E_r;
// Power_r = I_r * L_r^2 (功率对面积的积分需考虑 L_r^2 因子)
const double Power_r = I_r * L_r * L_r;

// --- 实际使用的常数定义 (Active Constants) ---

// 根据选择的单位制设置常数
const double C = USE_NORMALIZED_UNITS ? 1.0 : C_SI;
// 在归一化单位下，设定 eps0=2.0 使得 E = sqrt(2*I) 的公式成立
// 推导: E_norm = E_SI / E_r -> C_norm * EPS_norm = 2
const double EPS0 = USE_NORMALIZED_UNITS ? 2.0 : EPS0_SI;

// --- 模拟参数 (Simulation Parameters) ---
// 根据单位制自动转换
const double LAMBDA0 = USE_NORMALIZED_UNITS ? (LAMBDA0_SI / L_r) : LAMBDA0_SI;
const double W0      = USE_NORMALIZED_UNITS ? (W0_SI / L_r)      : W0_SI;
const double P_PEAK  = USE_NORMALIZED_UNITS ? (P_PEAK_SI / Power_r) : P_PEAK_SI;
const double WIDTH   = USE_NORMALIZED_UNITS ? (WIDTH_SI / L_r)   : WIDTH_SI;
const double F0      = USE_NORMALIZED_UNITS ? (F0_SI / L_r)      : F0_SI;
const double D       = USE_NORMALIZED_UNITS ? (D_SI / L_r)       : D_SI;
const double HOLE    = USE_NORMALIZED_UNITS ? (HOLE_SI / L_r)    : HOLE_SI;
const double RANGE   = USE_NORMALIZED_UNITS ? (RANGE_SI / L_r)   : RANGE_SI;

// 这一项是无量纲的，保持不变
const int ORDER = 7;
const int RADIAL = 1;

// --- 模拟精度与网格参数 (Simulation Grid Parameters) ---
// 决定了计算的精度和内存消耗
 
const int NI = 8; // 几何积分的辛普森网格阶数 (2^NI)
const int N_FREQ = 101; // 频谱采样点数

// 焦平面观察窗口网格数
const int local_ycells = 100;
const int local_zcells = 100;
const int num_x = local_ycells;
const int num_y = local_zcells;

// --- 数据结构定义 ---

/**
 * 几何结构体 (Geometry)
 * 存储抛物面反射镜表面的离散化信息
 */
struct Geometry {
  std::vector<double> x, y, z;    // 反射镜表面网格点的坐标
  std::vector<double> Nx, Ny, Nz; // 表面单位法向量
  std::vector<double> dS;         // 面积微元 (已包含辛普森权重)
  std::vector<double> r, theta;   // 对应的极坐标
  std::vector<double> psi_map;    // 径向偏振的 Mosaic 映射角
  int total_points;               // 总积分点数
};

/**
 * 频谱结构体 (Spectrum)
 * 存储离散化的频域信息
 */
struct Spectrum {
  std::vector<double> omega_list; // 角频率列表
  std::vector<double> weights;    // 频域积分权重 (辛普森法)
  std::vector<double> E_amp;      // 对应频率的电场振幅
  double d_omega;                 // 频率间隔 d_omega
};

// --- 辅助函数声明 ---

// 线生成函数
std::vector<double> linspace(double start, double end, int num);

// 构建几何对象
Geometry build_geometry();

// 计算光谱分布
Spectrum calculate_spectrum();

#endif // CONSTANTS_H