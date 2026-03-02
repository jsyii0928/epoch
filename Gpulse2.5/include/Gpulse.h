#ifndef GPULSE_H
#define GPULSE_H

/**
 * Gpulse - 紧聚焦脉冲场求解器 (Tight Focusing Pulse Solver)
 * 
 * 模块功能：
 * 计算径向偏振光 (Radially Polarized Light) 经过抛物面镜聚焦后的矢量电磁场。
 * 使用基于Stratton-Chu衍射积分的物理模型，支持超短脉冲宽光谱特性。
 * 
 * 核心算法特性：
 * 1. 矢量衍射积分：考虑完整的电磁矢量特性 (Ex, Ey, Ez, Bx, By, Bz)。
 * 2. 宽光谱叠加：模拟飞秒脉冲的时空聚焦特性。
 * 3. 内存驻留设计：预计算整个焦场的时空分布并驻留内存，查询速度极快。
 * 4. 双线性插值：支持任意空间坐标 (x,y) 的连续查询。
 */

/**
 * 初始化 Gpulse 求解器。
 * 
 * 功能：
 * 1. 构建聚焦几何结构 (抛物面镜)。
 * 2. 计算超高斯频谱分布。
 * 3. 执行Stratton-Chu矢量积分，生成焦场数据表。
 * 
 * 注意：
 * - 首次调用 Gpulse() 时会自动触发此函数，无需手动在该函数前调用。
 * - 只有在需要显式控制初始化时机（如预加载避免运行时卡顿）时使用。
 * - 该过程可能耗时数秒 (取决于 spectral 和 grid 分辨率)。
 */
void Gpulse_Init();

/**
 * 释放 Gpulse 求解器占用的内存。
 * 
 * 功能：
 * 清空内部预计算的场数据表 (g_field)，释放数百MB级别的内存。
 */
void Gpulse_Free();

/**
 * 获取特定时刻和位置的电磁场矢量。
 * 
 * 核心查询函数。支持任意坐标查询，内部会自动进行网格插值。
 * 
 * @param t  模拟时间 (秒). t=0 对应脉冲峰值到达焦平面时刻。
 * @param x  x 坐标 (米). 焦平面坐标系。
 * @param y  y 坐标 (米). 焦平面坐标系。
 * @param z  z 坐标 (米). 
 *           注意: 当前版本算法采用 2.5D 近似：
 *           假设场在焦深范围内横向分布形状不变，仅附加传播相位 exp(i(kz - \omega t))。
 *           在焦深(Rayleigh range)范围内该近似是准确的。
 * 
 * @param Ex [输出] x方向电场分量 (V/m)
 * @param Ey [输出] y方向电场分量 (V/m)
 * @param Ez [输出] z方向电场分量 (V/m)
 * @param Bx [输出] x方向磁场分量 (Tesla)
 * @param By [输出] y方向磁场分量 (Tesla)
 * @param Bz [输出] z方向磁场分量 (Tesla)
 */
// This header doesn't know about constants.h by default if only included by tools.
// But we need Complex type.
#include "constants.h" 

// Calculate fields at t, x, y, z (Returns Amplitude/Envelope)
void Gpulse(double t, double x, double y, double z,
            double& Ex, double& Ey, double& Ez,
            double& Bx, double& By, double& Bz);

// Calculate fields at t, x, y, z (Returns Complex Field e.g. for Projection)
void GpulseComplex(double t, double x, double y, double z,
                   Complex& Ex, Complex& Ey, Complex& Ez,
                   Complex& Bx, Complex& By, Complex& Bz);

#endif // GPULSE_H
