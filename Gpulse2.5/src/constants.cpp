/**
 * sctr - Constants Implementation
 * 常量和辅助函数实现
 */

#include "constants.h"

// --- 复数单位虚部 定义 ---
const Complex I(0, 1);

// --- 辅助函数实现 ---

/*
 * 生成线性间隔向量 (Linspace)
 * MATLAB style linspace function
 */
std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> ret(num);
    double step = (end - start) / (num - 1);
    for(int i = 0; i < num; ++i) {
        ret[i] = start + i * step;
    }
    return ret;
}

/**
 * 构建抛物面镜几何 (Build Geometry)
 * 核心功能：
 * 1. 在极坐标系 (r, theta) 下划分网格。
 * 2. 将极坐标映射到抛物面上的三维坐标 (x,y,z)。
 * 3. 计算每一点的法向量 n 和面积微元 dS。
 * 4. 应用 Simpson 积分法则计算权重。
 * 5. 应用 Mosaic 扇区偏振映射 (将线偏振转换为近似径向偏振)。
 */
Geometry build_geometry() {
    Geometry geom;
    // 计算 Simpson 积分所需的网格点数 (必须是奇数)
    int Nr = static_cast<int>(pow(2, NI)) + 1;
    int Nth = static_cast<int>(pow(2, NI + 1)) + 1;
    geom.total_points = Nr * Nth;

    // 分配内存
    geom.x.resize(geom.total_points);
    geom.y.resize(geom.total_points);
    geom.z.resize(geom.total_points);
    geom.Nx.resize(geom.total_points);
    geom.Ny.resize(geom.total_points);
    geom.Nz.resize(geom.total_points);
    geom.dS.resize(geom.total_points);
    geom.r.resize(geom.total_points);
    geom.theta.resize(geom.total_points);
    geom.psi_map.resize(geom.total_points);

    // 生成径向 (r) 和 角度 (theta) 的基础网格
    // r 从 HOLE/2 (内径) 到 D/2 (外径)
    std::vector<double> r_vec = linspace(HOLE / 2, D / 2, Nr);
    std::vector<double> th_vec = linspace(-PI, PI, Nth);
    double dr = r_vec[1] - r_vec[0];
    double dth = th_vec[1] - th_vec[0];

    // --- Simpson 积分权重生成 (1-4-2-4...-1) ---
    // 径向权重
    std::vector<double> wr(Nr, 1.0);
    for(int i = 1; i < Nr - 1; i += 2) wr[i] = 4;
    for(int i = 2; i < Nr - 2; i += 2) wr[i] = 2;

    // 角向权重
    std::vector<double> wth(Nth, 1.0);
    for(int i = 1; i < Nth - 1; i += 2) wth[i] = 4;
    for(int i = 2; i < Nth - 2; i += 2) wth[i] = 2;

    // --- 生成表面网格 ---
    int idx = 0;
    for(int j = 0; j < Nth; ++j) {
        for(int i = 0; i < Nr; ++i) {
            double r = r_vec[i], th = th_vec[j];
            geom.r[idx] = r;
            geom.theta[idx] = th;

            // 1. 抛物面坐标变换
            // 抛物面方程: z = r^2 / (4*f) - f
            // (焦点位于原点，顶点位于 -f)
            geom.x[idx] = r * cos(th);
            geom.y[idx] = r * sin(th);
            geom.z[idx] = r * r / (4 * F0) - F0;

            // 2. 计算法向量 (Normal Vector)
            // 抛物面法向量公式推导所得
            double dzdx = geom.x[idx] / (2 * F0);
            double dzdy = geom.y[idx] / (2 * F0);
            double Norm = sqrt(dzdx * dzdx + dzdy * dzdy + 1);
            
            // 指向焦点方向的法向量
            geom.Nx[idx] = -dzdx / Norm;
            geom.Ny[idx] = -dzdy / Norm;
            geom.Nz[idx] = 1.0 / Norm;

            // 3. 计算面积微元 (Area Element)
            // dS = r * dr * dtheta * sqrt(1 + (dz/dr)^2) ... 
            // 这里 Norm 即为度规因子。并将 Simpson 权重 (wr * wth / 9) 合并入 dS
            geom.dS[idx] = (r * dr * dth) * (wth[j] * wr[i] / 9.0) * Norm;

            // 4. 计算 Mosaic 偏振映射 (Mosaic Polarization Mapping)
            // 将反射镜分为 4 个象限 (Quadrants)，每个象限赋予不同的偏振旋转角 psi
            // 这种技术用于通过拼接波片将线偏振近似转换为径向偏振
            // Q1: psi=pi/4, Q2: psi=pi/2, Q3: psi=3pi/4, Q4: psi=0
            if (RADIAL) {
                double trot = arg(exp(I * (th + PI / 4.0))); // 旋转参考系
                
                if(trot >= -PI / 2 && trot < 0) {
                    geom.psi_map[idx] = 0;
                } else if(trot >= 0 && trot < PI / 2) {
                    geom.psi_map[idx] = PI / 4;
                } else if(trot >= PI / 2 || trot <= -PI) {
                    geom.psi_map[idx] = PI / 2;
                } else {
                    geom.psi_map[idx] = 3 * PI / 4;
                }
            } else {
                // 线偏振模式: 不进行偏振旋转
                geom.psi_map[idx] = 0;
            }
            idx++;
        }
    }
    return geom;
}

/**
 * 计算光谱分布 (Calculate Spectrum)
 * 生成超高斯 (Super-Gaussian) 频谱分布并归一化能量
 */
Spectrum calculate_spectrum() {
    Spectrum spec;

    // 1. 计算中心频率与带宽
    double w0 = 2 * PI * C / LAMBDA0; // 中心角频率
    // 将波长带宽转换为频率带宽近似
    double dw_Range = WIDTH * w0 * w0 / (2 * PI * C);
    double m = ORDER; // 超高斯阶数 (m=2 为普通高斯)
    
    // 计算特征宽度参数 (Width Parameter)
    double width_p = (dw_Range / 2.0) / pow(log(2.0), 1.0 / m);

    // 2. 生成频率列表
    spec.omega_list = linspace(w0 - dw_Range, w0 + dw_Range, N_FREQ);
    spec.d_omega = spec.omega_list[1] - spec.omega_list[0];

    // 3. 生成频域 Simpson 权重
    spec.weights.resize(N_FREQ);
    for(int i = 0; i < N_FREQ; ++i) {
        // 首尾点权重1，奇数点4，偶数点2，最后除以3
        double w_simp = (i == 0 || i == N_FREQ - 1) ? 1 : ((i % 2 == 1) ? 4 : 2);
        spec.weights[i] = (w_simp / 3.0) * spec.d_omega;
    }

    // 4. 计算振幅分布与能量归一化
    // 目标总能量 (Intensity -> E_field)
    double I_peak = 2 * P_PEAK / (PI * pow(W0, 2));
    double E_tgt = sqrt(2 * I_peak / (C * EPS0));
    
    double sum_num = 0;
    spec.E_amp.resize(N_FREQ);
    
    for(int i = 0; i < N_FREQ; ++i) {
        // 超高斯函数 A(w)
        double val = exp(-0.5 * pow(abs((spec.omega_list[i] - w0) / width_p), m));
        spec.E_amp[i] = val;
        
        // 累积积分值用于归一化
        sum_num += val * spec.weights[i];
    }

    // 执行归一化：使得所有频率分量积分后的峰值等于目标场强 E_tgt
    double ratio = E_tgt / sum_num;
    for(auto& v : spec.E_amp) {
        v *= ratio;
    }

    return spec;
}