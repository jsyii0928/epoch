#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "Gpulse.h"
#include "constants.h"

using namespace std;

// 将指定网格内的六个场分量写入原始二进制文件 (.dat)
void dump_field_to_binary(double t, 
                          int nx_global, int ny_global, int nz_global,
                          double x_min, double x_max,
                          double y_min, double y_max,
                          double z_min, double z_max) {
    cout << "\nDumping field data to binary files..." << endl;
    cout << "Grid domain: [" << x_min << ", " << x_max << "] x [" 
         << y_min << ", " << y_max << "] x [" 
         << z_min << ", " << z_max << "]" << endl;
    cout << "Grid size: " << nx_global << " x " << ny_global << " x " << nz_global << " points" << endl;

    std::ofstream out_Ex("Ex.dat", std::ios::binary);
    std::ofstream out_Ey("Ey.dat", std::ios::binary);
    std::ofstream out_Ez("Ez.dat", std::ios::binary);
    std::ofstream out_Bx("Bx.dat", std::ios::binary);
    std::ofstream out_By("By.dat", std::ios::binary);
    std::ofstream out_Bz("Bz.dat", std::ios::binary);

    double dx = nx_global > 1 ? (x_max - x_min) / (nx_global - 1) : 0.0;
    double dy = ny_global > 1 ? (y_max - y_min) / (ny_global - 1) : 0.0;
    double dz = nz_global > 1 ? (z_max - z_min) / (nz_global - 1) : 0.0;

    std::vector<double> buf_Ex(nx_global);
    std::vector<double> buf_Ey(nx_global);
    std::vector<double> buf_Ez(nx_global);
    std::vector<double> buf_Bx(nx_global);
    std::vector<double> buf_By(nx_global);
    std::vector<double> buf_Bz(nx_global);

    for (int iz = 0; iz < nz_global; ++iz) {
        double z = z_min + iz * dz;
        for (int iy = 0; iy < ny_global; ++iy) {
            double y = y_min + iy * dy;
            for (int ix = 0; ix < nx_global; ++ix) {
                double x = x_min + ix * dx;
                
                double Ex, Ey, Ez, Bx, By, Bz;
                Gpulse(t, x, y, z, Ex, Ey, Ez, Bx, By, Bz);
                
                buf_Ex[ix] = Ex;
                buf_Ey[ix] = Ey;
                buf_Ez[ix] = Ez;
                buf_Bx[ix] = Bx;
                buf_By[ix] = By;
                buf_Bz[ix] = Bz;
            }
            out_Ex.write(reinterpret_cast<const char*>(buf_Ex.data()), nx_global * sizeof(double));
            out_Ey.write(reinterpret_cast<const char*>(buf_Ey.data()), nx_global * sizeof(double));
            out_Ez.write(reinterpret_cast<const char*>(buf_Ez.data()), nx_global * sizeof(double));
            out_Bx.write(reinterpret_cast<const char*>(buf_Bx.data()), nx_global * sizeof(double));
            out_By.write(reinterpret_cast<const char*>(buf_By.data()), nx_global * sizeof(double));
            out_Bz.write(reinterpret_cast<const char*>(buf_Bz.data()), nx_global * sizeof(double));
        }
    }
    cout << "Finished writing binary .dat files." << endl;
}

int main() {

    double t = 0.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double Ex, Ey, Ez, Bx, By, Bz;

    cout << "Testing Gpulse Library..." << endl;
    cout << "Polarization Mode: " << (RADIAL ? "Radial (Mosaic)" : "Linear") << endl;
    
    // Center point (should work without error)
    Gpulse(t, x, y, z, Ex, Ey, Ez, Bx, By, Bz);

    cout << "Field at center (t=0): Ez = " << Ez << endl;

    // Pick a point halfway between grid points (conceptually)
    // Grid range 4e-6, 40 points -> step ~ 0.2 um
    double x_off = 8.61e-7; // 0.861 um
    Gpulse(t, x_off, y, z, Ex, Ey, Ez, Bx, By, Bz);
    cout << "Field at x=0.861um (t=0): Ez = " << Ez << endl;
    
    // 使用预设的计算域和测试网格数来输出完整的场数据
    // 我们在这个测试程序中设定网格数为 num_x (100) 的一个合理子集/匹配值以便于测试
    int nx_global = num_x; // 从 constants.h 读取的参数（如 100）
    int ny_global = num_y; // 如 100
    int nz_global = 100;    // 
    
    dump_field_to_binary(t, 
                         nx_global, ny_global, nz_global,
                         -RANGE, RANGE,
                         -RANGE, RANGE,
                         -RANGE/2, RANGE/2);

    return 0;
}
