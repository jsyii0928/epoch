"""
激光场分布可视化工具
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt

def plot_fields():
    """
    读取并绘制激光场分布 (6个子图: Ex, Ey, Ez, Bx, By, Bz)
    """
    filename = 'field_snapshot.h5'
    
    # Configuration
    PLOT_ABSOLUTE = 0  # Set to True for Magnitude, False for Real part

    try:
        with h5py.File(filename, 'r') as f:
            # 读取并重构复数电场
            Ex = np.array(f['Ex_r']) + 1j * np.array(f['Ex_i'])
            Ey = np.array(f['Ey_r']) + 1j * np.array(f['Ey_i'])
            Ez = np.array(f['Ez_r']) + 1j * np.array(f['Ez_i'])
            
            # 读取并重构复数磁场
            Bx = np.array(f['Bx_r']) + 1j * np.array(f['Bx_i'])
            By = np.array(f['By_r']) + 1j * np.array(f['By_i'])
            Bz = np.array(f['Bz_r']) + 1j * np.array(f['Bz_i'])
            
            range_val = f.attrs['range']
            
            # 读取时间属性
            try: t_val = f.attrs['t']
            except KeyError: t_val = 0.0
                
            # 读取偏振属性
            try:
                radial = f.attrs['radial']
                pol_str = "Radial" if radial else "Linear"
            except KeyError:
                pol_str = "Unknown"

            N = Ex.shape[0]
            
            # Generate grids (Meshgrid)
            x_vec = np.linspace(-range_val, range_val, N)
            y_vec = np.linspace(-range_val, range_val, N)
            X, Y = np.meshgrid(x_vec, y_vec, indexing='ij')
            
            # Helper for processing
            def process_field(field_c):
                return np.abs(field_c) if PLOT_ABSOLUTE else np.real(field_c)

            # 投影逻辑
            if radial:
                # 计算方位角
                Theta = np.arctan2(Y, X)
                CosT = np.cos(Theta)
                SinT = np.sin(Theta)
                
                # 投影复数场到极坐标基
                # Er = Ex cos + Ey sin
                # Eth = -Ex sin + Ey cos
                Er_c = Ex * CosT + Ey * SinT
                Eth_c = -Ex * SinT + Ey * CosT
                
                Bth_c = -Bx * SinT + By * CosT
                Br_c = Bx * CosT + By * SinT # Oops previously I might have missed defining Br_c correctly if it wasn't there, checking...
                # Wait, original code had:
                # Br_c = Bx * CosT + By * SinT
                # Bth_c = -Bx * SinT + By * CosT

                E1_data, E1_name = process_field(Er_c), 'Er'
                E2_data, E2_name = process_field(Eth_c), 'Etheta'
                B1_data, B1_name = process_field(Br_c), 'Br'
                B2_data, B2_name = process_field(Bth_c), 'Btheta'
            else:
                E1_data, E1_name = process_field(Ex), 'Ex'
                E2_data, E2_name = process_field(Ey), 'Ey'
                B1_data, B1_name = process_field(Bx), 'Bx'
                B2_data, B2_name = process_field(By), 'By'
                
            # Ez Bz are scalars
            Ez_data = process_field(Ez)
            Bz_data = process_field(Bz)

            t_fs = t_val * 1e15
            mode_str = "Abs" if PLOT_ABSOLUTE else "Real"

            fig, ax = plt.subplots(2, 3, figsize=(15, 9))
            fig.suptitle(f'Laser Field ({mode_str}, {pol_str}) Distribution @ t = {t_fs:.1f} fs', fontsize=16)
            
            def plot_subplot(axes_obj, data, title, unit):
                # Pass data as is because X,Y match indexing='ij'
                im = axes_obj.pcolormesh(X*1e6, Y*1e6, data, shading='auto', cmap='turbo')
                axes_obj.set_title(f'{title} ({unit})')
                plt.colorbar(im, ax=axes_obj)
                axes_obj.set_xlabel('x (um)')
                axes_obj.set_ylabel('y (um)')
                axes_obj.axis('scaled')

            plot_subplot(ax[0,0], E1_data, E1_name, 'V/m')
            plot_subplot(ax[0,1], E2_data, E2_name, 'V/m')
            plot_subplot(ax[0,2], Ez_data, 'Ez', 'V/m')

            plot_subplot(ax[1,0], B1_data, B1_name, 'T')
            plot_subplot(ax[1,1], B2_data, B2_name, 'T')
            plot_subplot(ax[1,2], Bz_data, 'Bz', 'T')

            plt.tight_layout()
            
            output_filename = f'field_plot_{t_fs:.1f}fs_{pol_str}_{mode_str}.png'
            plt.savefig(output_filename, dpi=300)
            print(f"图形已保存到 {output_filename}")

    except FileNotFoundError:
        print(f"错误: 找不到文件 {filename}。请先运行 C++ 导出程序生成数据文件。")
    except Exception as e:
        print(f"发生错误: {e}")

if __name__ == "__main__":
    plot_fields()
