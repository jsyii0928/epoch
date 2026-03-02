"""
激光脉冲时域功率可视化工具
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt

def plot_temporal_power():
    filename = 'temporal_intensity.h5'
    
    try:
        with h5py.File(filename, 'r') as f:
            time = np.array(f['time']) * 1e15  # Convert to fs
            # Intensity is stored in W/m^2 (SI).
            # Convert to W/cm^2: 1 W/m^2 = 1e-4 W/cm^2
            intensity_si = np.array(f['intensity'])
            intensity = intensity_si * 1e-4 
            
            # Find peak
            peak_idx = np.argmax(intensity)
            peak_val = intensity[peak_idx]
            peak_time = time[peak_idx]
            
            # --- FWHM Calculation ---
            half_max = peak_val / 2.0
            
            # Find left crossing (before peak)
            if peak_idx > 0:
                left_idx = np.argmin(np.abs(intensity[:peak_idx] - half_max))
                # Linear interpolation for better precision
                t1 = time[left_idx]
                p1 = intensity[left_idx]
                # Use neighbor to interpolate
                if p1 < half_max: 
                    t2, p2 = time[left_idx+1], intensity[left_idx+1]
                else:
                    t2, p2 = time[left_idx-1], intensity[left_idx-1]
                t_left = t1 + (half_max - p1) * (t2 - t1) / (p2 - p1)
            else:
                t_left = time[0]

            # Find right crossing (after peak)
            if peak_idx < len(intensity) - 1:
                right_idx = np.argmin(np.abs(intensity[peak_idx:] - half_max)) + peak_idx
                t1 = time[right_idx]
                p1 = intensity[right_idx]
                if p1 < half_max:
                    t2, p2 = time[right_idx-1], intensity[right_idx-1]
                else:
                    t2, p2 = time[right_idx+1], intensity[right_idx+1]
                t_right = t1 + (half_max - p1) * (t2 - t1) / (p2 - p1)
            else:
                t_right = time[-1]
            
            fwhm = t_right - t_left
            print(f"Calculated FWHM: {fwhm:.2f} fs")

            # Plot
            plt.figure(figsize=(10, 6))
            plt.plot(time, intensity, 'b-', linewidth=2, label=f'Intensity (Peak={peak_val:.2e} W/cm$^2$)')
            plt.fill_between(time, intensity, alpha=0.3, color='blue')
            
            # Plot FWHM
            plt.hlines(half_max, t_left, t_right, colors='r', linestyles='dashed', linewidth=2)
            plt.text(peak_time, half_max * 1.05, f'FWHM = {fwhm:.1f} fs', 
                     ha='center', va='bottom', color='r', fontweight='bold', fontsize=12)
            plt.plot([t_left, t_right], [half_max, half_max], 'ro') # Mark points
            
            # --- Simulation Duration Recommendation ---
            # Criterion: Time span where Power > 1% of Peak (-30 dB intensity contrast)
            threshold = peak_val * 1e-2
            valid_indices = np.where(intensity > threshold)[0]
            if len(valid_indices) > 0:
                t_start_rec = time[valid_indices[0]]
                t_end_rec = time[valid_indices[-1]]
                duration_rec = t_end_rec - t_start_rec
                
                print("\n" + "="*40)
                print(f"    Simulation Duration Recommendation")
                print("="*40)
                print(f"Based on 1% Peak Intensity Threshold (-20dB):")
                print(f"  - Active window: [{t_start_rec:.1f}, {t_end_rec:.1f}] fs")
                print(f"  - Minimum Duration: {duration_rec:.1f} fs")
                print(f"  - Suggested Simulation Time: > {duration_rec * 1.2:.0f} fs (with margin)")
                print("="*40 + "\n")
            
            # 读取偏振属性 (如果存在)
            try:
                radial = f.attrs['radial']
                pol_str = "Radial" if radial else "Linear"
            except KeyError:
                pol_str = "Unknown"
            
            plt.title(f'On-Axis Intensity Envelope vs Time ({pol_str} Pol.)', fontsize=14)
            plt.xlabel('Time (fs)', fontsize=12)
            plt.ylabel('Intensity (W/cm$^2$)', fontsize=12)
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.legend()
            
            output_file = f'temporal_intensity_{pol_str}.png'
            plt.savefig(output_file, dpi=300)
            print(f"图表已保存到 {output_file}")
            
    except FileNotFoundError:
        print(f"找不到文件 {filename}。请先运行 ./export_temporal")
    except Exception as e:
        print(f"发生错误: {e}")

if __name__ == "__main__":
    plot_temporal_power()
