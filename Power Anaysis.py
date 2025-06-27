import tkinter as tk
from tkinter import ttk
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class PowerCalculatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Aircraft Power Calculator - Consistent CL")
        self.root.geometry("1300x1000")
        
        # Variables with your specified values
        self.rho = tk.DoubleVar(value=0.589)       # kg/m³ at 7000m
        self.velocity = tk.DoubleVar(value=40.0)   # m/s
        self.wing_area = tk.DoubleVar(value=7.5)   # m²
        self.weight = tk.DoubleVar(value=650.0)    # kg
        self.CL = tk.DoubleVar(value=1.3821)       # Fixed CL value
        self.CD0 = tk.DoubleVar(value=0.07866)     # Parasite drag
        self.e = tk.DoubleVar(value=0.8)           # Oswald efficiency
        self.AR = tk.DoubleVar(value=13.33)        # Aspect ratio
        self.eta = tk.DoubleVar(value=0.75)        # Propulsion efficiency
        self.climb_sin_theta = tk.DoubleVar(value=0.0625)  # sin(3.58°)
        self.k_vtol = tk.DoubleVar(value=1.5)      # VTOL factor
        self.safety_margin = tk.DoubleVar(value=15.0)  # %
        
        self.create_inputs()
        self.create_outputs()
        self.create_plot_frame()
        
    def create_inputs(self):
        input_frame = ttk.LabelFrame(self.root, text="Flight Parameters", padding=10)
        input_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=5)
        
        params = [
            ("Air Density (kg/m³):", self.rho),
            ("Cruise Velocity (m/s):", self.velocity),
            ("Wing Area (m²):", self.wing_area),
            ("Weight (kg):", self.weight),
            ("Lift Coefficient (CL):", self.CL),
            ("Parasite Drag (CD0):", self.CD0),
            ("Oswald Efficiency:", self.e),
            ("Aspect Ratio:", self.AR),
            ("Propulsion Efficiency:", self.eta),
            ("sinθ (from climb):", self.climb_sin_theta),
            ("VTOL Drag Factor (k):", self.k_vtol),
            ("Safety Margin (%):", self.safety_margin)
        ]
        
        for i, (text, var) in enumerate(params):
            ttk.Label(input_frame, text=text).grid(row=i, column=0, sticky=tk.W, pady=2)
            ttk.Entry(input_frame, textvariable=var, width=10).grid(row=i, column=1, pady=2)
        
        ttk.Button(input_frame, text="Calculate Power", 
                 command=self.calculate_power).grid(row=len(params), column=0, columnspan=2, pady=10)
        
    def create_outputs(self):
        output_frame = ttk.LabelFrame(self.root, text="Calculation Results", padding=10)
        output_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)
        
        self.results_text = tk.Text(output_frame, height=18, wrap=tk.WORD)
        self.results_text.pack(fill=tk.BOTH, expand=True)
        self.results_text.insert(tk.END, "Parameters will be calculated here...")
        self.results_text.config(state=tk.DISABLED)
        
    def create_plot_frame(self):
        plot_frame = ttk.LabelFrame(self.root, text="Power vs Velocity", padding=10)
        plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.fig, self.ax = plt.subplots(figsize=(10, 6))
        self.ax.set_xlabel('Velocity (m/s)')
        self.ax.set_ylabel('Power Required (kW)')
        self.ax.grid(True, linestyle='--', alpha=0.7)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
    def calculate_power(self):
        try:
            # Get input values
            rho = self.rho.get()
            V = self.velocity.get()
            S = self.wing_area.get()
            W = self.weight.get() * 9.81
            CL = self.CL.get()  # Using fixed CL throughout
            CD0 = self.CD0.get()
            e = self.e.get()
            AR = self.AR.get()
            eta = self.eta.get()
            sin_theta = self.climb_sin_theta.get()
            k = self.k_vtol.get()
            margin = self.safety_margin.get() / 100
            
            # 1. Cruise Power Calculation (using fixed CL)
            CDi = (CL**2) / (math.pi * e * AR)
            CD = CD0 + CDi
            D_cruise = 0.5 * rho * V**2 * S * CD
            P_cruise = (D_cruise * V) / eta
            
            # 2. Climb Power Calculation (using fixed CL)
            D_climb = D_cruise * (1 + k * sin_theta)
            weight_component = W * sin_theta
            T_climb = D_climb + weight_component
            P_climb = (T_climb * V) / eta
            
            # 3. Peak Power with Safety Margin
            P_peak = P_climb * (1 + margin)
            
            # Generate velocity array that includes exact V
            velocities = [x for x in range(20, 81, 1)]  # 1 m/s steps
            if V not in velocities:
                velocities.append(V)
                velocities.sort()
            
            # Calculate power curves using FIXED CL for all velocities
            cruise_powers = []
            climb_powers = []
            
            for v in velocities:
                # Use the same fixed CL for all calculations
                cd = CD0 + (CL**2) / (math.pi * e * AR)
                drag = 0.5 * rho * v**2 * S * cd
                cruise_powers.append((drag * v) / eta)
                
                # Climb power with fixed CL
                drag_climb = drag * (1 + k * sin_theta)
                thrust_climb = drag_climb + (W * sin_theta)
                climb_powers.append((thrust_climb * v) / eta)
            
            # Get exact index for design velocity
            v_idx = velocities.index(V)
            
            # Update results
            self.results_text.config(state=tk.NORMAL)
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END,
                f"=== CRUISE POWER @ {V} m/s ===\n"
                f"Using fixed CL: {CL:.4f}\n"
                f"CDi = {CDi:.6f}\n"
                f"CD = {CD:.6f}\n"
                f"Drag = {D_cruise:.2f} N\n"
                f"Power = {P_cruise/1000:.3f} kW\n\n"
                
                f"=== CLIMB POWER @ {V} m/s ===\n"
                f"Drag (climb) = {D_climb:.2f} N\n"
                f"Weight component = {weight_component:.2f} N\n"
                f"Thrust required = {T_climb:.2f} N\n"
                f"Power = {P_climb/1000:.3f} kW\n\n"
                
                f"=== PEAK POWER ===\n"
                f"With {self.safety_margin.get():.0f}% margin: {P_peak/1000:.3f} kW\n\n"
                
                f"Note: Using fixed CL={CL:.4f} for all calculations"
            )
            self.results_text.config(state=tk.DISABLED)
            
            # Update plot
            self.ax.clear()
            
            # Plot curves
            self.ax.plot(velocities, [p/1000 for p in cruise_powers], 'b-', linewidth=2, label=f'Cruise Power (CL={CL:.2f})')
            self.ax.plot(velocities, [p/1000 for p in climb_powers], 'g-', linewidth=2, label='Climb Power')
            self.ax.axhline(y=P_peak/1000, color='r', linestyle='--', linewidth=2, 
                           label=f'Peak (+{self.safety_margin.get():.0f}%)')
            
            # Mark design points
            self.ax.plot(V, cruise_powers[v_idx]/1000, 'bo', markersize=8)
            self.ax.plot(V, climb_powers[v_idx]/1000, 'go', markersize=8)
            
            # Annotate with exact values
            self.ax.annotate(f'{cruise_powers[v_idx]/1000:.3f} kW', 
                           (V, cruise_powers[v_idx]/1000), 
                           xytext=(10, -10), textcoords='offset points',
                           bbox=dict(boxstyle="round", fc="white", alpha=0.8))
            self.ax.annotate(f'{climb_powers[v_idx]/1000:.3f} kW', 
                           (V, climb_powers[v_idx]/1000), 
                           xytext=(-60, 7), textcoords='offset points',
                           bbox=dict(boxstyle="round", fc="white", alpha=0.8))
            self.ax.annotate(f'{P_peak/1000:.3f} kW', 
                           (V, P_peak/1000), 
                           xytext=(10, 10), textcoords='offset points',
                           bbox=dict(boxstyle="round", fc="white", alpha=0.8))
            
            # Formatting
            self.ax.set_xlabel('Velocity (m/s)')
            self.ax.set_ylabel('Power Required (kW)')
            self.ax.set_title('Aircraft Power Requirements (Fixed CL)')
            self.ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
            self.fig.tight_layout()
            self.canvas.draw()
            
        except Exception as e:
            self.results_text.config(state=tk.NORMAL)
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, f"Calculation Error: {str(e)}")
            self.results_text.config(state=tk.DISABLED)

if __name__ == "__main__":
    root = tk.Tk()
    app = PowerCalculatorApp(root)
    root.mainloop()
