import tkinter as tk
from tkinter import ttk
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class PowerCalculatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Aircraft Power Requirement Calculator")
        self.root.geometry("1000x800")
        
        # Variables with default values
        self.rho = tk.DoubleVar(value=0.589)     # Air density (kg/m³)
        self.velocity = tk.DoubleVar(value=40.0)  # Velocity (m/s)
        self.wing_area = tk.DoubleVar(value=10.0) # Wing area (m²)
        self.CL = tk.DoubleVar(value=1.3821)      # Lift coefficient (user input)
        self.CD0 = tk.DoubleVar(value=0.07866)    # Parasite drag coefficient
        self.e = tk.DoubleVar(value=0.8)          # Oswald efficiency
        self.AR = tk.DoubleVar(value=13.33)       # Aspect ratio
        self.eta = tk.DoubleVar(value=0.8)        # Propulsion efficiency
        self.use_CL_type = tk.StringVar(value="CL_3D")  # CL type selector
        
        # Create GUI
        self.create_inputs()
        self.create_outputs()
        self.create_plot_frame()
        self.create_info_panel()
        
    def create_inputs(self):
        input_frame = ttk.LabelFrame(self.root, text="Flight Parameters", padding=10)
        input_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=5)
        
        params = [
            ("Air Density (kg/m³):", self.rho),
            ("Velocity (m/s):", self.velocity),
            ("Wing Area (m²):", self.wing_area),
            ("Lift Coefficient (CL):", self.CL),
            ("Parasite Drag (CD0):", self.CD0),
            ("Oswald Efficiency:", self.e),
            ("Aspect Ratio:", self.AR),
            ("Propulsion Efficiency:", self.eta)
        ]
        
        for i, (text, var) in enumerate(params):
            ttk.Label(input_frame, text=text).grid(row=i, column=0, sticky=tk.W, pady=2)
            ttk.Entry(input_frame, textvariable=var, width=10).grid(row=i, column=1, pady=2)
        
        # CL type selector
        ttk.Label(input_frame, text="CL Type:").grid(row=len(params), column=0, sticky=tk.W, pady=2)
        cl_type_frame = ttk.Frame(input_frame)
        cl_type_frame.grid(row=len(params), column=1, sticky=tk.W)
        ttk.Radiobutton(cl_type_frame, text="CL_3D", variable=self.use_CL_type, value="CL_3D").pack(side=tk.LEFT)
        ttk.Radiobutton(cl_type_frame, text="CL_max_3D", variable=self.use_CL_type, value="CL_max_3D").pack(side=tk.LEFT)
        
        ttk.Button(input_frame, text="Calculate Power", 
                  command=self.calculate_power).grid(row=len(params)+1, column=0, columnspan=2, pady=10)
        
    def create_outputs(self):
        output_frame = ttk.LabelFrame(self.root, text="Calculation Results", padding=10)
        output_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)
        
        self.results_text = tk.Text(output_frame, height=10, wrap=tk.WORD)
        self.results_text.pack(fill=tk.BOTH, expand=True)
        self.results_text.insert(tk.END, "Enter parameters and click 'Calculate Power' to see results.")
        self.results_text.config(state=tk.DISABLED)
        
    def create_plot_frame(self):
        plot_frame = ttk.LabelFrame(self.root, text="Power vs Velocity", padding=10)
        plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.fig, self.ax = plt.subplots(figsize=(8, 4))
        self.ax.set_xlabel('Velocity (m/s)')
        self.ax.set_ylabel('Power Required (kW)')
        self.ax.grid(True)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
    def create_info_panel(self):
        info_frame = ttk.LabelFrame(self.root, text="CL Type Guidance", padding=10)
        info_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)
        
        info_text = """
CL Type Selection Guidance:
- CL_3D: Use this for normal cruise flight calculations (typical operational CL)
- CL_max_3D: Use this for maximum performance analysis (near stall conditions)

For most power calculations during cruise, use CL_3D (typically 70-90% of CL_max_3D)
"""
        ttk.Label(info_frame, text=info_text, justify=tk.LEFT).pack()
        
    def calculate_power(self):
        try:
            # Get input values
            rho = self.rho.get()
            V = self.velocity.get()
            S = self.wing_area.get()
            CL = self.CL.get()
            CD0 = self.CD0.get()
            e = self.e.get()
            AR = self.AR.get()
            eta = self.eta.get()
            cl_type = self.use_CL_type.get()
            
            # Calculate induced drag coefficient
            CDi = (CL**2) / (math.pi * e * AR)
            
            # Total drag coefficient
            CD = CD0 + CDi
            
            # Calculate drag force
            D = 0.5 * rho * V**2 * S * CD
            
            # Calculate power required
            P_cruise = (D * V) / eta
            
            # Generate velocity sweep for plot
            velocities = range(20, 81, 5)  # From 20 to 80 m/s in steps of 5
            powers = []
            for v in velocities:
                # For the plot, we'll assume CL changes with velocity to maintain lift
                # (L = 0.5*rho*V²*S*CL must equal weight)
                W = 0.5 * rho * V**2 * S * CL  # Reference weight at design point
                cl = (2 * W) / (rho * v**2 * S) if v > 5 else 0  # Avoid division by zero
                cdi = (cl**2) / (math.pi * e * AR) if v > 5 else 0
                cd = CD0 + cdi
                drag = 0.5 * rho * v**2 * S * cd
                power = (drag * v) / eta
                powers.append(power)
            
            # Update results text
            self.results_text.config(state=tk.NORMAL)
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, 
                                   f"=== Power Required Calculation ===\n\n"
                                   f"Using {cl_type}: {CL:.4f}\n"
                                   f"Induced Drag Coefficient (CDi): {CDi:.6f}\n"
                                   f"Total Drag Coefficient (CD): {CD:.6f}\n"
                                   f"Drag Force (D): {D:.2f} N\n"
                                   f"Power Required (P_cruise): {P_cruise/1000:.2f} kW\n\n"
                                   f"Assumptions:\n"
                                   f"- Propulsion efficiency (η): {eta:.2f}\n"
                                   f"- Oswald efficiency (e): {eta:.2f}\n"
                                   f"- Aspect ratio (AR): {AR:.2f}")
            self.results_text.config(state=tk.DISABLED)
            
            # Update plot
            self.ax.clear()
            self.ax.plot(velocities, [p/1000 for p in powers], 'b-', linewidth=2)
            self.ax.plot(V, P_cruise/1000, 'ro', markersize=8, label=f'Design Point ({V} m/s, {CL:.2f} CL)')
            self.ax.set_xlabel('Velocity (m/s)')
            self.ax.set_ylabel('Power Required (kW)')
            self.ax.set_title(f'Power Required vs Velocity (Using {cl_type})')
            self.ax.grid(True)
            self.ax.legend()
            self.canvas.draw()
            
        except Exception as e:
            self.results_text.config(state=tk.NORMAL)
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, f"Error in calculation: {str(e)}")
            self.results_text.config(state=tk.DISABLED)

if __name__ == "__main__":
    root = tk.Tk()
    app = PowerCalculatorApp(root)
    root.mainloop()
