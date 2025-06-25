import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class ConstraintAnalysisApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Advanced UAV Constraint Analysis")
        self.root.geometry("1400x900")
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)

        self.defaults = {
            'mtow': 650,          # kg
            'cruise_speed': 40,   # m/s
            'climb_speed': 35,    # m/s
            'climb_rate': 2.5,    # m/s
            'stall_speed': 30,    # m/s
            'cl_max_2d': 2.2,     # Airfoil's max CL (2D)
            'cd0': 0.07,
            'wing_span': 12,      # m
            'wing_chord': 0.8,    # m
            'e': 0.8,             # Oswald efficiency
            'bank_angle': 20,     # deg
            'ld_transition': 4,
            'altitude': 7000,     # m
            'sweep_angle': 0,     # deg
            'mu': 2.64e-5         # Dynamic viscosity
        }

        self.create_widgets()
        self.set_defaults()
        self.run_analysis()

    def on_close(self):
        plt.close('all')
        self.root.destroy()

    def create_widgets(self):
        # Input frame
        input_frame = ttk.LabelFrame(self.root, text="Aircraft Parameters", padding=10)
        input_frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10, pady=10)

        self.entries = {}
        params = [
            ('MTOW (kg)', 'mtow'),
            ('Cruise Speed (m/s)', 'cruise_speed'),
            ('Climb Speed (m/s)', 'climb_speed'),
            ('Climb Rate (m/s)', 'climb_rate'),
            ('Stall Speed (m/s)', 'stall_speed'),
            ('Airfoil CL_max (2D)', 'cl_max_2d'),
            ('CD0', 'cd0'),
            ('Wing Span (m)', 'wing_span'),
            ('Wing Chord (m)', 'wing_chord'),
            ('Oswald Efficiency', 'e'),
            ('Bank Angle (deg)', 'bank_angle'),
            ('L/D Transition', 'ld_transition'),
            ('Ceiling Altitude (m)', 'altitude'),
            ('Wing Sweep (deg)', 'sweep_angle'),
            ('Dynamic Viscosity', 'mu')
        ]

        for i, (label, key) in enumerate(params):
            ttk.Label(input_frame, text=label).grid(row=i, column=0, sticky='w', pady=2)
            self.entries[key] = ttk.Entry(input_frame)
            self.entries[key].grid(row=i, column=1, pady=2)

        ttk.Button(input_frame, text="Run Analysis", command=self.run_analysis).grid(
            row=len(params), column=0, columnspan=2, pady=10)

        # Constraint results table
        self.constraint_table = ttk.Treeview(input_frame, columns=('Constraint', 'T/W'), show='headings', height=7)
        self.constraint_table.heading('Constraint', text='Constraint')
        self.constraint_table.heading('T/W', text='T/W Required')
        self.constraint_table.column('Constraint', width=120)
        self.constraint_table.column('T/W', width=80)
        self.constraint_table.grid(row=len(params)+1, column=0, columnspan=2, pady=10, sticky='nsew')

        # Results and plot frame
        result_frame = ttk.Frame(self.root)
        result_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.fig, self.ax = plt.subplots(figsize=(8, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=result_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Expanded results display
        results_text = ttk.LabelFrame(result_frame, text="Detailed Analysis Results", padding=10)
        results_text.pack(fill=tk.BOTH, expand=True, pady=10)
        
        self.results_display = tk.Text(results_text, height=20, wrap=tk.WORD, font=('Courier', 10))
        self.results_display.pack(fill=tk.BOTH, expand=True)

    def set_defaults(self):
        for key, entry in self.entries.items():
            entry.insert(0, str(self.defaults[key]))

    def calculate_cl_3d_max(self, cl_max_2d, AR, sweep_deg, Re):
        """Calculate realistic 3D max CL with sweep and Re effects"""
        cl_3d_max = 0.9 * cl_max_2d * np.cos(np.radians(sweep_deg))
        cl_3d_max *= (1 - 0.5/AR**0.7)
        if Re < 1e6:
            cl_3d_max *= 0.8
        return cl_3d_max

    def get_density_at_altitude(self, altitude_m):
        """Standard atmosphere density calculation"""
        return 1.225 * np.exp(-altitude_m / 8500)

    def run_analysis(self):
        try:
            # Get all parameters
            params = {key: float(entry.get()) for key, entry in self.entries.items()}
            g = 9.81
            rho_sl = 1.225
            rho_alt = self.get_density_at_altitude(params['altitude'])
            
            # Calculate wing geometry
            S = params['wing_span'] * params['wing_chord']
            AR = params['wing_span']**2 / S
            c = S / params['wing_span']  # Mean chord
            
            # Calculate Reynolds number for CL_max correction
            Re_climb = (rho_sl * params['climb_speed'] * c) / params['mu']
            
            # Calculate 3D CL_max for stall condition
            cl_max_3d = self.calculate_cl_3d_max(
                params['cl_max_2d'], AR, params['sweep_angle'], Re_climb)
            
            # Wing loading range (N/m² and kg/m²)
            W_S_N = np.linspace(100, 1750, 250)
            W_S_kg = W_S_N / g
            
            # 1. Stall Speed Constraint (uses CL_max_3d)
            W_S_stall_N = 0.5 * rho_sl * params['stall_speed']**2 * cl_max_3d
            W_S_stall_kg = W_S_stall_N / g
            
            # 2. Cruise Constraint (uses operational CL)
            CL_cruise = W_S_N / (0.5 * rho_alt * params['cruise_speed']**2)
            CD_cruise = params['cd0'] + (CL_cruise**2) / (np.pi * AR * params['e'])
            T_W_cruise = CD_cruise / CL_cruise
            
            # 3. Climb Constraint
            CL_climb = W_S_N / (0.5 * rho_sl * params['climb_speed']**2)
            CD_climb = params['cd0'] + (CL_climb**2) / (np.pi * AR * params['e'])
            T_W_climb = (params['climb_rate'] / params['climb_speed']) + (CD_climb / CL_climb)
            
            # 4. Loiter (Endurance) Constraint
            CL_loiter = np.sqrt(params['cd0'] * np.pi * AR * params['e'])
            CD_loiter = params['cd0'] + (CL_loiter**2) / (np.pi * AR * params['e'])
            T_W_loiter = (CD_loiter / CL_loiter) * np.ones_like(W_S_N)
            
            # 5. Turn Constraint
            bank_rad = np.radians(params['bank_angle'])
            CL_turn = W_S_N / (0.5 * rho_alt * params['cruise_speed']**2)
            CD_turn = params['cd0'] + (CL_turn**2) / (np.pi * AR * params['e'])
            n = 1 / np.cos(bank_rad)
            T_W_turn = n * (CD_turn / CL_turn)
            
            # 6. Transition Constraint
            T_W_transition = 1 / params['ld_transition'] * np.ones_like(W_S_N)
            
            # 7. Ceiling Constraint
            R_c_ceiling = 0.5  # Excess power at ceiling
            CL_ceiling = W_S_N / (0.5 * rho_alt * params['cruise_speed']**2)
            CD_ceiling = params['cd0'] + (CL_ceiling**2) / (np.pi * AR * params['e'])
            T_W_ceiling = (R_c_ceiling / params['cruise_speed']) + (CD_ceiling / CL_ceiling)
            
            # Design point selection (90% of stall limit)
            W_S_design = 0.9 * W_S_stall_kg
            
            # Find required T/W at design point
            def interp(x, y, x0):
                return np.interp(x0, x, y)
            
            # Get T/W values at design point
            cruise_tw = interp(W_S_kg, T_W_cruise, W_S_design)
            climb_tw = interp(W_S_kg, T_W_climb, W_S_design)
            loiter_tw = interp(W_S_kg, T_W_loiter, W_S_design)
            turn_tw = interp(W_S_kg, T_W_turn, W_S_design)
            transition_tw = interp(W_S_kg, T_W_transition, W_S_design)
            ceiling_tw = interp(W_S_kg, T_W_ceiling, W_S_design)
            
            T_W_design = max([cruise_tw, climb_tw, loiter_tw, turn_tw, transition_tw, ceiling_tw])
            
            # Update constraint table
            self.constraint_table.delete(*self.constraint_table.get_children())
            constraints = [
                ('Stall Limit', f"{W_S_stall_kg:.2f} kg/m²"),
                ('Cruise', f"{cruise_tw:.4f}"),
                ('Climb', f"{climb_tw:.4f}"),
                ('Loiter', f"{loiter_tw:.4f}"),
                ('Turn', f"{turn_tw:.4f}"),
                ('Transition', f"{transition_tw:.4f}"),
                ('Ceiling', f"{ceiling_tw:.4f}"),
                ('Design Point', f"{T_W_design:.4f}")
            ]
            
            for constraint, value in constraints:
                self.constraint_table.insert('', 'end', values=(constraint, value))
            
            # Plotting
            self.ax.clear()
            self.ax.set_title("Aircraft Constraint Analysis")
            self.ax.set_xlabel("Wing Loading [kg/m²]")
            self.ax.set_ylabel("Thrust-to-Weight Ratio (T/W)")
            
            # Secondary x-axis for N/m²
            ax2 = self.ax.secondary_xaxis('top', functions=(lambda x: x * g, lambda x: x / g))
            ax2.set_xlabel("Wing Loading [N/m²]")
            
            # Plot all constraints
            self.ax.plot(W_S_kg, T_W_cruise, label='Cruise')
            self.ax.plot(W_S_kg, T_W_climb, label='Climb')
            self.ax.plot(W_S_kg, T_W_loiter, label='Loiter (Endurance)')
            self.ax.plot(W_S_kg, T_W_turn, label=f'Turn ({params["bank_angle"]}° bank)')
            self.ax.plot(W_S_kg, T_W_transition, label='Transition')
            self.ax.plot(W_S_kg, T_W_ceiling, label=f'Ceiling @ {params["altitude"]}m')
            self.ax.axvline(x=W_S_stall_kg, color='purple', linestyle='--', label='Stall Limit')
            
            # Mark design point
            self.ax.plot(W_S_design, T_W_design, 'ro', markersize=8)
            self.ax.annotate(
                f"Design Point\nW/S = {W_S_design:.1f} kg/m²\nT/W = {T_W_design:.3f}",
                (W_S_design, T_W_design), xytext=(10, 10), textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.9))
            
            self.ax.grid(True)
            self.ax.legend(loc='upper right')
            self.ax.set_ylim(0, 0.5)  # Adjusted for typical UAV values
            
            # Update results display with all calculated values
            self.results_display.delete(1.0, tk.END)
            results_text = f"""
=== WING GEOMETRY ===
Wing Span:       {params['wing_span']:.2f} m
Wing Chord:      {params['wing_chord']:.2f} m
Wing Area:       {S:.2f} m²
Aspect Ratio:    {AR:.2f}

=== LIFT COEFFICIENTS ===
2D CL_max:       {params['cl_max_2d']:.2f}
3D CL_max:       {cl_max_3d:.2f}
Operational CL (Cruise): {np.interp(W_S_design, W_S_kg, CL_cruise):.2f}
Operational CL (Climb):  {np.interp(W_S_design, W_S_kg, CL_climb):.2f}

=== WING LOADING ===
Stall Limit:     {W_S_stall_kg:.2f} kg/m² ({W_S_stall_N:.1f} N/m²)
Design Point:    {W_S_design:.2f} kg/m²

=== THRUST REQUIREMENTS ===
Cruise T/W:      {cruise_tw:.4f}
Climb T/W:       {climb_tw:.4f}
Loiter T/W:      {loiter_tw:.4f}
Turn T/W:        {turn_tw:.4f}
Transition T/W:  {transition_tw:.4f}
Ceiling T/W:     {ceiling_tw:.4f}

=== DESIGN POINT ===
Selected W/S:    {W_S_design:.2f} kg/m²
Required T/W:    {T_W_design:.4f}
"""
            self.results_display.insert(tk.END, results_text)
            
            self.canvas.draw()
            
        except ValueError as e:
            messagebox.showerror("Input Error", f"Invalid input value: {str(e)}")
        except Exception as e:
            messagebox.showerror("Error", f"Analysis failed: {str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = ConstraintAnalysisApp(root)
    root.mainloop()
