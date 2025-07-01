import tkinter as tk
from tkinter import ttk, messagebox
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import Cursor

class AerospaceDesignSuite:
    def __init__(self, root):
        self.root = root
        self.root.title("Aerospace Design Suite")
        self.root.geometry("1400x900")
        
        # Create notebook for tabs
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        # Create each tool as a tab
        self.create_aerodynamic_tab()
        self.create_constraint_tab()
        self.create_weight_tab()
        self.create_power_tab()
        self.create_tail_tab()
        
    def create_aerodynamic_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Aerodynamic Analysis")
        
        # Copy content from Aerodynamic_Analysis.py
        class AerodynamicAnalysisApp:
            def __init__(self, parent):
                self.parent = parent
                
                # Input variables with expanded parameters
                self.inputs = {
                    'weight': tk.DoubleVar(value=650),        # kg
                    'density': tk.DoubleVar(value=0.589),     # kg/m³
                    'wing_span': tk.DoubleVar(value=10),      # m
                    'wing_chord': tk.DoubleVar(value=0.75),    # m
                    'airspeed': tk.DoubleVar(value=40),       # m/s
                    'cd0': tk.DoubleVar(value=0.07),
                    'oswald': tk.DoubleVar(value=0.8),       # Oswald efficiency
                    'cl_max_2d': tk.DoubleVar(value=2.2),    # Airfoil's max CL
                    'sweep_angle': tk.DoubleVar(value=0),    # Wing sweep (deg)
                    'mu': tk.DoubleVar(value=2.64e-5),       # Dynamic viscosity
                    'sound_speed': tk.DoubleVar(value=303),  # Speed of sound (m/s)
                    'alpha_range': tk.DoubleVar(value=15)     # Alpha range for plots
                }
                
                self.create_widgets()
                self.run_analysis()
            
            def create_widgets(self):
                # Configure grid weights for layout
                self.parent.grid_columnconfigure(0, weight=1)
                self.parent.grid_columnconfigure(1, weight=4)
                
                # Input frame (left panel - compact)
                input_frame = ttk.LabelFrame(self.parent, text="Design Parameters", padding=5)
                input_frame.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
                
                # Create all input fields with reduced row height
                params = [
                    ('Weight (kg)', 'weight'),
                    ('Air Density (kg/m³)', 'density'),
                    ('Wing Span (m)', 'wing_span'),
                    ('Wing Chord (m)', 'wing_chord'),
                    ('Airspeed (m/s)', 'airspeed'),
                    ('CD0', 'cd0'),
                    ("Oswald Eff.", 'oswald'),
                    ('Airfoil CL_max', 'cl_max_2d'),
                    ('Wing Sweep (deg)', 'sweep_angle'),
                    ('Viscosity', 'mu'),
                    ('Speed of Sound (m/s)', 'sound_speed'),
                    ('Alpha Range (deg)', 'alpha_range')
                ]
                
                for row, (label, param) in enumerate(params):
                    ttk.Label(input_frame, text=label).grid(row=row, column=0, sticky='w', pady=1)
                    ttk.Entry(input_frame, textvariable=self.inputs[param], width=10).grid(row=row, column=0, pady=1)
                
                # Calculate button with black color
                ttk.Button(input_frame, text="Calculate", command=self.run_analysis,
                          style='White.TButton').grid(row=len(params), column=0, columnspan=2, pady=5)
                
                # Compact results display
                results_frame = ttk.LabelFrame(input_frame, text="Design Point Results", padding=3)
                results_frame.grid(row=len(params)+1, column=0, columnspan=1, sticky='ew', pady=3)
                
                self.results_text = tk.Text(results_frame, height=20, wrap=tk.WORD, font=('Courier', 10))
                self.results_text.pack(fill=tk.BOTH)
                
                # Formulas display
                formulas_frame = ttk.LabelFrame(input_frame, text="Formulas Used", padding=3)
                formulas_frame.grid(row=len(params)+2, column=0, columnspan=1, sticky='ew', pady=3)
                
                formulas_text = tk.Text(formulas_frame, height=12, wrap=tk.WORD, font=('Courier', 9))
                formulas_text.pack(fill=tk.BOTH)
                formulas = """
        Key Formulas:
        1. Wing Area (S) = Span × Chord
        2. Aspect Ratio (AR) = Span² / S
        3. CL_3D = (2 × Weight × g) / (ρ × V² × S)
        4. CD_total = CD0 + (CL²)/(π × AR × e)
        5. L/D = CL / CD_total
        6. Re = (ρ × V × c) / μ
        7. CL_max_3D = 0.9 × CL_max_2D × cos(sweep) × (1 - 0.5/AR^0.7)
        8. Mach = Airspeed / Speed of Sound
        """
                formulas_text.insert(tk.END, formulas)
                formulas_text.config(state='disabled')
                
                # Create notebook for different analysis tabs
                self.notebook = ttk.Notebook(self.parent)
                self.notebook.grid(row=0, column=1, sticky="nsew", padx=5, pady=5)
                
                # Create frames for each tab
                self.wing_area_frame = ttk.Frame(self.notebook)
                self.alpha_frame = ttk.Frame(self.notebook)
                self.mach_frame = ttk.Frame(self.notebook)
                
                self.notebook.add(self.wing_area_frame, text="Analysis vs Wing Area")
                self.notebook.add(self.alpha_frame, text="Analysis vs Alpha")
                self.notebook.add(self.mach_frame, text="Analysis vs Mach")
                
                # Create figures for each tab
                self.fig_wing = plt.Figure(figsize=(7, 5), tight_layout=True)
                self.axs_wing = self.fig_wing.subplots(2, 2)
                
                self.fig_alpha = plt.Figure(figsize=(7, 5), tight_layout=True)
                self.axs_alpha = self.fig_alpha.subplots(2, 2)
                
                self.fig_mach = plt.Figure(figsize=(7, 5), tight_layout=True)
                self.axs_mach = self.fig_mach.subplots(2, 2)
                
                # Create canvases for each tab
                self.canvas_wing = FigureCanvasTkAgg(self.fig_wing, master=self.wing_area_frame)
                self.canvas_wing.get_tk_widget().pack(fill=tk.BOTH, expand=True)
                
                self.canvas_alpha = FigureCanvasTkAgg(self.fig_alpha, master=self.alpha_frame)
                self.canvas_alpha.get_tk_widget().pack(fill=tk.BOTH, expand=True)
                
                self.canvas_mach = FigureCanvasTkAgg(self.fig_mach, master=self.mach_frame)
                self.canvas_mach.get_tk_widget().pack(fill=tk.BOTH, expand=True)
                
                # Style configuration
                style = ttk.Style()
                style.configure('Black.TButton', foreground='white', background='black')
            
            def calculate_cl_3d_max(self, cl_max_2d, AR, sweep_deg, Re):
                """Calculate realistic 3D max CL with sweep and Re effects"""
                cl_3d_max = 0.9 * cl_max_2d * np.cos(np.radians(sweep_deg))
                cl_3d_max *= (1 - 0.5/AR**0.7)
                if Re < 1e6:
                    cl_3d_max *= 0.8  # Reduction for low Reynolds numbers
                return cl_3d_max
            
            def calculate_aero(self, S, W, rho, V, b, CD0, e, cl_max_2d, sweep_deg, mu):
                """Calculate all aerodynamic parameters for given wing area"""
                # Calculate derived parameters
                AR = b**2 / S
                c = S / b  # Mean aerodynamic chord
                Re = (rho * V * c) / mu
                
                # Calculate maximum achievable CL (3D)
                cl_3d_max = self.calculate_cl_3d_max(cl_max_2d, AR, sweep_deg, Re)
                
                # Current flight CL (limited to 95% of max for safety)
                CL_required = (W * 9.81) / (0.5 * rho * V**2 * S)
                CL_3D = min(CL_required / (1 + (CL_required / (np.pi * AR * e))), cl_3d_max * 0.95)
                
                # Drag calculation
                CD_induced = (CL_3D**2) / (np.pi * AR * e)
                CD_total = CD0 + CD_induced
                LD = CL_3D / CD_total
                
                return {
                    'AR': AR,
                    'CL_3D': CL_3D,
                    'CL_3D_max': cl_3d_max,
                    'CD_total': CD_total,
                    'CD_induced': CD_induced,
                    'CD_parasitic': CD0,
                    'LD': LD,
                    'Re': Re,
                    'S': S
                }
            
            def calculate_aero_alpha(self, alpha, S, W, rho, V, b, CD0, e, cl_max_2d, sweep_deg, mu):
                """Calculate aerodynamic parameters at specific angle of attack"""
                AR = b**2 / S
                c = S / b
                Re = (rho * V * c) / mu
                
                # Simplified CL vs alpha relationship (linear until stall)
                alpha_0 = -4  # Zero-lift angle of attack (degrees)
                cl_alpha = 0.1  # Lift curve slope per degree
                
                # Calculate CL for this alpha (limited by max CL)
                CL_alpha = min(cl_alpha * (alpha - alpha_0), cl_max_2d * 0.95)
                
                # Drag calculation
                CD_induced = (CL_alpha**2) / (np.pi * AR * e)
                CD_total = CD0 + CD_induced
                LD = CL_alpha / CD_total
                
                return {
                    'AR': AR,
                    'CL_3D': CL_alpha,
                    'CD_total': CD_total,
                    'CD_induced': CD_induced,
                    'CD_parasitic': CD0,
                    'LD': LD,
                    'alpha': alpha
                }
            
            def calculate_aero_mach(self, mach, S, W, rho, V, b, CD0, e, cl_max_2d, sweep_deg, mu, sound_speed):
                """Calculate aerodynamic parameters at specific Mach number"""
                # Convert Mach to airspeed
                V_mach = mach * sound_speed
                
                # Recalculate parameters with new airspeed
                AR = b**2 / S
                c = S / b
                Re = (rho * V_mach * c) / mu
                
                # Calculate CL required for lift
                CL_required = (W * 9.81) / (0.5 * rho * V_mach**2 * S)
                
                # Calculate maximum achievable CL (3D)
                cl_3d_max = self.calculate_cl_3d_max(cl_max_2d, AR, sweep_deg, Re)
                CL_3D = min(CL_required / (1 + (CL_required / (np.pi * AR * e))), cl_3d_max * 0.95)
                
                # Drag calculation with Mach effects
                CD_induced = (CL_3D**2) / (np.pi * AR * e)
                
                # Simple Mach drag rise approximation
                if mach > 0.6:
                    CD0_mach = CD0 * (1 + 0.2 * (mach - 0.6)**2)
                else:
                    CD0_mach = CD0
                    
                CD_total = CD0_mach + CD_induced
                LD = CL_3D / CD_total
                
                return {
                    'AR': AR,
                    'CL_3D': CL_3D,
                    'CD_total': CD_total,
                    'CD_induced': CD_induced,
                    'CD_parasitic': CD0_mach,
                    'LD': LD,
                    'mach': mach
                }
            
            def run_analysis(self):
                try:
                    # Get all input values
                    W = self.inputs['weight'].get()
                    rho = self.inputs['density'].get()
                    b = self.inputs['wing_span'].get()
                    c = self.inputs['wing_chord'].get()
                    V = self.inputs['airspeed'].get()
                    CD0 = self.inputs['cd0'].get()
                    e = self.inputs['oswald'].get()
                    cl_max_2d = self.inputs['cl_max_2d'].get()
                    sweep_deg = self.inputs['sweep_angle'].get()
                    mu = self.inputs['mu'].get()
                    sound_speed = self.inputs['sound_speed'].get()
                    alpha_range = self.inputs['alpha_range'].get()
                    
                    # Calculate nominal wing area
                    S_nominal = b * c
                    
                    # ===== Wing Area Analysis =====
                    S_min = max(0.1, 0.7 * S_nominal)
                    S_max = 1.3 * S_nominal
                    S_values = np.linspace(S_min, S_max, 25)
                    
                    wing_results = [self.calculate_aero(S, W, rho, V, b, CD0, e, cl_max_2d, sweep_deg, mu) 
                                  for S in S_values]
                    
                    nominal_result = next(r for r in wing_results if abs(r['S'] - S_nominal) < 1e-6)
                    self.plot_wing_area_results(S_values, wing_results, S_nominal)
                    
                    # ===== Alpha Analysis =====
                    alpha_min = -5
                    alpha_max = alpha_range
                    alpha_values = np.linspace(alpha_min, alpha_max, 25)
                    
                    alpha_results = [self.calculate_aero_alpha(alpha, S_nominal, W, rho, V, b, CD0, e, 
                                                             cl_max_2d, sweep_deg, mu) 
                                   for alpha in alpha_values]
                    
                    self.plot_alpha_results(alpha_values, alpha_results)
                    
                    # ===== Mach Analysis =====
                    current_mach = V / sound_speed
                    mach_min = max(0.1, current_mach - 0.2)  # 4 steps before (assuming 0.05 step)
                    mach_max = current_mach + 0.2            # 4 steps after
                    mach_values = np.linspace(mach_min, mach_max, 25)
                    
                    mach_results = [self.calculate_aero_mach(mach, S_nominal, W, rho, V, b, CD0, e,
                                                           cl_max_2d, sweep_deg, mu, sound_speed)
                                  for mach in mach_values]
                    
                    self.plot_mach_results(mach_values, mach_results, current_mach)
                    
                    # Update results text
                    self.update_results_text(nominal_result)
                    
                    # Draw all canvases
                    self.canvas_wing.draw()
                    self.canvas_alpha.draw()
                    self.canvas_mach.draw()
                    
                except Exception as e:
                    tk.messagebox.showerror("Error", f"Calculation failed:\n{str(e)}")
            
            def plot_wing_area_results(self, S_values, results, S_nominal):
                """Create wing area performance plots"""
                for ax in self.axs_wing.flat:
                    ax.clear()
                
                AR = [r['AR'] for r in results]
                CL = [r['CL_3D'] for r in results]
                CL_max = [r['CL_3D_max'] for r in results]
                CD_total = [r['CD_total'] for r in results]
                CD_induced = [r['CD_induced'] for r in results]
                CD_parasitic = [r['CD_parasitic'] for r in results]
                LD = [r['LD'] for r in results]
                
                legend_props = {'loc': 'upper right', 'fontsize': 8, 'framealpha': 0.7}
                
                # Plot 1: Lift Coefficients
                self.axs_wing[0,0].plot(S_values, CL, 'b-', label='CL (operational)')
                self.axs_wing[0,0].plot(S_values, CL_max, 'b--', label='CL_max (3D)')
                self.axs_wing[0,0].axvline(S_nominal, color='k', linestyle='--', label='Design Point')
                self.axs_wing[0,0].set_xlabel('Wing Area (m²)')
                self.axs_wing[0,0].set_ylabel('Lift Coefficient')
                self.axs_wing[0,0].set_title('Lift Coefficient vs Wing Area')
                self.axs_wing[0,0].legend(**legend_props)
                self.axs_wing[0,0].grid(True)
                
                # Plot 2: Drag Breakdown
                self.axs_wing[0,1].plot(S_values, CD_total, 'r-', label='Total CD')
                self.axs_wing[0,1].plot(S_values, CD_induced, 'g--', label='Induced CD')
                self.axs_wing[0,1].plot(S_values, CD_parasitic, 'b--', label='Parasitic CD')
                self.axs_wing[0,1].axvline(S_nominal, color='k', linestyle='--')
                self.axs_wing[0,1].set_xlabel('Wing Area (m²)')
                self.axs_wing[0,1].set_ylabel('Drag Coefficient')
                self.axs_wing[0,1].set_title('Drag Breakdown')
                self.axs_wing[0,1].legend(**legend_props)
                self.axs_wing[0,1].grid(True)
                
                # Plot 3: L/D Ratio
                self.axs_wing[1,0].plot(S_values, LD, 'm-')
                self.axs_wing[1,0].axvline(S_nominal, color='k', linestyle='--')
                self.axs_wing[1,0].set_xlabel('Wing Area (m²)')
                self.axs_wing[1,0].set_ylabel('L/D Ratio')
                self.axs_wing[1,0].set_title('Aerodynamic Efficiency')
                self.axs_wing[1,0].grid(True)
                
                # Plot 4: Aspect Ratio
                self.axs_wing[1,1].plot(S_values, AR, 'c-')
                self.axs_wing[1,1].axvline(S_nominal, color='k', linestyle='--')
                self.axs_wing[1,1].set_xlabel('Wing Area (m²)')
                self.axs_wing[1,1].set_ylabel('Aspect Ratio')
                self.axs_wing[1,1].set_title('Aspect Ratio Variation')
                self.axs_wing[1,1].grid(True)
                
                self.fig_wing.tight_layout()
            
            def plot_alpha_results(self, alpha_values, results):
                """Create angle of attack performance plots"""
                for ax in self.axs_alpha.flat:
                    ax.clear()
                
                AR = [r['AR'] for r in results]
                CL = [r['CL_3D'] for r in results]
                CD_total = [r['CD_total'] for r in results]
                CD_induced = [r['CD_induced'] for r in results]
                CD_parasitic = [r['CD_parasitic'] for r in results]
                LD = [r['LD'] for r in results]
                
                legend_props = {'loc': 'upper right', 'fontsize': 8, 'framealpha': 0.7}
                
                # Plot 1: Lift Coefficient vs Alpha
                self.axs_alpha[0,0].plot(alpha_values, CL, 'b-')
                self.axs_alpha[0,0].set_xlabel('Angle of Attack (deg)')
                self.axs_alpha[0,0].set_ylabel('Lift Coefficient')
                self.axs_alpha[0,0].set_title('Lift Coefficient vs Alpha')
                self.axs_alpha[0,0].grid(True)
                
                # Plot 2: Drag Breakdown vs Alpha
                self.axs_alpha[0,1].plot(alpha_values, CD_total, 'r-', label='Total CD')
                self.axs_alpha[0,1].plot(alpha_values, CD_induced, 'g--', label='Induced CD')
                self.axs_alpha[0,1].plot(alpha_values, CD_parasitic, 'b--', label='Parasitic CD')
                self.axs_alpha[0,1].set_xlabel('Angle of Attack (deg)')
                self.axs_alpha[0,1].set_ylabel('Drag Coefficient')
                self.axs_alpha[0,1].set_title('Drag Breakdown vs Alpha')
                self.axs_alpha[0,1].legend(**legend_props)
                self.axs_alpha[0,1].grid(True)
                
                # Plot 3: L/D Ratio vs Alpha
                self.axs_alpha[1,0].plot(alpha_values, LD, 'm-')
                self.axs_alpha[1,0].set_xlabel('Angle of Attack (deg)')
                self.axs_alpha[1,0].set_ylabel('L/D Ratio')
                self.axs_alpha[1,0].set_title('Aerodynamic Efficiency vs Alpha')
                self.axs_alpha[1,0].grid(True)
                
                # Plot 4: Aspect Ratio vs Alpha (constant in this case)
                self.axs_alpha[1,1].plot(alpha_values, AR, 'c-')
                self.axs_alpha[1,1].set_xlabel('Angle of Attack (deg)')
                self.axs_alpha[1,1].set_ylabel('Aspect Ratio')
                self.axs_alpha[1,1].set_title('Aspect Ratio (Constant)')
                self.axs_alpha[1,1].grid(True)
                
                self.fig_alpha.tight_layout()
            
            def plot_mach_results(self, mach_values, results, current_mach):
                """Create Mach number performance plots"""
                for ax in self.axs_mach.flat:
                    ax.clear()
                
                AR = [r['AR'] for r in results]
                CL = [r['CL_3D'] for r in results]
                CD_total = [r['CD_total'] for r in results]
                CD_induced = [r['CD_induced'] for r in results]
                CD_parasitic = [r['CD_parasitic'] for r in results]
                LD = [r['LD'] for r in results]
                
                legend_props = {'loc': 'upper right', 'fontsize': 8, 'framealpha': 0.7}
                
                # Plot 1: Lift Coefficient vs Mach
                self.axs_mach[0,0].plot(mach_values, CL, 'b-')
                self.axs_mach[0,0].axvline(current_mach, color='k', linestyle='--', label='Current Mach')
                self.axs_mach[0,0].set_xlabel('Mach Number')
                self.axs_mach[0,0].set_ylabel('Lift Coefficient')
                self.axs_mach[0,0].set_title('Lift Coefficient vs Mach')
                self.axs_mach[0,0].legend(**legend_props)
                self.axs_mach[0,0].grid(True)
                
                # Plot 2: Drag Breakdown vs Mach
                self.axs_mach[0,1].plot(mach_values, CD_total, 'r-', label='Total CD')
                self.axs_mach[0,1].plot(mach_values, CD_induced, 'g--', label='Induced CD')
                self.axs_mach[0,1].plot(mach_values, CD_parasitic, 'b--', label='Parasitic CD')
                self.axs_mach[0,1].axvline(current_mach, color='k', linestyle='--')
                self.axs_mach[0,1].set_xlabel('Mach Number')
                self.axs_mach[0,1].set_ylabel('Drag Coefficient')
                self.axs_mach[0,1].set_title('Drag Breakdown vs Mach')
                self.axs_mach[0,1].legend(**legend_props)
                self.axs_mach[0,1].grid(True)
                
                # Plot 3: L/D Ratio vs Mach
                self.axs_mach[1,0].plot(mach_values, LD, 'm-')
                self.axs_mach[1,0].axvline(current_mach, color='k', linestyle='--')
                self.axs_mach[1,0].set_xlabel('Mach Number')
                self.axs_mach[1,0].set_ylabel('L/D Ratio')
                self.axs_mach[1,0].set_title('Aerodynamic Efficiency vs Mach')
                self.axs_mach[1,0].grid(True)
                
                # Plot 4: Aspect Ratio vs Mach (constant in this case)
                self.axs_mach[1,1].plot(mach_values, AR, 'c-')
                self.axs_mach[1,1].axvline(current_mach, color='k', linestyle='--')
                self.axs_mach[1,1].set_xlabel('Mach Number')
                self.axs_mach[1,1].set_ylabel('Aspect Ratio')
                self.axs_mach[1,1].set_title('Aspect Ratio (Constant)')
                self.axs_mach[1,1].grid(True)
                
                self.fig_mach.tight_layout()
            
            def update_results_text(self, design_point):
                """Update the results text box with design point information"""
                text = f"""DESIGN POINT @ {design_point['S']:.2f} m²:

        Wing Geometry:
          Span = {self.inputs['wing_span'].get():.2f} m
          Chord = {self.inputs['wing_chord'].get():.2f} m
          AR = {design_point['AR']:.2f}

        Aerodynamics:
          CL = {design_point['CL_3D']:.4f}
          CL_max = {design_point['CL_3D_max']:.4f}
          Re = {design_point['Re']:.1e}

        Drag:
          Total = {design_point['CD_total']:.4f}
            Parasitic = {design_point['CD_parasitic']:.4f}
            Induced = {design_point['CD_induced']:.4f}

        Efficiency:
          L/D = {design_point['LD']:.2f}"""

                self.results_text.delete(1.0, tk.END)
                self.results_text.insert(tk.END, text)

        # Initialize the aerodynamic analysis app in this tab
        AerodynamicAnalysisApp(tab)

    def create_constraint_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Constraint Analysis")
        
        # Copy content from Latest_Constraint_analysis_with_results.py
        class ConstraintAnalysisApp:
            def __init__(self, parent):
                self.parent = parent
                self.defaults = {
                    'mtow': 650,          # kg
                    'cruise_speed': 40,   # m/s
                    'climb_speed': 35,    # m/s
                    'climb_rate': 2.5,    # m/s
                    'stall_speed': 30,    # m/s
                    'cl_max_2d': 2.2,     # Airfoil's max CL (2D)
                    'cd0': 0.07,
                    'wing_span': 10,      # m
                    'wing_chord': 0.75,    # m
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
            
            def create_widgets(self):
                # Input frame
                input_frame = ttk.LabelFrame(self.parent, text="Aircraft Parameters", padding=10)
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
                result_frame = ttk.Frame(self.parent)
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

        # Initialize the constraint analysis app in this tab
        ConstraintAnalysisApp(tab)

    def create_weight_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Weight Estimation")
        
        # Copy content from Latest_Final_Weight_Fraction.py
        class WeightEstimationApp:
            def __init__(self, parent):
                self.parent = parent
                
                # Configure main frames
                self.left_frame = tk.Frame(parent)
                self.left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)
                
                self.right_frame = tk.Frame(parent)
                self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=False, padx=10, pady=10)
                
                # Initialize variables
                self.entry_vars = {
                    'gross_weight': tk.StringVar(value='650'),
                    'cruise_speed': tk.StringVar(value='40'),
                    'sfc': tk.StringVar(value='0.3'),
                    'ld': tk.StringVar(value='10.88'),
                    'prop_efficiency': tk.StringVar(value='0.7'),
                    'cruise_range': tk.StringVar(value='200'),
                    'loiter_time': tk.StringVar(value='5'),
                    'climb_factor': tk.StringVar(value='0.04')  # 4% climb fuel factor
                }
                
                self.setup_input_section()
                self.setup_results_section()
                self.setup_formulas_section()
                self.setup_graph_section()

            def setup_input_section(self):
                input_frame = tk.LabelFrame(self.left_frame, text="Input Parameters")
                input_frame.pack(fill=tk.X, padx=5, pady=5)
                
                fields = [
                    ("Gross Weight (kg)", 'gross_weight', "Total weight of the aircraft including payload, fuel, and structure."),
                    ("Cruise Speed (m/s)", 'cruise_speed', "Constant velocity at which the aircraft cruises during mission."),
                    ("SFC (kg/kWh)", 'sfc', "Specific Fuel Consumption η between 0.1 and 0.9 (for propeller aircraft)."),
                    ("L/D Ratio", 'ld', "Lift-to-Drag ratio for both cruise and loiter phases."),
                    ("Propeller Efficiency", 'prop_efficiency', "Propeller efficiency η between 0.1 and 0.9."),
                    ("Cruise Range (km)", 'cruise_range', "Total one-way cruise range in kilometers."),
                    ("Loiter Time (hr)", 'loiter_time', "Loiter duration in hours."),
                    ("Climb Fuel Factor", 'climb_factor', "Fuel fraction used during climb (typically 3-5%).")
                ]

                for idx, (label_text, var_key, tooltip_text) in enumerate(fields):
                    lbl = tk.Label(input_frame, text=label_text)
                    lbl.grid(row=idx, column=0, padx=5, pady=3, sticky='e')
                    entry = tk.Entry(input_frame, textvariable=self.entry_vars[var_key], width=15)
                    entry.grid(row=idx, column=1, padx=5, pady=3)
                    self.create_tooltip(entry, tooltip_text)

                btn = tk.Button(input_frame, text="Run Estimation", command=self.run_estimation)
                btn.grid(row=len(fields), column=0, columnspan=2, pady=10)

            def setup_results_section(self):
                results_frame = tk.LabelFrame(self.left_frame, text="Results")
                results_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
                
                columns = ("Component", "Weight (kg)", "Fraction")
                self.tree = ttk.Treeview(results_frame, columns=columns, show='headings', height=12)
                
                for col in columns:
                    self.tree.heading(col, text=col)
                    self.tree.column(col, anchor='center', width=100)
                    
                self.tree.pack(fill=tk.BOTH, expand=True)

            def setup_formulas_section(self):
                formulas_frame = tk.LabelFrame(self.left_frame, text="Formulas Used")
                formulas_frame.pack(fill=tk.BOTH, padx=5, pady=5)
                
                formulas = (
                    "Raymer's Empty Weight Formula:\n"
                    "  Wₑ/W₀ = A × W₀^C × Kᵥₛ\n"
                    "  (A = 1.53, C = -0.16, Kᵥₛ = 1.0)\n\n"
                    "Breguet Range Equation:\n"
                    "  W₂/W₁ = exp[–(R × SFC) / (η × V × (L/D))]\n\n"
                    "Breguet Endurance Equation:\n"
                    "  W₃/W₂ = exp[–(E × SFC) / (η × (L/D))]\n\n"
                    "Mission Segments:\n"
                    "  W1/W0 = 1.0 (Start)\n"
                    "  W2/W1 = 0.96 (Climb, 4% fuel)\n"
                    "  W3/W2 = Outbound cruise\n"
                    "  W4/W3 = Loiter\n"
                    "  W5/W4 = Return cruise\n\n"
                    "Final Payload Equation:\n"
                    "  Wₚₐyₗₒₐd = W₀ × (1 – Wₑ/W₀ – Wf/W₀)"
                )
                
                formula_label = tk.Label(formulas_frame, text=formulas, justify='left', font=('Courier', 9), anchor='nw')
                formula_label.pack(fill=tk.BOTH, expand=True)

            def setup_graph_section(self):
                # Top graph frame (Range sensitivity)
                self.top_graph_frame = tk.LabelFrame(self.right_frame, text="Range Sensitivity", height=300)
                self.top_graph_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
                
                # Bottom graph frame (Endurance sensitivity)
                self.bottom_graph_frame = tk.LabelFrame(self.right_frame, text="Endurance Sensitivity", height=300)
                self.bottom_graph_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

            def create_tooltip(self, widget, text):
                tip = tk.Toplevel(widget)
                tip.withdraw()
                tip.overrideredirect(True)
                lbl = tk.Label(tip, text=text, background='#ffffe0', relief='solid', borderwidth=1, justify='left')
                lbl.pack(ipadx=1)

                def enter(event):
                    tip.deiconify()
                    x = widget.winfo_rootx() + widget.winfo_width() + 5
                    y = widget.winfo_rooty()
                    tip.geometry(f"+{x}+{y}")

                def leave(event):
                    tip.withdraw()

                widget.bind('<Enter>', enter)
                widget.bind('<Leave>', leave)

            def calculate_fuel_fraction(self, params):
                w1_w0 = 1.0
                w2_w1 = 1 - float(self.entry_vars['climb_factor'].get())  # Climb fuel factor (4%)
                
                # Outbound cruise
                range_m = params['cruise_range'] * 1000
                exponent = (range_m * params['sfc']) / (params['prop_efficiency'] * params['cruise_speed'] * params['ld'] * 3600)
                w3_w2 = math.exp(-exponent)
                
                # Loiter
                exponent = (params['loiter_time'] * params['sfc']) / (params['prop_efficiency'] * params['ld'])
                w4_w3 = math.exp(-exponent)
                
                # Return cruise (same as outbound)
                w5_w4 = w3_w2
                
                w5_w0 = w5_w4 * w4_w3 * w3_w2 * w2_w1 * w1_w0
                
                return {
                    'total_fuel_fraction': 1 - w5_w0,
                    'mission_segments': {
                        'w1_w0': w1_w0,
                        'w2_w1': w2_w1,
                        'w3_w2': w3_w2,
                        'w4_w3': w4_w3,
                        'w5_w4': w5_w4
                    }
                }

            def estimate_weights(self, params):
                empty_weight_fraction = calculate_empty_weight_fraction(params['gross_weight'])
                fuel_result = self.calculate_fuel_fraction(params)
                fuel_fraction = fuel_result['total_fuel_fraction']
                denominator = 1 - fuel_fraction - empty_weight_fraction
                
                if denominator <= 0:
                    raise ValueError("Cannot achieve mission with given parameters - negative payload")
                    
                payload_weight = params['gross_weight'] * denominator
                return {
                    'gross_weight': params['gross_weight'],
                    'empty_weight': params['gross_weight'] * empty_weight_fraction,
                    'fuel_weight': params['gross_weight'] * fuel_fraction,
                    'payload_weight': payload_weight,
                    'empty_weight_fraction': empty_weight_fraction,
                    'fuel_fraction': fuel_fraction,
                    'payload_fraction': payload_weight / params['gross_weight'],
                    'mission_segments': fuel_result['mission_segments']
                }

            def create_sensitivity_plot(self, frame, x_values, y1_values, y2_values, x_label, current_value):
                for widget in frame.winfo_children():
                    widget.destroy()
                    
                fig, ax = plt.subplots(figsize=(5, 3))
                ax.plot(x_values, y1_values, 'r-', label='Fuel Weight')
                ax.plot(x_values, y2_values, 'b-', label='Payload Weight')
                ax.axvline(x=current_value, color='k', linestyle='--', label='Current Value')
                ax.set_xlabel(x_label)
                ax.set_ylabel('Weight (kg)')
                ax.legend()
                ax.grid(True)
                
                canvas = FigureCanvasTkAgg(fig, master=frame)
                canvas.draw()
                canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

            def run_estimation(self):
                try:
                    params = {
                        'gross_weight': float(self.entry_vars['gross_weight'].get()),
                        'cruise_speed': float(self.entry_vars['cruise_speed'].get()),
                        'sfc': float(self.entry_vars['sfc'].get()),
                        'ld': float(self.entry_vars['ld'].get()),
                        'prop_efficiency': float(self.entry_vars['prop_efficiency'].get()),
                        'cruise_range': float(self.entry_vars['cruise_range'].get()),
                        'loiter_time': float(self.entry_vars['loiter_time'].get())
                    }

                    if not (0.1 <= params['sfc'] <= 0.9):
                        raise ValueError("SFC must be between 0.1 and 0.9")
                    if not (0.1 <= params['prop_efficiency'] <= 0.9):
                        raise ValueError("Propeller Efficiency must be between 0.1 and 0.9")

                    results = self.estimate_weights(params)

                    # Clear previous results
                    for row in self.tree.get_children():
                        self.tree.delete(row)

                    # Main weight components
                    self.tree.insert('', 'end', values=("Empty Weight (We)", f"{results['empty_weight']:.2f}", f"{results['empty_weight_fraction']:.4f}"))
                    self.tree.insert('', 'end', values=("Payload Weight", f"{results['payload_weight']:.2f}", f"{results['payload_fraction']:.4f}"))
                    self.tree.insert('', 'end', values=("Fuel Weight (Wfuel)", f"{results['fuel_weight']:.2f}", f"{results['fuel_fraction']:.4f}"))
                    self.tree.insert('', 'end', values=("Total Weight (W0)", f"{results['gross_weight']:.2f}", "1.0000"))
                    
                    # Mission segment fuel fractions
                    self.tree.insert('', 'end', values=("", "", ""))
                    self.tree.insert('', 'end', values=("Mission Segment Fuel Fractions", "", ""))
                    self.tree.insert('', 'end', values=("Start (W1/W0)", "", f"{results['mission_segments']['w1_w0']:.4f}"))
                    self.tree.insert('', 'end', values=("Climb (W2/W1)", "", f"{results['mission_segments']['w2_w1']:.4f}"))
                    self.tree.insert('', 'end', values=("Outbound Cruise (W3/W2)", "", f"{results['mission_segments']['w3_w2']:.4f}"))
                    self.tree.insert('', 'end', values=("Loiter (W4/W3)", "", f"{results['mission_segments']['w4_w3']:.4f}"))
                    self.tree.insert('', 'end', values=("Return Cruise (W5/W4)", "", f"{results['mission_segments']['w5_w4']:.4f}"))

                    # Create sensitivity plots
                    original_range = params['cruise_range']
                    range_values = [original_range * (1 + 0.2*i) for i in [-1.5, -1, -0.5, 0, 0.5, 1, 1.5]]
                    fuel_weights = []
                    payload_weights = []
                    
                    for r in range_values:
                        temp_params = params.copy()
                        temp_params['cruise_range'] = r
                        try:
                            temp_results = self.estimate_weights(temp_params)
                            fuel_weights.append(temp_results['fuel_weight'])
                            payload_weights.append(temp_results['payload_weight'])
                        except:
                            fuel_weights.append(float('nan'))
                            payload_weights.append(float('nan'))
                    
                    self.create_sensitivity_plot(
                        self.top_graph_frame,
                        range_values,
                        fuel_weights,
                        payload_weights,
                        'Cruise Range (km)',
                        original_range
                    )
                    
                    original_endurance = params['loiter_time']
                    endurance_values = [original_endurance * (1 + 0.2*i) for i in [-1.5, -1, -0.5, 0, 0.5, 1, 1.5]]
                    fuel_weights = []
                    payload_weights = []
                    
                    for e in endurance_values:
                        temp_params = params.copy()
                        temp_params['loiter_time'] = e
                        try:
                            temp_results = self.estimate_weights(temp_params)
                            fuel_weights.append(temp_results['fuel_weight'])
                            payload_weights.append(temp_results['payload_weight'])
                        except:
                            fuel_weights.append(float('nan'))
                            payload_weights.append(float('nan'))
                    
                    self.create_sensitivity_plot(
                        self.bottom_graph_frame,
                        endurance_values,
                        fuel_weights,
                        payload_weights,
                        'Loiter Time (hr)',
                        original_endurance
                    )

                except ValueError as e:
                    messagebox.showerror("Input Error", str(e))
                except Exception as e:
                    messagebox.showerror("Unexpected Error", str(e))

        def calculate_empty_weight_fraction(gross_weight):
            A = 1.53
            C = -0.16
            K_vs = 1.0
            return A * (gross_weight ** C) * K_vs

        # Initialize the weight estimation app in this tab
        WeightEstimationApp(tab)

    def create_power_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Power Analysis")
        
        # Copy content from Power_Anaysis.py
        class PowerCalculatorApp:
            def __init__(self, parent):
                self.parent = parent
                
                # Variables with your specified values
                self.rho = tk.DoubleVar(value=0.589)       # kg/m³ at 7000m
                self.velocity = tk.DoubleVar(value=40.0)   # m/s
                self.wing_area = tk.DoubleVar(value=7.5)   # m²
                self.weight = tk.DoubleVar(value=650.0)    # kg
                self.CL = tk.DoubleVar(value=1.3821)       # Fixed CL value
                self.CD0 = tk.DoubleVar(value=0.07866)     # Parasite drag
                self.e = tk.DoubleVar(value=0.8)           # Oswald efficiency
                self.AR = tk.DoubleVar(value=13.33)        # Aspect ratio
                self.eta = tk.DoubleVar(value=0.7)        # Propulsion efficiency
                self.climb_sin_theta = tk.DoubleVar(value=0.0625)  # sin(3.58°)
                self.k_vtol = tk.DoubleVar(value=1.5)      # VTOL factor
                self.safety_margin = tk.DoubleVar(value=15.0)  # %
                
                self.create_inputs()
                self.create_outputs()
                self.create_plot_frame()
            
            def create_inputs(self):
                input_frame = ttk.LabelFrame(self.parent, text="Flight Parameters", padding=10)
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
                output_frame = ttk.LabelFrame(self.parent, text="Calculation Results", padding=10)
                output_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)
                
                self.results_text = tk.Text(output_frame, height=18, wrap=tk.WORD)
                self.results_text.pack(fill=tk.BOTH, expand=True)
                self.results_text.insert(tk.END, "Parameters will be calculated here...")
                self.results_text.config(state=tk.DISABLED)
            
            def create_plot_frame(self):
                plot_frame = ttk.LabelFrame(self.parent, text="Power vs Velocity", padding=10)
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

        # Initialize the power calculator app in this tab
        PowerCalculatorApp(tab)

    def create_tail_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Tail Sizing")
        
        # Copy content from Tail_Sizing.py
        class TailSizingApp:
            def __init__(self, parent):
                self.parent = parent
                
                # Variables
                self.wingspan = tk.DoubleVar(value=10.0)
                self.fuselage_ratio = tk.DoubleVar(value=0.55)
                self.wing_area = tk.DoubleVar(value=7.5)
                self.mac = tk.DoubleVar(value=0.75)
                self.h_tail_coeff = tk.DoubleVar(value=1.1)
                self.h_aspect_ratio = tk.DoubleVar(value=7)
                self.v_tail_coeff = tk.DoubleVar(value=0.04)
                self.v_aspect_ratio = tk.DoubleVar(value=2.0)
                
                # Create GUI
                self.create_inputs()
                self.create_outputs()
            
            def create_inputs(self):
                input_frame = ttk.LabelFrame(self.parent, text="Input Parameters", padding=10)
                input_frame.pack(fill=tk.X, padx=10, pady=5)
                
                # Wing parameters
                ttk.Label(input_frame, text="Wingspan (m):").grid(row=0, column=0, sticky=tk.W)
                ttk.Entry(input_frame, textvariable=self.wingspan).grid(row=0, column=1)
                
                ttk.Label(input_frame, text="Fuselage Length Ratio:").grid(row=1, column=0, sticky=tk.W)
                ttk.Entry(input_frame, textvariable=self.fuselage_ratio).grid(row=1, column=1)
                
                ttk.Label(input_frame, text="Wing Area (m²):").grid(row=2, column=0, sticky=tk.W)
                ttk.Entry(input_frame, textvariable=self.wing_area).grid(row=2, column=1)
                
                ttk.Label(input_frame, text="Mean Aerodynamic Chord (m):").grid(row=3, column=0, sticky=tk.W)
                ttk.Entry(input_frame, textvariable=self.mac).grid(row=3, column=1)
                
                # Horizontal tail parameters
                ttk.Label(input_frame, text="Horizontal Tail Volume Coefficient:").grid(row=4, column=0, sticky=tk.W)
                ttk.Entry(input_frame, textvariable=self.h_tail_coeff).grid(row=4, column=1)
                
                ttk.Label(input_frame, text="Horizontal Tail Aspect Ratio:").grid(row=5, column=0, sticky=tk.W)
                ttk.Entry(input_frame, textvariable=self.h_aspect_ratio).grid(row=5, column=1)
                
                # Vertical tail parameters
                ttk.Label(input_frame, text="Vertical Tail Volume Coefficient:").grid(row=6, column=0, sticky=tk.W)
                ttk.Entry(input_frame, textvariable=self.v_tail_coeff).grid(row=6, column=1)
                
                ttk.Label(input_frame, text="Vertical Tail Aspect Ratio:").grid(row=7, column=0, sticky=tk.W)
                ttk.Entry(input_frame, textvariable=self.v_aspect_ratio).grid(row=7, column=1)
                
                # Calculate button
                ttk.Button(input_frame, text="Calculate", command=self.calculate).grid(row=8, column=0, columnspan=2, pady=10)
            
            def create_outputs(self):
                output_frame = ttk.LabelFrame(self.parent, text="Tail Sizing Results", padding=10)
                output_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
                
                # Horizontal tail results
                ttk.Label(output_frame, text="Horizontal Tail:").grid(row=0, column=0, sticky=tk.W)
                self.h_results = ttk.Label(output_frame, text="", wraplength=700)
                self.h_results.grid(row=1, column=0, sticky=tk.W)
                
                # Vertical tail results
                ttk.Label(output_frame, text="Vertical Tail:").grid(row=2, column=0, sticky=tk.W, pady=(10,0))
                self.v_results = ttk.Label(output_frame, text="", wraplength=700)
                self.v_results.grid(row=3, column=0, sticky=tk.W)
                
                # General info
                ttk.Label(output_frame, text="Note: For twin boom configuration, vertical tail area is split between two surfaces.").grid(row=4, column=0, sticky=tk.W, pady=(10,0))
            
            def calculate(self):
                try:
                    # Calculate fuselage length
                    fuselage_length = self.fuselage_ratio.get() * self.wingspan.get()
                    
                    # Calculate tail arm (70% of fuselage length)
                    tail_arm = 0.7 * fuselage_length
                    
                    # Horizontal tail calculations
                    h_area = (self.h_tail_coeff.get() * self.wing_area.get() * self.mac.get()) / tail_arm
                    h_span = math.sqrt(h_area * self.h_aspect_ratio.get())
                    h_chord = h_area / h_span
                    
                    # Vertical tail calculations (total area)
                    v_total_area = (self.v_tail_coeff.get() * self.wingspan.get() * self.wing_area.get()) / tail_arm
                    v_area_per_boom = v_total_area / 2
                    v_span = math.sqrt(v_area_per_boom * self.v_aspect_ratio.get())
                    v_chord = v_area_per_boom / v_span
                    
                    # Update results
                    h_text = f"Area: {h_area:.3f} m² | Span: {h_span:.3f} m | Chord: {h_chord:.3f} m | Aspect Ratio: {self.h_aspect_ratio.get():.1f}"
                    v_text = f"Total Area: {v_total_area:.3f} m² | Per Boom Area: {v_area_per_boom:.3f} m² | Span: {v_span:.3f} m | Chord: {v_chord:.3f} m | Aspect Ratio: {self.v_aspect_ratio.get():.1f}"
                    
                    self.h_results.config(text=h_text)
                    self.v_results.config(text=v_text)
                    
                except Exception as e:
                    self.h_results.config(text=f"Error: {str(e)}")
                    self.v_results.config(text="")

        # Initialize the tail sizing app in this tab
        TailSizingApp(tab)

if __name__ == "__main__":
    root = tk.Tk()
    app = AerospaceDesignSuite(root)
    root.mainloop()
