import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
from matplotlib.widgets import Cursor

class AerodynamicAnalysisApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Advanced Aerodynamic Performance Analyzer")
        self.root.geometry("1400x900")
        
        # Input variables with expanded parameters
        self.inputs = {
            'weight': tk.DoubleVar(value=650),        # kg
            'density': tk.DoubleVar(value=0.589),     # kg/m³
            'wing_span': tk.DoubleVar(value=12),      # m
            'wing_chord': tk.DoubleVar(value=0.8),    # m
            'airspeed': tk.DoubleVar(value=40),       # m/s
            'cd0': tk.DoubleVar(value=0.07),
            'oswald': tk.DoubleVar(value=0.8),       # Oswald efficiency
            'cl_max_2d': tk.DoubleVar(value=1.8),    # Airfoil's max CL
            'sweep_angle': tk.DoubleVar(value=0),    # Wing sweep (deg)
            'mu': tk.DoubleVar(value=2.64e-5),       # Dynamic viscosity
            'sound_speed': tk.DoubleVar(value=343),  # Speed of sound (m/s)
            'alpha_range': tk.DoubleVar(value=15)     # Alpha range for plots
        }
        
        self.create_widgets()
        self.run_analysis()
    
    def create_widgets(self):
        # Configure grid weights for layout
        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_columnconfigure(1, weight=4)
        
        # Input frame (left panel - compact)
        input_frame = ttk.LabelFrame(self.root, text="Design Parameters", padding=5)
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
            ttk.Entry(input_frame, textvariable=self.inputs[param], width=10).grid(row=row, column=1, pady=1)
        
        # Calculate button with black color
        ttk.Button(input_frame, text="Calculate", command=self.run_analysis,
                  style='Black.TButton').grid(row=len(params), column=0, columnspan=2, pady=5)
        
        # Compact results display
        results_frame = ttk.LabelFrame(input_frame, text="Design Point Results", padding=3)
        results_frame.grid(row=len(params)+1, column=0, columnspan=2, sticky='ew', pady=3)
        
        self.results_text = tk.Text(results_frame, height=12, wrap=tk.WORD, font=('Courier', 8))
        self.results_text.pack(fill=tk.BOTH)
        
        # Formulas display
        formulas_frame = ttk.LabelFrame(input_frame, text="Formulas Used", padding=3)
        formulas_frame.grid(row=len(params)+2, column=0, columnspan=2, sticky='ew', pady=3)
        
        formulas_text = tk.Text(formulas_frame, height=8, wrap=tk.WORD, font=('Courier', 7))
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
        self.notebook = ttk.Notebook(self.root)
        self.notebook.grid(row=0, column=1, sticky="nsew", padx=5, pady=5)
        
        # Create frames for each tab
        self.wing_area_frame = ttk.Frame(self.notebook)
        self.alpha_frame = ttk.Frame(self.notebook)
        self.mach_frame = ttk.Frame(self.notebook)
        
        self.notebook.add(self.wing_area_frame, text="Analysis vs Wing Area")
        self.notebook.add(self.alpha_frame, text="Analysis vs Alpha")
        self.notebook.add(self.mach_frame, text="Analysis vs Mach")
        
        # Create figures for each tab
        self.fig_wing = plt.Figure(figsize=(9, 7), tight_layout=True)
        self.axs_wing = self.fig_wing.subplots(2, 2)
        
        self.fig_alpha = plt.Figure(figsize=(9, 7), tight_layout=True)
        self.axs_alpha = self.fig_alpha.subplots(2, 2)
        
        self.fig_mach = plt.Figure(figsize=(9, 7), tight_layout=True)
        self.axs_mach = self.fig_mach.subplots(2, 2)
        
        # Create canvases for each tab
        self.canvas_wing = FigureCanvasTkAgg(self.fig_wing, master=self.wing_area_frame)
        self.canvas_wing.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.canvas_alpha = FigureCanvasTkAgg(self.fig_alpha, master=self.alpha_frame)
        self.canvas_alpha.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.canvas_mach = FigureCanvasTkAgg(self.fig_mach, master=self.mach_frame)
        self.canvas_mach.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Add interactive cursors to all plots
        self.setup_cursors_and_hover()
        
        # Style configuration
        style = ttk.Style()
        style.configure('Black.TButton', foreground='white', background='black')
    
    def setup_cursors_and_hover(self):
        """Setup cursors and hover annotations for all plots"""
        self.cursors = []
        self.annot = []
        
        # Wing area plots
        for ax in self.axs_wing.flat:
            cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
            annot = ax.annotate("", xy=(0,0), xytext=(10,10), textcoords="offset points",
                             bbox=dict(boxstyle="round", fc="w"),
                             arrowprops=dict(arrowstyle="->"))
            annot.set_visible(False)
            self.cursors.append(cursor)
            self.annot.append(annot)
        
        # Alpha plots
        for ax in self.axs_alpha.flat:
            cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
            annot = ax.annotate("", xy=(0,0), xytext=(10,10), textcoords="offset points",
                             bbox=dict(boxstyle="round", fc="w"),
                             arrowprops=dict(arrowstyle="->"))
            annot.set_visible(False)
            self.cursors.append(cursor)
            self.annot.append(annot)
        
        # Mach plots
        for ax in self.axs_mach.flat:
            cursor = Cursor(ax, useblit=True, color='red', linewidth=1)
            annot = ax.annotate("", xy=(0,0), xytext=(10,10), textcoords="offset points",
                             bbox=dict(boxstyle="round", fc="w"),
                             arrowprops=dict(arrowstyle="->"))
            annot.set_visible(False)
            self.cursors.append(cursor)
            self.annot.append(annot)
        
        # Connect hover event to all canvases
        self.canvas_wing.mpl_connect("motion_notify_event", lambda e: self.hover(e, self.axs_wing, self.annot[:4]))
        self.canvas_alpha.mpl_connect("motion_notify_event", lambda e: self.hover(e, self.axs_alpha, self.annot[4:8]))
        self.canvas_mach.mpl_connect("motion_notify_event", lambda e: self.hover(e, self.axs_mach, self.annot[8:12]))
    
    def hover(self, event, axs, annots):
        """Handle mouse hover events for given subplots"""
        for i, ax in enumerate(axs.flat):
            if event.inaxes == ax:
                x, y = event.xdata, event.ydata
                annots[i].xy = (x, y)
                annots[i].set_text(f"x={x:.2f}\ny={y:.2f}")
                annots[i].set_visible(True)
                event.canvas.draw_idle()
            else:
                annots[i].set_visible(False)
                event.canvas.draw_idle()
    
    def calculate_cl_3d_max(self, cl_max_2d, AR, sweep_deg, Re):
        """Calculate realistic 3D max CL with sweep and Re effects"""
        # Base 3D correction
        cl_3d_max = 0.9 * cl_max_2d * np.cos(np.radians(sweep_deg))
        
        # Aspect ratio effect
        cl_3d_max *= (1 - 0.5/AR**0.7)
        
        # Reynolds number effect
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

if __name__ == "__main__":
    root = tk.Tk()
    app = AerodynamicAnalysisApp(root)
    root.mainloop()