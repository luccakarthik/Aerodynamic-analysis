import tkinter as tk
from tkinter import messagebox, ttk
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class WeightEstimationApp:
    def __init__(self, root):
        self.root = root
        self.root.title("VTOL Aircraft Weight Estimation Tool")
        
        # Configure main frames
        self.left_frame = tk.Frame(root)
        self.left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.right_frame = tk.Frame(root)
        self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=False, padx=10, pady=10)
        
        # Initialize variables
        self.entry_vars = {
            'gross_weight': tk.StringVar(value='650'),
            'cruise_speed': tk.StringVar(value='40'),
            'sfc': tk.StringVar(value='0.3'),
            'ld': tk.StringVar(value='10.88'),
            'prop_efficiency': tk.StringVar(value='0.75'),
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

if __name__ == "__main__":
    root = tk.Tk()
    app = WeightEstimationApp(root)
    root.mainloop()
