import tkinter as tk
from tkinter import ttk
import math

class TailSizingApp:
    def __init__(self, root):
        self.root = root
        self.root.title("VTOL UAV Tail Sizing Calculator")
        self.root.geometry("800x600")
        
        # Variables
        self.wingspan = tk.DoubleVar(value=10.0)
        self.fuselage_ratio = tk.DoubleVar(value=0.55)
        self.wing_area = tk.DoubleVar(value=10.0)
        self.mac = tk.DoubleVar(value=1.0)
        self.h_tail_coeff = tk.DoubleVar(value=0.7)
        self.h_aspect_ratio = tk.DoubleVar(value=6.0)
        self.v_tail_coeff = tk.DoubleVar(value=0.04)
        self.v_aspect_ratio = tk.DoubleVar(value=2.0)
        
        # Create GUI
        self.create_inputs()
        self.create_outputs()
        
    def create_inputs(self):
        input_frame = ttk.LabelFrame(self.root, text="Input Parameters", padding=10)
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
        output_frame = ttk.LabelFrame(self.root, text="Tail Sizing Results", padding=10)
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

if __name__ == "__main__":
    root = tk.Tk()
    app = TailSizingApp(root)
    root.mainloop()
