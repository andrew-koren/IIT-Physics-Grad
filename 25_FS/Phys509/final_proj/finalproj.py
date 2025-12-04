import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interactive, FloatSlider, IntSlider, VBox, Layout
from IPython.display import display

def run_rabi_viz():
    """
    Creates and displays an interactive widget for Rabi Oscillations.
    Returns the widget object.
    """
    
    # --- 1. SETUP THE FIGURE ---
    # We use plt.ioff() to prevent the figure from showing up twice 
    # (once when created, once when the widget displays)
    plt.ioff()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4.5))
    fig.canvas.header_visible = False  # Hide the "Figure 1" header to save space
    plt.tight_layout()

    # --- Initialize Empty Plots: Subplot 1 (Energy) ---
    ln_dressed_plus, = ax1.plot([], [], 'b-', linewidth=2, label=r'Dressed $|1\rangle$')
    ln_dressed_minus, = ax1.plot([], [], 'r-', linewidth=2, label=r'Dressed $|1\rangle$')
    ln_bare_1, = ax1.plot([], [], 'k--', alpha=0.3, label='Bare States')
    ln_bare_2, = ax1.plot([], [], 'k--', alpha=0.3)
    ln_current_dot = ax1.scatter([], [], color='green', s=100, zorder=5)
    ln_current_line = ax1.axvline(x=0, color='green', linestyle=':', alpha=0.8, label=r'Current Delta')

    ax1.set_xlabel(r'Detuning (Delta)')
    ax1.set_ylabel('Energy')
    ax1.set_title('Energy Spectrum')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right', fontsize=8)
    ax1.set_xlim(-5, 5)
    ax1.set_ylim(-6, 6)

    # --- Initialize Empty Plots: Subplot 2 (Time) ---
    ln_prob_e, = ax2.plot([], [], 'purple', label=r'$P_e (|e, n\rangle)$', linewidth=2)
    ln_prob_g, = ax2.plot([], [], 'orange', label=r'$P_g (|g, n\rangle)$', linewidth=2, linestyle='--')
    txt_contrast = ax2.text(0.05, 0.5, "", transform=ax2.transAxes, bbox=dict(facecolor='white', alpha=0.8))

    ax2.set_xlabel('Time (t)')
    ax2.set_ylabel('Probability')
    ax2.set_title('Rabi Oscillations')
    ax2.set_ylim(-0.05, 1.05)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right')

    # Turn interactive mode back on so updates work
    plt.ion()

    # --- 2. UPDATE LOGIC ---
    def update_plot(detuning, g, n, time_range):
        # Math
        delta_scan = np.linspace(-10, 10, 200)
        coupling_term = 2 * g * np.sqrt(n)
        rabi_freq = np.sqrt(detuning**2 + coupling_term**2)
        
        # Energy Data
        E_plus = 0.5 * np.sqrt(delta_scan**2 + coupling_term**2)
        E_minus = -0.5 * np.sqrt(delta_scan**2 + coupling_term**2)
        current_E_plus = 0.5 * rabi_freq
        current_E_minus = -0.5 * rabi_freq
        
        # Update Energy Plots
        ln_dressed_plus.set_data(delta_scan, E_plus)
        ln_dressed_minus.set_data(delta_scan, E_minus)
        ln_bare_1.set_data(delta_scan, 0.5 * delta_scan)
        ln_bare_2.set_data(delta_scan, -0.5 * delta_scan)
        
        ln_current_dot.set_offsets(np.c_[[detuning, detuning], [current_E_plus, current_E_minus]])
        ln_current_line.set_xdata([detuning, detuning])
        
        # Time Data
        t = np.linspace(0, time_range, 500)
        if rabi_freq == 0:
            amp = 0
        else:
            amp = (coupling_term / rabi_freq)**2
        
        P_g = amp * np.sin(rabi_freq * t / 2)**2
        P_e = 1 - P_g
        
        # Update Time Plots
        ln_prob_e.set_data(t, P_e)
        ln_prob_g.set_data(t, P_g)
        ax2.set_xlim(0, time_range)
        
        txt_contrast.set_text(f"Contrast: {amp:.2f}\n$\Omega_R$: {rabi_freq:.2f}")
        
        fig.canvas.draw_idle()

    # --- 3. CREATE WIDGET ---
    widget = interactive(update_plot, 
             detuning=FloatSlider(value=0.5, min=-5.0, max=5.0, step=0.1, description='Detuning'),
             g=FloatSlider(value=1.0, min=0.0, max=2.0, step=0.1, description='Coupling g'),
             n=IntSlider(value=1, min=1, max=10, step=1, description='Photons n'),
             time_range=FloatSlider(value=10.0, min=1.0, max=50.0, step=1.0, description='Max Time')
            )
    
    # Arrange: Figure on top, sliders on bottom
    # widget.children[-1] is the output (the figure, partially) but matplotlib widget handles fig separately
    # so we construct a clean VBox with the figure object and the controls
    
    controls = VBox(widget.children[:-1])
    layout = VBox([fig.canvas, controls])
    
    return layout

import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interactive, FloatSlider, IntSlider, VBox, HTML
from IPython.display import display

def default_ramp(t, alpha):
    """Default linear ramp: Omega = alpha * t"""
    return alpha * t

def run_adiabatic_viz(omega_func=default_ramp, alpha_max=2.0):
    """
    Creates an interactive widget for Adiabatic Evolution.
    
    Parameters:
    omega_func: A python function f(t, alpha) returning the coupling strength.
    alpha_max: The maximum range for the alpha slider.
    """
    
    # --- 1. SETUP THE FIGURE ---
    plt.ioff()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.5))
    fig.canvas.header_visible = False
    plt.tight_layout()

    # --- Plot 1: Energy Levels vs Time ---
    # We plot the instantaneous eigenvalues
    ln_E_plus, = ax1.plot([], [], 'b-', linewidth=2, label=r'Upper Eigenstate $|1\rangle_t$')
    ln_E_minus, = ax1.plot([], [], 'r-', linewidth=2, label=r'Lower Eigenstate $|2\rangle_t$')
    # Fill between to show the gap
    poly_gap = ax1.fill_between([], [], [], color='gray', alpha=0.1, label='Energy Gap')

    ax1.set_xlabel('Time (t)')
    ax1.set_ylabel('Energy')
    ax1.set_title('Instantaneous Energy Levels')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left', fontsize=8)
    
    # --- Plot 2: State Composition (Probabilities) ---
    # Probability of finding the system in Bare States if we track the Lower Eigenstate
    ln_prob_g, = ax2.plot([], [], 'orange', linewidth=2, label=r'$|\langle g, n | 2 \rangle_t|^2$ (Ground char)')
    ln_prob_e, = ax2.plot([], [], 'purple', linewidth=2, linestyle='--', label=r'$|\langle e, n | 2 \rangle_t|^2$ (Excited char)')
    
    txt_info = ax2.text(0.5, 0.5, "", transform=ax2.transAxes, 
                        bbox=dict(facecolor='white', alpha=0.8), ha='center')

    ax2.set_xlabel('Time (t)')
    ax2.set_ylabel('Probability (Overlap)')
    ax2.set_title(r'Composition of Adiabatic State $|\psi(t)\rangle$')
    ax2.set_ylim(-0.05, 1.05)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='center right', fontsize=9)

    plt.ion()

    # --- 2. UPDATE LOGIC ---
    def update_plot(alpha, detuning, n, time_range):
        # Time array
        t = np.linspace(0, time_range, 500)
        
        # 1. Calculate Coupling Omega(t) using the input function
        # We ensure it handles numpy arrays
        try:
            Omega_t = omega_func(t, alpha)
        except Exception as e:
            # Fallback for scalar-only functions
            Omega_t = np.array([omega_func(ti, alpha) for ti in t])

        # 2. Calculate Instantaneous Rabi Frequency (The Gap)
        # R_n(t) = sqrt(Delta^2 + 4 * Omega(t)^2 * n)
        # We take absolute value of detuning for safe math, though Delta^2 handles it
        gap_t = np.sqrt(detuning**2 + 4 * Omega_t**2 * n)
        
        # 3. Energies E = +/- Gap/2
        # (Assuming center of mass energy is 0 for simplicity)
        E_plus = 0.5 * gap_t
        E_minus = -0.5 * gap_t
        
        # 4. Mixing Angle theta(t)
        # tan(2*theta) = 2 * Omega * sqrt(n) / Delta
        # We use arctan2 to handle Delta=0 correctly (avoid div by zero)
        # theta ranges from 0 (product state) to pi/4 (maximally entangled)
        theta_t = 0.5 * np.arctan2(2 * Omega_t * np.sqrt(n), detuning)
        
        # 5. Probabilities (assuming we follow the Lower state |->)
        # |-> = cos(theta)|g,n> - sin(theta)|e,n>
        P_ground = np.cos(theta_t)**2
        P_excited = np.sin(theta_t)**2

        # --- UPDATE PLOTS ---
        # Plot 1
        ln_E_plus.set_data(t, E_plus)
        ln_E_minus.set_data(t, E_minus)
        
        # Update fill_between (requires a bit of hack in matplotlib interactive)
        # We can't update data easily, simpler to clear and redraw collection or just ignore for speed
        # For simplicity in widgets, we often skip updating fill_between or use PathCollection hacks.
        # Let's just update limits.
        ax1.set_xlim(0, time_range)
        ax1.set_ylim(np.min(E_minus)*1.1, np.max(E_plus)*1.1)

        # Plot 2
        ln_prob_g.set_data(t, P_ground)
        ln_prob_e.set_data(t, P_excited)
        ax2.set_xlim(0, time_range)
        
        # Adiabaticity Check (Scalar estimate at t=0)
        # Condition: alpha << Delta^2
        if detuning != 0:
            adiabatic_factor = alpha / (detuning**2)
            cond_text = f"Adiabatic Factor $\\alpha / \Delta^2$: {adiabatic_factor:.2f}"
            if adiabatic_factor > 0.5:
                cond_text += "\n(Warning: Too fast!)"
        else:
            cond_text = "Resonant (Gap defined by $\Omega$)"

        txt_info.set_text(cond_text)
        
        fig.canvas.draw_idle()

    # --- 3. CREATE WIDGET ---
    widget = interactive(update_plot, 
             alpha=FloatSlider(value=0.5, min=0.0, max=alpha_max, step=0.01, description=r'Rate alpha'),
             detuning=FloatSlider(value=2.0, min=0.1, max=10.0, step=0.1, description=r'Detuning Delta'),
             n=IntSlider(value=1, min=1, max=10, step=1, description='Photons n'),
             time_range=FloatSlider(value=10.0, min=1.0, max=50.0, step=1.0, description='Max Time')
            )
    
    controls = VBox(widget.children[:-1])
    layout = VBox([fig.canvas, controls])
    
    return layout


import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interactive, FloatSlider, IntSlider, VBox
from scipy.integrate import solve_ivp

import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interactive, FloatSlider, IntSlider, VBox
from scipy.integrate import solve_ivp

import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interactive, FloatSlider, IntSlider, VBox
from scipy.integrate import solve_ivp

# Default linear ramp for backward compatibility
def default_ramp(t, alpha):
    return alpha * t

def run_dynamic_viz(omega_func=default_ramp, alpha_max=5.0, time_max=50.0):
    """
    Compares Ideal Adiabatic vs. Real Schrödinger Dynamics for ARBITRARY Omega(t).
    
    Parameters:
    omega_func: Function f(t, alpha) returning coupling strength.
    alpha_max: Range for the alpha slider.
    time_max: Range for the time slider.
    """
    
    # --- 1. SETUP FIGURE ---
    plt.ioff()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 4.5))
    fig.canvas.header_visible = False
    plt.tight_layout()

    # --- Plot 1: The Coupling Profile Omega(t) ---
    ln_omega, = ax1.plot([], [], 'k-', linewidth=2, label=r'$\Omega(t)$')
    ax1.set_xlabel('Time (t)')
    ax1.set_ylabel('Coupling Strength $\Omega$')
    ax1.set_title('Coupling Profile')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left')

    # --- Plot 2: Probability Dynamics ---
    # The "Ideal" line: The projection of |g,n> onto the INSTANTANEOUS eigenstate
    ln_ideal, = ax2.plot([], [], 'k--', linewidth=1.5, alpha=0.6, label='Adiabatic Prediction')
    
    # The "Real" line: Solving Schrödinger Eq
    ln_real, = ax2.plot([], [], 'orange', linewidth=2, label='Real Dynamics')
    
    ax2.set_xlabel('Time (t)')
    ax2.set_ylabel(r'Ground State Probability $P_g$')
    ax2.set_title('Dynamics: Ideal vs. Real')
    ax2.set_ylim(-0.05, 1.05)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='lower left', fontsize=9)

    plt.ion()

    # --- 2. SOLVER LOGIC ---
    def solve_schrodinger(t_eval, alpha, detuning, n):
        """
        Solves i d(psi)/dt = H(t) psi
        State vector psi = [c_g, c_e]
        """
        def hamiltonian_deriv(t, psi):
            cg = psi[0]
            ce = psi[1]
            
            # Evaluate user function
            Omega = omega_func(t, alpha)
            coupling = Omega * np.sqrt(n)
            
            # Schrödinger Eq: i * d_psi = H * psi
            d_cg = -1j * coupling * ce
            d_ce = -1j * (coupling * cg + detuning * ce)
            
            return [d_cg, d_ce]

        # Initial State: |g,n> -> [1, 0]
        y0 = [1+0j, 0+0j]
        
        # Integrate with high precision to catch fast oscillations
        sol = solve_ivp(hamiltonian_deriv, [t_eval[0], t_eval[-1]], y0, t_eval=t_eval, rtol=1e-6, atol=1e-7)
        return sol.y

    # --- 3. UPDATE FUNCTION ---
    def update_plot(alpha, detuning, n, time_range):
        t = np.linspace(0, time_range, 600)
        
        # 1. Evaluate Omega(t) (Handle scalar or vector input)
        try:
            Omega = omega_func(t, alpha)
        except:
            Omega = np.array([omega_func(ti, alpha) for ti in t])

        # Update Plot 1
        ln_omega.set_data(t, Omega)
        ax1.set_xlim(0, time_range)
        
        # Dynamic Y-limits for Omega plot
        y_max = np.max(Omega) if np.size(Omega) > 0 else 1.0
        y_min = np.min(Omega) if np.size(Omega) > 0 else 0.0
        margin = (y_max - y_min) * 0.1
        if margin == 0: margin = 0.1
        ax1.set_ylim(y_min - margin, y_max + margin)
        
        # 2. Calculate Ideal Adiabatic Path (Projection onto |->)
        # tan(2theta) = 2*Omega*sqrt(n) / Delta
        # Note: arctan2 handles the signs correctly
        theta = 0.5 * np.arctan2(2 * Omega * np.sqrt(n), detuning)
        P_ideal = np.cos(theta)**2
        
        # 3. Calculate Real Dynamics
        psi_sol = solve_schrodinger(t, alpha, detuning, n)
        P_real = np.abs(psi_sol[0])**2
        
        # 4. Update Plot 2
        ln_ideal.set_data(t, P_ideal)
        ln_real.set_data(t, P_real)
        ax2.set_xlim(0, time_range)
        
        fig.canvas.draw_idle()

    # --- 4. WIDGET ---
    widget = interactive(update_plot, 
             alpha=FloatSlider(value=0.5, min=0.0, max=alpha_max, step=0.01, description=r'Param alpha'),
             detuning=FloatSlider(value=2.0, min=0.1, max=10.0, step=0.1, description=r'Detuning Delta'),
             n=IntSlider(value=1, min=1, max=10, step=1, description='Photons n'),
             time_range=FloatSlider(value=15.0, min=1.0, max=time_max, step=1.0, description='Max Time')
            )
    
    controls = VBox(widget.children[:-1])
    layout = VBox([fig.canvas, controls])
    
    return layout