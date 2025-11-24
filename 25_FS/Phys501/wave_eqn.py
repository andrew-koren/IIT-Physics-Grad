import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.special import jv, jn_zeros

# ==========================================
# CONFIGURATION VARIABLES
# ==========================================
FILENAME = "animations/100_mode_3.3.6drumhead_vibration.webp"
FPS = 30              # Frames per second
DURATION = 4*np.pi    # Length of animation in seconds
N_MODES = 100         # Number of terms (n) in the summation
RESOLUTION = 150      # Grid resolution (lower = faster rendering)
RADIUS_A = 1.0        # Radius of the drum (a)
WAVE_SPEED_C = 1.0    # Speed of the wave (c)
Z_LIMITS = (-5.0, 5.0)# Fixed vertical axis limits

# ==========================================
# MATHEMATICAL SETUP
# ==========================================

# 1. Pre-calculate Bessel Zeros (roots) and Coefficients
# We need the first N zeros of J0.
roots = jn_zeros(0, N_MODES)

# Calculate the coefficients: 2 / (x_n * J1(x_n))
# These are constant for the simulation.
coeffs = 2 / (roots * jv(1, roots))

# 2. Create the Grid (Polar -> Cartesian)
r = np.linspace(0, RADIUS_A, RESOLUTION)
theta = np.linspace(0, 2 * np.pi, RESOLUTION)
R, THETA = np.meshgrid(r, theta)

# Convert to Cartesian for plotting
X = R * np.cos(THETA)
Y = R * np.sin(THETA)

# Pre-calculate Spatial terms J0(r * xn / a) to save time in the loop
# Shape: (N_MODES, RESOLUTION, RESOLUTION)
spatial_modes = np.array([jv(0, R * xn / RADIUS_A) for xn in roots])

# ==========================================
# ANIMATION SETUP
# ==========================================

def calculate_displacement(t):
    """
    u(r, t) = sum [ coeff * spatial_mode * sin(c * t * xn / a) ]
    """
    z_total = np.zeros_like(R)
    
    for n in range(N_MODES):
        xn = roots[n]
        # Temporal term: sin(c * t * xn / a)
        temporal = np.sin(WAVE_SPEED_C * t * xn / RADIUS_A)
        
        # Combine: coeff * J0(...) * sin(...)
        z_total += coeffs[n] * spatial_modes[n] * temporal
        
    return z_total

# Create Figure
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

def update(frame):
    t = frame / FPS
    Z = calculate_displacement(t)
    C = np.abs(Z)
    if C.max() == C.min():
        colors = plt.cm.viridis(np.full_like(C, 0.0))
    else:
        norm = (C - C.min()) / (C.max() - C.min())
        colors = plt.cm.viridis(norm)

    # 3D Surface animations in Matplotlib require clearing and redrawing
    ax.clear()
    ax.set_axis_off() # Hide axes for cleaner look
    ax.set_zlim(Z_LIMITS)
    ax.set_xlim(-RADIUS_A, RADIUS_A)
    ax.set_ylim(-RADIUS_A, RADIUS_A)
    
    # Plot the surface
    # cmap 'viridis', 'plasma', 'inferno', 'magma', or 'coolwarm' work well
    surf = ax.plot_surface(X, Y, Z, facecolors=colors, rstride=2, cstride=2, antialiased=True)
    
    ax.set_title(f"Drumhead Vibration\nTime t={t:.2f}", fontsize=14)
    return surf,

# ==========================================
# RENDER AND SAVE
# ==========================================
print(f"Generating animation ({int(DURATION*FPS)} frames)...")

anim = FuncAnimation(
    fig, 
    update, 
    frames=int(DURATION * FPS), 
    interval=1000/FPS
)

# Save using Pillow writer for WebP
anim.save(FILENAME, writer='pillow', fps=FPS)

print(f"Done! Saved as {FILENAME}")