import numpy as np
from scipy.special import jv, jn_zeros
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline
from scipy.ndimage import map_coordinates
import matplotlib.pyplot as plt
from dataclasses import dataclass

def naive_f_expansion(f, R, N=5, M=5, Nr=200, Ntheta=360, fixed_order=None):
    r = np.linspace(0, R, Nr)
    theta = np.linspace(0, 2*np.pi, Ntheta)
    
    rr, tt = np.meshgrid(r, theta, indexing='ij')
    
    # samples Nradius*Ntheta times
    f_vals = f(rr, tt, R)  
    a_nm = np.zeros((N, M), dtype=complex)

    for n in range(N):
        b_order = n if fixed_order is None else fixed_order

        alpha_nm = jn_zeros(b_order, M) 
        k_nm = alpha_nm / R      
        J_vals = jv(b_order, np.outer(r, k_nm)) 
        
        # Angular Integration
        exp_n_theta = np.exp(-1j * n * theta)
        f_theta = f_vals * exp_n_theta[None, :] 
        int_theta = simpson(f_theta, x=theta, axis=1)
        
        # Radial Integration
        integrand = int_theta[:, None] * J_vals * r[:, None]
        integral = simpson(integrand, x=r, axis=0)

        # Normalization 
        norm = (np.pi * R**2 / 2) * (jv(b_order + 1, alpha_nm)**2)
        a_nm[n, :] = integral / norm

    return a_nm

def recon_f_from_naive(a_nm, R, r, theta, fixed_order=None):
    rr, tt = np.meshgrid(r, theta, indexing='ij')
    f_recon = np.zeros_like(rr, dtype=complex)

    N, M = a_nm.shape

    for n in range(N):
        b_order = n if fixed_order is None else fixed_order

        alpha_nm = jn_zeros(b_order, M)
        k_nm = alpha_nm / R
        
        J_vals = jv(b_order, np.outer(r, k_nm))
        exp_n_theta = np.exp(1j * n * theta)
        
        f_n = J_vals @ a_nm[n, :] 
        f_recon += f_n[:, None] * exp_n_theta[None, :]

    return f_recon


def plot_3_spheres(f, R, Nr, Ntheta, N, M, fixed_order):

    # parameters
    r = np.linspace(0, R, Nr)
    theta = np.linspace(0, 2*np.pi, Ntheta)

    a_nm = naive_f_expansion(f, R, N, M, Nr, Ntheta)
    f_recon = recon_f_from_naive(a_nm, R, r, theta)

    a_nm_fixed = naive_f_expansion(f, R, N, M, Nr, Ntheta, fixed_order=fixed_order)
    f_recon_fixed = recon_f_from_naive(a_nm, R, r, theta,  fixed_order=fixed_order)

    rr, tt = np.meshgrid(r, theta, indexing='ij')
    f_true = f(rr, tt, R)

    alldata = [f_true, f_recon, f_recon_fixed]
    vmin = np.real(np.min(alldata))
    vmax = np.real(np.max(alldata))

    fig, ax = plt.subplots(1, 3, figsize=(15, 4), subplot_kw={'projection': 'polar'})
    ax[0].pcolormesh(tt, rr, np.real(f_true), shading='auto', vmin=vmin, vmax=vmax)
    ax[0].set_title("Original f(r,θ)")
    ax[0].set_axis_off()
    ax[1].pcolormesh(tt, rr, np.real(f_recon), shading='auto', vmin=vmin, vmax=vmax)
    ax[1].set_title(r"Textbook reconstruction $J_n$ corresponds to $e^{in\theta}$")
    ax[1].set_axis_off()
    ax[2].pcolormesh(tt, rr, np.real(f_recon_fixed), shading='auto', vmin=vmin, vmax=vmax)
    ax[2].set_title(rf"Fixed $J{fixed_order}$")
    ax[2].set_axis_off()
    plt.show()

### --- Best DFT implementation --- ###

# Zhou & Grisouard implementation 2210.09736v4

# python likes storing variables in a class
@dataclass
class DHT:
    q:      int            # order of bessel functions
    N:      int            # radial resolution
    k:      np.ndarray     # length N (first N zeros)
    kp1:    float          # (N+1)-th zero
    rq:     np.ndarray     # pseudospectral grid r_{q,i} = k_i / k_{N+1}
    w:      np.ndarray     # quadrature weights w_{q,i}
    M:      np.ndarray     # DHT matrix (NxN)

    def forward_transform(self, f_rq): # f_rq is the samples f(rq)
        return self.M @ f_rq
    
    def reverse_transform(self, F):
        return self.kp1**2 * self.M @ F
    
    def cardinal_basis(self, r): # interpolate using bessel functions
        q, k, kp1 = self.q, self.k, self.kp1
        r = np.asarray(r)
        factor = (2.0 * k[:,None]) / (k[:, None]**2 - (r[None,:]**2) * kp1**2)
        return factor * ( jv(q, r*kp1)[None,:] / jv(q+1,k)[:,None] )

# 1D pseudospectral Hankel transform, to be paired with Fourier Transform.
def make_DHT(order, N): 
    bessel_zeros=  jn_zeros(order, N+1)
    k           = bessel_zeros[:N].copy()
    kp1         = bessel_zeros[N]
    rq          = k / kp1  # pseudospectral grid
    w           = 2.0 / (kp1**2 * jv(order+1, k)**2) # weights
    
    # DHT transform matrix M_ij = J(k_i r_j) w_i, where J(k_i) = 0
    k_ij = np.outer(k,rq)
    M = jv(order, k_ij) * w[None, :]

    return DHT(q = order, N = N, k = k, kp1 = kp1,  rq = rq, w = w, M = M)

class Double_Transform:
    def __init__(self, Nr, Nth):
        """
        Nr (N): Number of Bessel modes.
        Nth (M): Number of Fourier modes (must be even).
        """
        self.Nr = Nr
        self.Nth = Nth
        
        # real space sampling - simple meshgrid
        self.r_uniform = np.linspace(0, 1, 2 * self.Nr) 
        self.th_uniform = np.linspace(0, 2*np.pi, self.Nth, endpoint=False)
        self.R, self.TH = np.meshgrid(self.r_uniform, self.th_uniform, indexing='xy')

        # Precompute DHTs for each frequency q (same as order q of DHT)
        self.freqs = np.fft.fftfreq(self.Nth, d=1/self.Nth).astype(int)
        self.dhts = {} # dict not set
        self.cardinal_matrices = {}
        
        for q in np.unique(np.abs(self.freqs)):
            dht = make_DHT(q, self.Nr)
            self.dhts[q] = dht
            self.cardinal_matrices[q] = dht.cardinal_basis(self.r_uniform)

    def phys_to_spec(self, f_data):
        """
        Forward Transform: Physical(r, theta) -> Spectral(q, j)
        """
        # FFT along theta axis
        f_hat_r = np.fft.fft(f_data, axis=0)

        a_nm = np.zeros((self.Nth, self.Nr), dtype=complex)
        # DHT of FFT
        # freqs matches q
        for i, q in enumerate(self.freqs): # must do in loop since q varies
            q_abs = abs(q)
            dht = self.dhts[q_abs]
            
            # i is the radial slice for mode q
            frequency_slice = f_hat_r[i, :]
            
            # Interpolate: Uniform Grid -> DHT Grid (r_q,i)
            # Sec 2.3.1 use cubic spline
            cs = CubicSpline(self.r_uniform, frequency_slice) # will magically allow sampling outside frequency 
            f_rq = cs(dht.rq)
            
            # Perform DHT
            a_nm[i, :] = dht.forward_transform(f_rq)
            
        return a_nm

    def spec_to_phys(self, a_nm):
        """
        Backward Transform: Spectral(q, j) -> Physical(r, theta)
        """
        # rebuild k-space from transform
        f_hat_r = np.zeros((self.Nth, len(self.r_uniform)), dtype=complex)
        for i, q in enumerate(self.freqs): # again, one frequency at a time
            dht = self.dhts[abs(q)]
            
            # reverse DHT 
            a_q = a_nm[i, :]
            f_on_dht_grid = dht.reverse_transform(a_q)
            
            # Interpolate: DHT Grid -> Uniform Grid
            C_matrix = self.cardinal_matrices[abs(q)]
            f_hat_r[i, :] = f_on_dht_grid @ C_matrix

        f_data = np.fft.ifft(f_hat_r, axis=0)
        
        return f_data
    
# got some help from gemini 3 for interpolation
def sample_image_on_polar_grid(img_arr, Double_Transform):
    """
    Maps a cartesian image onto the solver's (N_theta, 2*N_radial) grid.
    """
    h, w = img_arr.shape
    cy, cx = h / 2, w / 2
    max_radius_pix = min(h, w) / 2

    # Transform grid
    R_grid = Double_Transform.R
    TH_grid = Double_Transform.TH

    # Get points from image via cubic spline interpolation, same as FFT step in paper
    X_polar = R_grid * np.cos(TH_grid) * max_radius_pix + cx
    Y_polar = R_grid * np.sin(TH_grid) * max_radius_pix + cy
    coords = np.vstack((Y_polar.flatten(), X_polar.flatten()))

    sampled_flat = map_coordinates(img_arr, coords, order=3, mode='constant', cval=0.0)  # cval = 0 when sampling outside the image
    sampled_polar = sampled_flat.reshape(R_grid.shape)
    sampled_polar[R_grid > 1.0] = 0 # just to be sure
    
    return sampled_polar

def spectrum_plot(ax, FFT_DHT, a_nm, phase = False):
    # Shift coefficients so q=0 is in the center of the x-axis
    coeffs_shifted = np.fft.fftshift(a_nm, axes=0)
    extent = [-FFT_DHT.Nth//2, FFT_DHT.Nth//2, FFT_DHT.Nr, 0] 


    if phase:
        # Plot Phase
        spec_phase = np.angle(coeffs_shifted)
        im = ax.imshow(spec_phase.T, cmap='twilight', aspect='auto', 
                       extent=extent, vmin=-np.pi, vmax=np.pi)
        ax.set_title("Spectrum Phase (rad)")
        ax.set_xlabel("Azimuthal Mode (q)")
        ax.set_yticks([]) # Hide y-axis as it's usually aligned with the magnitude plot
    else:
        # Plot Magnitude
        spec_mag = np.log1p(np.abs(coeffs_shifted))
        vmax = np.percentile(spec_mag, 99.5) 
        im = ax.imshow(spec_mag.T, cmap='inferno', aspect='auto', 
                       extent=extent, vmax=vmax)
        ax.set_title(r"Spectrum Amplitude $\log(1+|F|)$")
        ax.set_xlabel("Azimuthal Mode (q)")
        ax.set_ylabel("Radial Mode (j)")

    # Unified colorbar call
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

def image_visualization(FFT_DHT, img_arr):
    polar_image = sample_image_on_polar_grid(img_arr, FFT_DHT)
    coeffs = FFT_DHT.phys_to_spec(polar_image)
    reconstructed_polar = FFT_DHT.spec_to_phys(coeffs) # in polar coordinates (r, theta)

    fig, ax = plt.subplots(2, 3, figsize=(18, 10))

    ax[0, 0].imshow(img_arr, cmap='gray')
    ax[0, 0].set_title("Original Image (Cartesian)")
    ax[0, 0].axis('off')

    ax[1,0].imshow(polar_image.T, cmap='gray', aspect='auto')
    ax[1,0].set_ylabel('Radius (j)')
    ax[1,0].set_xlabel('Theta (q)')
    ax[1,0].set_title("Unwrapped Input (Polar Domain)")

    ax[0, 2].imshow(reconstructed_polar.real.T, cmap='gray', aspect='auto')
    ax[0, 2].set_title("Unwrapped Reconstruction")
    ax[0, 2].set_xlabel('Theta (q)')


    spectrum_plot(ax[0, 1], FFT_DHT, coeffs, phase=False)
    spectrum_plot(ax[1, 1], FFT_DHT, coeffs, phase=True)

    ax[1, 2].remove() # need to make polar projection
    ax[1, 2] = fig.add_subplot(2, 3, 6, projection='polar')

    r_grid = FFT_DHT.r_uniform
    th_grid = FFT_DHT.th_uniform
    R_mesh, TH_mesh = np.meshgrid(r_grid, th_grid)

    # Plot
    ax[1, 2].pcolormesh(TH_mesh, R_mesh, reconstructed_polar.real, cmap='gray', shading='auto')
    ax[1, 2].set_theta_direction(-1)
    ax[1, 2].set_theta_zero_location("E")
    ax[1, 2].set_title("Reconstruction Projected")
    ax[1, 2].set_axis_off()

    plt.tight_layout()
    plt.show()