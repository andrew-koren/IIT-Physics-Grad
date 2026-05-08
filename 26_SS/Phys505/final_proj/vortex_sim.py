from typing import Optional
from dataclasses import dataclass, field
import numpy as np
from numpy import pi
from numpy.fft import fftfreq, ifft, fftshift
from scipy.linalg import solve_banded    # for B field solver
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.special import k1, k0
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

from scipy.constants import mu_0, physical_constants    # mu_0 in H/m
phi0 = physical_constants['mag. flux quantum'][0]       # flux V*s (Wb)

import matplotlib.pyplot as plt


# 1. set up material, all parameters as a function of depth
# 2. compute forces
# 3. compute dynamics
# 4. make some cool animations or smth idk

def ppma_dist(x, curve='O'):
    '''
    Crude approximation of the Fig. 1, N-Infused 
    from Hannah Hu et al. arXiv:2509.18564
    '''
    if curve=='O':
        x1, y1 = 8.4, 4.4
        x2, y2 = 83, np.log10(200)
    else:
        x1, y1 = 7.2, 3
        x2, y2 = 16, 2     
    
    m = (y2 - y1) / (x2 - x1)
    c = y1 - m * x1
    
    curve = np.where(x > x2, 10**(y2), 10**(m * x + c))
        
    return curve

def ell_extracted(x):
    O_dist, N_dist = ppma_dist(x), ppma_dist(x, curve='N')
    c = (O_dist + N_dist)/1e6 * 5.56e22 
    r = 14.6 * (1e-9*1e2)  # nm * (m/nm) * (cm/m) -> cm
    return 1/(c*np.pi*r*r)



@dataclass
class nonhomogeneous_sc:
    '''
    sets up material values and simulation grid
    '''
    # bulk params x >> lambda

    lambda_L    : float = 39e-9     # m
    xi0         : float = 38e-9     # m
    rho_n       : float = 1e-9      # ohm-meters # from paper
    omega       : float = 2*pi*1.3e9# GHz
    mfp_bulk    : float = 200e-9    # m
    v0_bulk     : float = 100       # m/s. Biggest guess work value 1-100 m/s.

    # sim params
    # maybe should be 
    B0          : float = 0.100     # T
    Bvp         : Optional[float] = None # T B0>=Bvp causes first flux entry
    M           : Optional[int]   = None # set y grid fidelity
    gaussian    : bool            = False
    interpolate : bool            = True
                                    
    valid = 0

    # depth-dependent vals are loaded in next
    x               : np.ndarray = field(init=False, repr=False) # m depth 
    # mfp             : np.ndarray = field(init=False, repr=False) # m
    xi_s            : np.ndarray = field(init=False, repr=False) # m
    lambda_s        : np.ndarray = field(init=False, repr=False) # m
    kappa_s         : np.ndarray = field(init=False, repr=False) # m
    B_profile       : np.ndarray = field(default=None, repr=False) # unitless profile
    flux_self_field : np.ndarray = field(default=None, repr=False) # T
    gibbs_profile   : np.ndarray = field(init=False, repr=False)

    # def __post_init__(self):
    #     if self.Bvp is not None:
    #         ...


    def load_depth_dependence(self, x:np.ndarray, xi_s:np.ndarray, lambda_s:np.ndarray):
        ''' assumes all math has been done to match boundaries'''
        self.x = x
        # self.mfp = mfp # not necessary here
        self.xi_s = xi_s
        self.lambda_s = lambda_s
        self.kappa_s = lambda_s/xi_s

        # if self.Bvp is not None:
        #     print('overwritting B_vp')

    # def create_depth_dependence(self, x, mfp):
    #     if self.Bvp is None:
    #         raise ValueError(
    #             'Normalization requires surface B_vp be provided'
    #         )

    def get_depth_index(self, depth: float) -> int: 
        return np.argmin(np.abs(self.x - depth))
    
    def solve_B_applied(self): # static d (λ^2 dB) - B = 0 
                                
        h = self.x[1] - self.x[0] # assumes constant step size
        N = len(self.x)
        L2 = self.lambda_s**2
        L2_halfstep = (L2[:-1] + L2[1:]) / 2 # size N-1

        ab = np.zeros((3,N))

        ab[1, 1:-1] = - (L2_halfstep[1:] + L2_halfstep[:-1])/h**2 - 1
        ab[1, 0]  = 1
        ab[1, -1] = 1

        ab[0, 2:]  = L2_halfstep[1:]  / h**2 # M top row is 1,0,...
        ab[2, :-2] = L2_halfstep[:-1] / h**2#  M bottom row is ...,0,1

        rhs = np.zeros(N)
        rhs[0] = 1.0 # B(0)=1, multiply by B(t)

        B_profile = solve_banded(
            (1,1),
            ab,
            rhs
        ) 

        self.B_profile = B_profile # normalized [0,1]

    def solve_B_imag(self, vortex_depth=None, deep=False):

        h = self.x[1] - self.x[0] # assumes constant step size
        if deep:
            L2 = np.concat((self.lambda_s[::-1]**2, self.lambda_s**2))
            offset = len(self.x)
        else:
            L2 = self.lambda_s**2
            offset = 0

        N = len(L2)

        if self.M is None:
            M = 2*len(self.x)
            self.M = M
        else:
            M = self.M

        # ygrid size is set so B(x0, y->ymax) = 0
        # means that cross-comparing y != 0 is inaccurate for now
        vortex_i = np.argmin(np.abs(self.x-vortex_depth)) + offset
        L_local = np.sqrt(L2[vortex_i])
        y_total_gridsize = 2*5*L_local
        dy = y_total_gridsize/M
    
        k = fftfreq(M, d=dy) * 2 * pi 
        
        L2_halfstep = (L2[:-1] + L2[1:]) / 2 # size N-1

        base_ab = np.zeros((3,N))

        base_ab[1, 1:-1] = -(L2_halfstep[1:] + L2_halfstep[:-1])/h**2 - 1.0

        base_ab[0, 1: ] = L2_halfstep / h**2 # M top row is 1,0,...
        base_ab[2, :-1] = L2_halfstep / h**2#  M bottom row is ...,0,1
        
        if self.gaussian:
            xi = self.xi_s[vortex_i-offset]
            if deep:
                dist_sq_x = (np.concat((-self.x[::-1], self.x)) - vortex_depth)**2
            else:
                dist_sq_x = (self.x - vortex_depth)**2

            # rhs is f(x)g(k)
            f_x = np.exp(-dist_sq_x / (2 * xi**2))
            # f_x /= np.sum(f_x * h)  # Normalize so the spatial integral is 1 (* phi_0)
            f_x /= np.sqrt(2*pi)*xi   # alternative norm
            rhs_base = -phi0 * f_x / dy

        else:
            rhs = np.zeros(N)
            rhs[vortex_i] = -phi0/(h*dy)# δ(y) source still contributes
                                        # 1/dy term even after fft       
 
        B_k = np.zeros((N, M), dtype=complex)
        for m in range(M):
            ab = base_ab.copy()
            ab[1, 1:-1] -= k[m]*k[m] * L2[1:-1]
            if self.gaussian:
                rhs = rhs_base * np.exp(-0.5 * (k[m] * xi)**2)

            
            # better than Dirichlet BC - leaks  like e^-kappa x
            kappa = np.sqrt(k[m]**2 + 1.0/L2[-1])
            ab[1, -1] =   1.0 + kappa * h   # Main  diagonal at index N-1
            ab[2, -2] = - 1.0               # Lower diagonal at index N-2

            if deep:
                ab[1, 0] =  1.0 + kappa * h # already mirrored λ(-x) = λ(x)
                ab[0, 1] = -1.0
            else:
                ab[1, 0] = 1.0
                ab[0, 1] = 0.0

            B_k[:, m] = solve_banded((1,1), ab, rhs)

        B_2D = fftshift(ifft(B_k, axis=1), axes=1)
        return B_2D, y_total_gridsize
    
    def get_self_field_profile(self): # biggest bottleneck
        N = len(self.x)

        self_field = np.zeros(N)
        for i, x0 in enumerate(self.x):
            ix = np.argmin(np.abs(self.x-x0))
            assert ix == i # just in case

            total_self_field = self.solve_B_imag(vortex_depth=x0           )[0].real
            vortex_field     = self.solve_B_imag(vortex_depth=x0, deep=True)[0].real
            image_at_x0      = total_self_field[ix, self.M//2] - vortex_field[ix+N, self.M//2]

            self_field[i] = image_at_x0

        self.flux_self_field = self_field


    def gibbs(self, Bt=None, components=False, recompute=True):
        if Bt is None:
            Bt = self.B0

        if recompute | (self.flux_self_field is None):
            self.get_self_field_profile()
        if recompute | (self.B_profile is None):
            self.solve_B_applied()

        applied_field_energy = phi0/mu_0 * self.B_profile*Bt 
        surface_energy = phi0/mu_0 * self.flux_self_field /2
        local_Bc1 = phi0 * np.log(self.kappa_s) / (4*pi*self.lambda_s**2)
        fluxon_energy = phi0/mu_0 * local_Bc1

        energies = (fluxon_energy, surface_energy, applied_field_energy)      
        return energies if components else np.sum(energies, axis=0)
    
    def forces(self, Bt=None, components=False, recompute=True):
        if Bt is None:
            Bt = self.B0

        fluxon_energy, surface_energy, applied_energy = self.gibbs(Bt, components=True, recompute=recompute)
        flux_force      = -np.gradient(fluxon_energy,  self.x)
        surface_force_r = -np.gradient(surface_energy, self.x)
        applied_force   = -np.gradient(applied_energy, self.x)

        if self.interpolate:
            F_interp = CubicSpline(self.x[1:], surface_force_r[1:], extrapolate=True)
            x_eff = np.sqrt(self.x**2 + self.xi_s**2)
            surface_force = F_interp(x_eff)
        else:
            surface_force = surface_force_r

        forces = (flux_force, surface_force, applied_force)
        return forces if components else np.sum(forces, axis=0)


# some ai generated code here
# mostly rewritting / packaging plots
# I didn't have time to pack into this 
# module by hand

class plotter:
    def __init__(self, sc: nonhomogeneous_sc):
        self.sc: nonhomogeneous_sc = sc

    def xscale(self, x=None, x_scale=None):
        """Return (x_display_array, default_xlabel_string)."""
        if x is None:
            x = self.sc.x
        if x_scale is None:
            return x / self.sc.lambda_L, r"Depth  $x\,/\,\lambda_L$"
        f    = float(x_scale)
        unit = {1e9: "nm", 1e6: "μm", 1e3: "mm", 1.0: "m"}.get(f, f"×{f:.2g} m")
        return x * f, f"Depth ({unit})"

    def x_norm(self, x_scale):
        return 1.0 / self.sc.lambda_L if x_scale is None else float(x_scale)

    def y_norm(self, arr, y_norm=False):
        if not y_norm:
            return arr
        finite = arr[np.isfinite(arr)]
        ref = finite[-1] if len(finite) else 1.0
        return arr / ref if ref != 0 else arr

    def setup_ax(self, ax, figsize=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
            return fig, ax
        return None, ax

    def draw_labels(self, ax, fig=None, *,
                  xlabel=None, ylabel=None, title=None,
                  xlim=None, ylim=None, legend=True, grid=True):
        if xlabel:  ax.set_xlabel(xlabel)
        if ylabel:  ax.set_ylabel(ylabel)
        if title:   ax.set_title(title)
        if xlim:    ax.set_xlim(*xlim)
        if ylim:    ax.set_ylim(*ylim)
        if grid:    ax.grid(True, alpha=0.3)
        if legend:  ax.legend()
        if fig is not None:
            plt.tight_layout()
            plt.show()

    def _analytical_B(self, x=None, Bt=None):
        sc = self.sc
        if x   is None: x  = sc.x
        if Bt  is None: Bt = sc.B0
        return Bt * np.exp(-x / sc.lambda_L)

    def _analytical_forces(self, x=None, Bt=None):
        sc = self.sc
        if x  is None: x  = sc.x
        if Bt is None: Bt = sc.B0
        lam, xi0 = sc.lambda_L, sc.xi0

        ana_applied = phi0 * Bt / (mu_0 * lam) * np.exp(-x / lam)
        ana_surface = -(phi0**2 /
                        (2 * pi * mu_0 * lam**3) *
                        k1(2 / lam * np.sqrt(x**2 + xi0**2)))
        ana_flux    = np.zeros_like(x)
        return ana_flux, ana_surface, ana_applied

    def plot_material(self, ax=None, x_scale=None, y_norm=True,
                      ell=None, xlabel=None, ylabel=None,
                      title="Material Parameters"):

        fig, ax = self.setup_ax(ax)
        xd, xl  = self.xscale(x_scale=x_scale)

        curves = {
            r"$\xi_s$":      self.sc.xi_s,
            r"$\lambda_s$":  self.sc.lambda_s,
            r"$\kappa_s$":   self.sc.kappa_s,
        }
        if ell is not None:
            curves[r"$\ell$"] = ell            

        for label, arr in curves.items():
            ax.plot(xd, self.y_norm(arr, y_norm), label=label)

        yl = ylabel or ("Normalised (bulk = 1)" if y_norm else "Value (SI)")
        self.draw_labels(ax, fig,
                       xlabel=xlabel or xl, ylabel=yl, title=title)
        return ax

    # ================================================================ B field

    def plot_B(self, ax=None, x_scale=None, components=False, analytical=False,
               Bt=None, y_norm=False, xlabel=None, ylabel=None,
               safe_slice=slice(None), title="Magnetic Field Profile"):
        sc = self.sc
        if sc.B_profile is None:
            sc.solve_B_applied()
        if Bt is None:
            Bt = sc.B0

        fig, ax = self.setup_ax(ax)
        xd, xl  = self.xscale(x_scale=x_scale)
        sl = safe_slice
        xds = xd[sl]

        def plot(y, label, ls='-', alpha=1.0):
            ax.plot(xds, self.y_norm(np.asarray(y)[sl], y_norm),
                    ls=ls, label=label, alpha=alpha)

        if components:
            plot(sc.B_profile * Bt, "applied (computed)", alpha=0.75)
            if sc.flux_self_field is not None:
                plot(sc.flux_self_field, "self-field (computed)", alpha=0.75)
        plot(sc.B_profile * Bt, "total (computed)", alpha=0.95)

        if analytical:
            plot(self._analytical_B(sc.x, Bt), "applied (analytical)", ls='--', alpha=0.6)

        yl = ylabel or ("B / B(bulk)" if y_norm else "B (T)")
        self.draw_labels(ax, fig,
                       xlabel=xlabel or xl, ylabel=yl, title=title)
        return ax

    # ================================================================ Gibbs free energy

    def plot_gibbs(self, ax=None, x_scale=None, components=True, analytical=False,
                   Bt=None, y_norm=False, xlabel=None, ylabel=None,
                   safe_slice=slice(None), title=None):

        sc = self.sc
        flux, surface, applied = sc.gibbs(Bt=Bt, components=True)
        total = flux + surface + applied

        fig, ax = self.setup_ax(ax)
        xd, xl  = self.xscale(x_scale=x_scale)
        sl  = safe_slice
        xds = xd[sl]

        def plot(y, label, ls='-', alpha=1.0):
            ax.plot(xds, self.y_norm(np.asarray(y)[sl], y_norm),
                    ls=ls, label=label, alpha=alpha)

        if components:
            plot(flux,    "flux (condensation)")
            plot(surface, "image (boundary)")
            plot(applied, "Lorentz (applied)")
        plot(total, "total Gibbs", alpha=0.95)

        if analytical:
            _, a_surf, a_app = self._analytical_forces(sc.x[sl], Bt)
            if components:
                plot(a_app,         "applied (analytical)", ls='--', alpha=0.55)
                plot(a_surf,        "image (analytical)",   ls='--', alpha=0.55)
            else:
                plot(a_app + a_surf,"total (analytical)",   ls='--', alpha=0.55)

        ttl = title or ("Gibbs Energy — components" if components else "Gibbs Free Energy")
        yl  = ylabel or ("G / G(bulk)" if y_norm else r"G  (J m$^{-1}$)")
        self.draw_labels(ax, fig,
                       xlabel=xlabel or xl, ylabel=yl, title=ttl)
        return ax

    # ================================================================ forces

    def plot_force(self, ax=None, x_scale=None, components=True, analytical=False,
                   Bt=None, y_norm=False, xlabel=None, ylabel=None,
                   safe_slice=slice(None), title=None):

        sc = self.sc
        flux, surface, applied = sc.forces(Bt=Bt, components=True, recompute=False)
        total = flux + surface + applied

        fig, ax = self.setup_ax(ax)
        xd, xl  = self.xscale(x_scale=x_scale)
        sl  = safe_slice
        xds = xd[sl]

        def plot(y, label, ls='-', alpha=1.0):
            ax.plot(xds, self.y_norm(np.asarray(y)[sl], y_norm),
                    ls=ls, label=label, alpha=alpha)

        if components:
            plot(flux,    "flux (condensation)")
            plot(surface, "image (boundary)")
            plot(applied, "Lorentz (applied)")
        plot(total, "total force", alpha=0.95)

        if analytical:
            _, a_surf, a_app = self._analytical_forces(sc.x[sl], Bt)
            if components:
                plot(a_app,         "applied (analytical)", ls='--', alpha=0.55)
                plot(a_surf,        "image (analytical)",   ls='--', alpha=0.55)
            else:
                plot(a_app + a_surf,"total (analytical)",   ls='--', alpha=0.55)

        ttl = title or ("Force — components" if components else "Net Vortex Force")
        yl  = ylabel or ("F / F(bulk)" if y_norm else r"Force  (N m$^{-2}$)")
        self.draw_labels(ax, fig,
                       xlabel=xlabel or xl, ylabel=yl, title=ttl)
        return ax

    # ================================================================ error vs analytical

    def plot_error(self, ax=None, x_scale=None, which='force',
                   log=True, Bt=None, safe_slice=slice(None),
                   xlabel=None, ylabel=None, title=None):
        sc = self.sc
        fig, ax = self.setup_ax(ax)
        xd, xl  = self.xscale(x_scale=x_scale)
        sl  = safe_slice
        xds = xd[sl]

        if which == 'B':
            if sc.B_profile is None:
                sc.solve_B_applied()
            if Bt is None: Bt = sc.B0
            computed  = sc.B_profile[sl] * Bt
            reference = self._analytical_B(sc.x[sl], Bt)
            tag = "B field"
        else:
            flux, surf, app = sc.forces(Bt=Bt, components=True)
            computed  = (flux + surf + app)[sl]
            _, as_, aa_ = self._analytical_forces(sc.x[sl], Bt)
            reference = (as_ + aa_)[sl]
            tag = "Force"

        rel_err = np.abs(computed - reference) / (np.abs(reference) + 1e-30) * 100
        ax.plot(xds, rel_err)
        if log:
            ax.set_yscale('log')

        ttl = title or f"{tag} — relative error vs analytical"
        yl  = ylabel or "Relative error (%)"
        self.draw_labels(ax, fig,
                       xlabel=xlabel or xl, ylabel=yl, title=ttl, legend=False)
        return ax

    # ================================================================ vortex 2-D imshow

    def plot_vortex_field(self, ax=None, x_scale=None, vortex_depth=None,
                          deep=True, xlim=None, ylim=None,
                          cmap='RdBu_r', title=None):
        sc = self.sc
        if vortex_depth is None:
            vortex_depth = 2 * sc.lambda_L

        result, yrange = sc.solve_B_imag(vortex_depth=vortex_depth, deep=deep)

        fig, ax = self.setup_ax(ax)
        xd, xl  = self.xscale(x_scale=x_scale)
        f       = self.x_norm(x_scale)

        x_max  = xd[-1]
        y_half = (yrange / 2) * f

        extent = [-x_max if deep else xd[0], x_max, -y_half, y_half]
        data   = result.real.T
        vb     = np.max(np.abs(data))

        im = ax.imshow(data, extent=extent, aspect='auto',
                       vmin=-vb, vmax=vb, cmap=cmap, origin='lower')
        if fig is not None:
            plt.colorbar(im, ax=ax, label="B (T)")

        # apply display limits (easily overridden by caller)
        if xlim is not None:
            ax.set_xlim(*xlim)
        elif deep:
            ax.set_xlim(-3 * sc.lambda_L * f, 3 * sc.lambda_L * f)

        if ylim is not None:
            ax.set_ylim(*ylim)

        yl_unit = xl.split("(")[-1].rstrip(")") if "(" in xl else r"\lambda_L"
        ttl = title or ("Vortex field — deep domain" if deep else "Vortex field — surface domain")
        self.draw_labels(ax, fig,
                       xlabel=xl, ylabel=f"y ({yl_unit})",
                       title=ttl, legend=False, grid=False)
        return ax

    # ================================================================ self-field (two-panel)

    def plot_self_field(self, axes=None, x_scale=None,
                        log=True, xlim_lambda=3, xlabel=None):
        sc = self.sc
        if sc.flux_self_field is None:
            sc.get_self_field_profile()

        lam = sc.lambda_L
        self_field_ana      = -(phi0 / (2 * pi * lam**2)) * k0(2 * sc.x / lam)
        self_field_ana[0]   = self_field_ana[1]   # patch singularity

        xd, xl  = self.xscale(x_scale=x_scale)
        f       = self.x_norm(x_scale)
        xmax    = xlim_lambda * lam * f   # display-unit xlim

        if axes is None:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9), sharex=True)
            show = True
        else:
            fig       = None
            ax1, ax2  = axes
            show      = False

        # ---- panel 1: field magnitudes
        ax1.plot(xd, np.abs(sc.flux_self_field), label="computed",       alpha=0.85)
        ax1.plot(xd, np.abs(self_field_ana),      label="analytical K₀",
                 ls='--', alpha=0.85)
        if log:
            ax1.set_yscale('log')
        ax1.set_xlim(0, xmax)
        ax1.set_ylabel("|B| (T)")
        ax1.set_title("Self-field: computed vs analytical K₀")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # ---- panel 2: percentage error
        rel_err = (np.abs(sc.flux_self_field - self_field_ana) /
                   (np.abs(self_field_ana) + 1e-30))
        ax2.plot(xd, rel_err)
        ax2.set_yscale('log')
        ax2.set_xlim(0, xmax)
        ax2.set_ylim(top=400)
        ax2.set_xlabel(xlabel or xl)
        ax2.set_ylabel("Relative error")
        ax2.set_title("Self-field error vs analytical")
        ax2.grid(True, alpha=0.3)

        if show:
            plt.tight_layout()
            plt.show()
        return ax1, ax2

    # ================================================================ fluxon on-axis cut (y=0)

    def plot_fluxon_cut(self, ax=None, x_scale=None, vortex_depth=None,
                        xlim_lambda=3, xlabel=None, title=None):
        """
        On-axis (y = 0) field decomposition for a single fluxon:
          total field, free-space vortex field, and their difference (image field).

        Calls sc.solve_B_imag twice (surface + deep) so it can be slow.
        """
        sc = self.sc
        if vortex_depth is None:
            vortex_depth = 2 * sc.lambda_L

        result,      _ = sc.solve_B_imag(vortex_depth=vortex_depth)
        result_free, _ = sc.solve_B_imag(vortex_depth=vortex_depth, deep=True)

        N            = len(sc.x)
        total        = result.real[:,      sc.M // 2]
        free_vortex  = result_free.real[N:, sc.M // 2]
        image_field  = total - free_vortex

        fig, ax = self.setup_ax(ax)
        xd, xl  = self.xscale(x_scale=x_scale)
        f       = self.x_norm(x_scale)

        ax.plot(xd, total,       label="total field",  alpha=0.85)
        ax.plot(xd, free_vortex, label="free vortex",  alpha=0.85, ls='--')
        ax.plot(xd, image_field, label="image field",  alpha=0.85, ls=':')

        x0_disp = vortex_depth * f
        ax.axvline(x0_disp, color='grey', alpha=0.25,
                   label=f"vortex @ {x0_disp:.2g}")
        ax.set_xlim(0, xlim_lambda * sc.lambda_L * f)

        ttl = title or "On-axis field decomposition  (y = 0)"
        self.draw_labels(ax, fig,
                       xlabel=xlabel or xl, ylabel="B (T)", title=ttl)
        return ax


