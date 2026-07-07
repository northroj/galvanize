import numpy as np
import matplotlib.pyplot as plt


# ----------------------------
# User inputs
# ----------------------------

case = 2

if case == 1: # cold relativistic rutherford
    # Target material inputs
    z_target = 74
    a_target = 184.0
    density_target = 19.25  # g/cc
else:
    z_target = 1
    z_projectile = 1
    a_target = 2.013553214
    a_projectile = 1.007276
    density_target = 5.0 # g/cc

# Projectile inputs
z_projectile = 1

# Choose projectile rest mass energy in keV.
# For a proton:
projectile_mass_energy_keV = 938272.0813

# Energy grid (kinetic energy, keV)
energy_grid = np.logspace(np.log10(10.0), np.log10(2.0e6), 400)


# ----------------------------
# Physics helpers
# ----------------------------

def beta_sq_from_kinetic_energy(kinetic_energy_keV, mass_energy_keV):
    """
    Relativistic beta^2 from kinetic energy T and rest mass energy mc^2.
    T, mc^2 in keV.
    """
    gamma = 1.0 + kinetic_energy_keV / mass_energy_keV
    return 1.0 - 1.0 / (gamma * gamma)


def mean_excitation_energy_approximation(z):

    if z == 1:
        mei = 19.0*1e-3
    elif z >= 2 and z <= 13:
        mei = (11.2 + 11.7*z) *1e-3
    else:
        mei = (52.8 + 8.71 * z) *1e-3

    return mei  # keV


def rutherford_ion_electron_constant(beta_sq, z_projectile, z_target, a_target, rho_target):
    """
    Matches your C++:
    0.1536*1e3*z_projectile^2*z_target*rho_target / (a_target*beta_sq)
    """
    return 0.1536 * 1.0e3 * z_projectile**2 * z_target * rho_target / (a_target * beta_sq)

def rutherford_ion_ion_constant(z_projectile, z_target, a_projectile, a_target, rho_target, energy_projectile):
    a_ratio = a_target/a_projectile
    alpha_hbar_c = 1.43996448e-4
    value = z_projectile * z_target * alpha_hbar_c * (1.0+a_ratio) / (2.0*energy_projectile*a_ratio)
    n_target = rho_target*6.022E23 / a_target
    return n_target * value * value

def cutoff_fp_mu_min(a_ratio):
    mu_min = 0.0
    eps = 1e-10
    if (a_ratio > 1.0 + eps):
        return -1.0
    elif (np.abs(a_ratio - 1.0)<= eps):
        return np.sqrt(2.0) / 2.0
    else:
        val = 1.0 - a_ratio*a_ratio
        if val < 0:
            val = 0
        return np.sqrt(val)



def rutherford_dcs_moment(kinetic_energy_keV, z_projectile, z_target, a_target, rho_target, order):
    """
    Mirrors your C++ rutherford_dcs_moment for the ion-electron case.
    """
    beta_sq = beta_sq_from_kinetic_energy(kinetic_energy_keV, projectile_mass_energy_keV)

    if case == 1:

        prefix = rutherford_ion_electron_constant(
            beta_sq, z_projectile, z_target, a_target, rho_target
        )

        q_min = mean_excitation_energy_approximation(z_target)  # keV
        q_max = 1.022 * 1.0e3 * beta_sq / (1.0 - beta_sq)      # keV

        # Guard against invalid log/division regions
        if q_max <= q_min or beta_sq <= 0.0:
            return np.nan

        if order == 0:  # scattering cross section [cm^-1]
            suffix = (1.0 / q_min - 1.0 / q_max) - beta_sq / q_max * np.log(q_max / q_min)
        elif order == 1:  # stopping power [keV/cm]
            suffix = np.log(q_max / q_min) - beta_sq * (1.0 - q_min / q_max)
        else:  # straggling n = 2 and higher
            n = float(order)
            suffix = (
                1.0 / (n - 1.0)
                * (
                    q_max**(n - 1.0) * (1.0 - beta_sq * ((n - 1.0) / n))
                    - q_min**(n - 1.0) * (1.0 - beta_sq * ((n - 1.0) / n) * q_min / q_max)
                )
            )

    elif case == 2:
        a_ratio = a_target/a_projectile
        mu_cut = 1.0-1e-4
        mu_min = cutoff_fp_mu_min(a_ratio)
        mu_min = 0.95
        t_factor = 2.0 * kinetic_energy_keV * a_ratio / ((1.0 + a_ratio)*(1.0 + a_ratio))
        q_min = (1.0 - mu_cut) * t_factor
        q_max = (1.0 - mu_min) * t_factor
        prefix = t_factor * rutherford_ion_ion_constant(z_projectile, z_target, a_projectile, a_target, rho_target, kinetic_energy_keV)
        if (order == 0):
            suffix = (1.0/q_min - 1.0/q_max)
        elif (order == 1):
            suffix = np.log(q_max/q_min)
        elif (order > 1):
            suffix = (1.0/(order - 1.0)) * (np.pow(q_max,order-1) - np.pow(q_min,order-1))

    return prefix * suffix


# ----------------------------
# GFP cross section builders
# ----------------------------

def xs_gfp2(kinetic_energy_keV, z_projectile, z_target, a_target, rho_target):
    """
    Matches your C++ select_scattering branch for rutherford_order == 2
    xs = alpha / beta
    where:
        alpha = Q1
        beta  = Q2 / (2 Q1)
    """
    Q1 = rutherford_dcs_moment(
        kinetic_energy_keV, z_projectile, z_target, a_target, rho_target, 1
    )
    Q2 = rutherford_dcs_moment(
        kinetic_energy_keV, z_projectile, z_target, a_target, rho_target, 2
    )

    if not np.isfinite(Q1) or not np.isfinite(Q2) or Q1 == 0.0:
        return np.nan

    alpha = Q1
    beta = Q2 / (2.0 * Q1)

    if beta == 0.0:
        return np.nan

    return alpha / beta


def xs_gfp4(kinetic_energy_keV, z_projectile, z_target, a_target, rho_target):
    """
    Matches your C++ select_scattering branch for rutherford_order == 4
    xs = alpha1/beta1 + alpha2/beta2
    """
    Q1 = rutherford_dcs_moment(
        kinetic_energy_keV, z_projectile, z_target, a_target, rho_target, 1
    )
    Q2 = rutherford_dcs_moment(
        kinetic_energy_keV, z_projectile, z_target, a_target, rho_target, 2
    )
    Q3 = rutherford_dcs_moment(
        kinetic_energy_keV, z_projectile, z_target, a_target, rho_target, 3
    )
    Q4 = rutherford_dcs_moment(
        kinetic_energy_keV, z_projectile, z_target, a_target, rho_target, 4
    )

    vals = [Q1, Q2, Q3, Q4]
    if not all(np.isfinite(v) for v in vals):
        return np.nan

    sqrt3 = np.sqrt(3.0)

    rad_alpha = (
        4.0 * Q3 * Q3 * (-9.0 * Q2 * Q2 + 8.0 * Q1 * Q3)
        + 36.0 * Q2 * (Q2 * Q2 - Q1 * Q3) * Q4
        + 3.0 * Q1 * Q1 * Q4 * Q4
    )

    rad_beta = (
        12.0 * Q3 * Q3 * (-9.0 * Q2 * Q2 + 8.0 * Q1 * Q3)
        + 108.0 * Q2 * (Q2 * Q2 - Q1 * Q3) * Q4
        + 9.0 * Q1 * Q1 * Q4 * Q4
    )

    denom_alpha = 2.0 * np.sqrt(rad_alpha) if rad_alpha > 0.0 else np.nan
    denom_beta = 36.0 * Q2 * Q2 - 24.0 * Q1 * Q3

    if not np.isfinite(denom_alpha) or denom_alpha == 0.0 or denom_beta == 0.0:
        return np.nan

    sqrt_rad_alpha = np.sqrt(rad_alpha)
    sqrt_rad_beta = np.sqrt(rad_beta) if rad_beta >= 0.0 else np.nan

    if not np.isfinite(sqrt_rad_beta):
        return np.nan

    alpha1 = (
        -6.0 * sqrt3 * Q2 * Q2 * Q2
        + 6.0 * sqrt3 * Q1 * Q2 * Q3
        + Q1 * (-sqrt3 * Q1 * Q4 + sqrt_rad_alpha)
    ) / denom_alpha

    alpha2 = (
        6.0 * sqrt3 * Q2 * Q2 * Q2
        - 6.0 * sqrt3 * Q1 * Q2 * Q3
        + Q1 * (sqrt3 * Q1 * Q4 + sqrt_rad_alpha)
    ) / denom_alpha

    beta1 = -(-6.0 * Q2 * Q3 + 3.0 * Q1 * Q4 + sqrt_rad_beta) / denom_beta
    beta2 = (6.0 * Q2 * Q3 - 3.0 * Q1 * Q4 + sqrt_rad_beta) / denom_beta

    if beta1 == 0.0 or beta2 == 0.0:
        return np.nan

    return alpha1 / beta1 + alpha2 / beta2


# ----------------------------
# Evaluate over the grid
# ----------------------------

xs_order2 = np.array([
    xs_gfp2(E, z_projectile, z_target, a_target, density_target)
    for E in energy_grid
])

xs_order4 = np.array([
    xs_gfp4(E, z_projectile, z_target, a_target, density_target)
    for E in energy_grid
])


# ----------------------------
# Plot
# ----------------------------

plt.figure(figsize=(8, 6))
plt.plot(energy_grid, xs_order2, label="GFP2")
plt.plot(energy_grid, xs_order4, label="GFP4")
plt.yscale("log")
plt.xscale("log")

plt.xlabel("Projectile kinetic energy [keV]")
plt.ylabel("Scattering cross section [cm$^{-1}$]")
if case == 1:
    plt.title("GFP Relativistic Rutherford cross section protons on tungsten")
elif case == 2:
    plt.title("GFP Cutoff cross section protons on deuterium")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()