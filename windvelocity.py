import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import interp1d

#read file in 

filename = "BoOST Data files//f020-500.izw18.dat"


data = np.genfromtxt(filename, autostrip=True)


class StellarModel:
        def __init__(self, data, filename):
            # Physical constant
            self.G = scipy.constants.G

            # Time in years from data
            self.time = data[:, 0]  # years

            # Initial mass from filename string
            self.initialmass = filename[1:4]

            # Convert time from years to Myr
            self.time_myr = self.time / 1e6  # Myr

            # Stellar parameters from columns
            self.mass = data[:, 1]  # M/Msun
            self.temp = data[:, 2]  # K
            self.luminosity = data[:, 3]  # log(L/Lsun)
            self.radius = data[:, 4]  # R/Rsun
            self.mdot = data[:, 5]  # log(Mdot[Msun/yr])
            self.logg = data[:, 6]  # surface gravity
            self.v_surf = data[:, 7]  # km/s
            self.v_crit = data[:, 8]  # km/s
            self.gamma = data[:, 9]  # Eddington factor

            # Hydrogen and helium mass fractions
            self.H_massfrac = data[:, 26] + data[:, 27]
            self.He_massfrac = data[:, 28] + data[:, 29]
            self.He_coremassfrac = data[:, 62] + data[:, 63]

            # Derived metallicity quantities
            self.metallicity = 1 - self.H_massfrac - self.He_massfrac
            self.initial_metallicity = 1 - data[0, 26] - data[0, 27]

            # Carbon and oxygen mass fractions
            self.C_massfrac = data[:, 37] + data[:, 38] + data[:, 39]
            self.O_massfrac = data[:, 43] + data[:, 44] + data[:, 45]

vars = StellarModel(data, filename)

def compute_vwind_vink(mass, radius, Gamma, temp):
    """
    Compute v_wind for an array of timesteps using escape velocity scaling.

    Parameters:
    - mass: array of stellar masses in solar masses
    - radius: array of stellar radii in solar radii
    - Gamma: array of Eddington factors unitless
    - temp: array of surface temperatures in kelvin
    Returns:
    - v_wind: array of wind velocities [km/s]
    """

    # Constants
    G = 6.67430e-8               # [cm^3 g^-1 s^-2]
    M_sun = 1.98847e33           # [g]
    R_sun = 6.957e10             # [cm]

    # Convert to CGS
    M_cgs = mass * M_sun
    R_cgs = radius * R_sun


    f_esc = np.where(temp >= 21000, 2.6, 1.3)

    # Calculate escape velocity
    v_esc = np.sqrt((2 * G * M_cgs * (1 - Gamma)) / R_cgs)  # [cm/s]
    v_esc_kms = v_esc / 1e5  # [km/s]

    # Final wind speed
    v_wind_vink = f_esc * v_esc_kms
    v_wind_vink_cgs = f_esc * v_esc
    
    return v_wind_vink, v_wind_vink_cgs

def compute_vwind_nugis(mass, radius, gamma, luminosity, He_massfrac, C_massfrac, O_massfrac):
    '''
    ***DOESN'T REALLY WORK, RETURNS UNREALISTICALLY LOW VALUES***
    
    Compute v_wind using Nugis & Lamers (2000) prescription for WR stars.

    Parameters:
    - mass: array of stellar masses [M_sun]
    - radius: array of stellar radii [R_sun]
    - gamma: array of Eddington factors
    - luminosity: array of log(L/L_sun)
    - He_massfrac: surface Helium mass fraction
    - C_massfrac: surface Carbon mass fraction
    - O_massfrac: surface Oxygen mass fraction

    Returns:
    - v_wind: array of wind velocities [km/s]
    '''

    # Physical constants
    G = 6.67430e-8               # cm^3 g^-1 s^-2
    M_sun = 1.98847e33           # g
    R_sun = 6.957e10             # cm

    # Convert to CGS
    M_cgs = mass * M_sun
    R_cgs = radius * R_sun

    # Escape velocity with Eddington factor
    v_esc = np.sqrt((2 * G * M_cgs * (1 - gamma)) / R_cgs)  # [cm/s]
    v_esc_kms = v_esc / 1e5  # [km/s]

    # Effective metallicity for WR wind scaling (only C and O matter)
    Z_eff = C_massfrac + O_massfrac
    Z_eff_clipped = np.clip(Z_eff, 0.005, 0.04)  # Valid range from Nugis & Lamers

    # WR phase classification: WN if (C+O)/He < 0.03
    CO_to_He = (C_massfrac + O_massfrac) / He_massfrac
    is_WN = CO_to_He < 0.03
    is_WC = ~is_WN

    # WN exponent
    WN_exponent = 0.61 - 0.13 * luminosity + 0.3 * np.log10(He_massfrac)

    # WC exponent (uses log Z_eff_clipped)
    WC_exponent = -2.37 + 0.13 * luminosity - 0.07 * np.log10(Z_eff_clipped)

    # Combine
    exponent = np.where(is_WN, WN_exponent, WC_exponent)
    v_wind = np.power(10, exponent) * v_esc_kms

    return v_wind

def plotting(time_myr, v_wind_vink, mdot,radius, initialmass):

             
    fig, (ax1, ax2, ax3 ) = plt.subplots(nrows=3, figsize=(8, 6), sharex=True)
   
    #finds late-time values
    totallife = time_myr[-1]
    cutoff_time = totallife * (2/3)
    start_index = np.searchsorted(time_myr, cutoff_time, side='left')

    time_late = time_myr[start_index:]
    v_wind_late_vink = v_wind_vink[start_index:]

    # Plot vwind
    ax1.plot(time_myr, v_wind_vink, color="steelblue", label = "vink")
    ax1.set_ylim()
    ax1.set_ylabel("Wind velocity (km/s)")
    ax1.set_title(f"Initial Mass: {initialmass} Solar Masses")
    ax1.grid(True)

    # Plot mdot
    ax2.plot(time_myr, mdot, color="darkred")
    ax2.set_ylabel("Mass loss rate (log(M☉/yr))")
    ax2.grid(True)

    #plot radius
    ax3.plot(time_myr, radius, color="green")
    ax3.set_xlabel("Time (Myrs)")
    ax3.set_ylabel("Radius (R/Rsun)")
    ax3.grid(True)

    plt.tight_layout()
    plt.savefig(f"Wind Velocity Plots/{initialmass}.png")


    #plot late-times wind velocity
    plt.figure()
    plt.title("Wind Velocity for the last 1/3 of lifetime")
    plt.plot(time_late, v_wind_late_vink, color = "steelblue")
    plt.ylabel("Wind Velocity (km/s)")
    plt.xlabel("Time (Myr)")
    plt.savefig(f"Wind Velocity Plots/latetimes/{initialmass}latetimes.png")
    return

def metal_plotting(time_myr, H_massfrac, He_massfrac, metallicity, initialmass):
    #metallicity plotting
    plt.figure()
    plt.plot(time_myr, H_massfrac, color = "Blue", label = "H Mass Fraction")
    plt.plot(time_myr, He_massfrac, color = "Green", label = "He Mass Fraction")
    plt.plot(time_myr, metallicity, color = "Red", label = "Metallicity")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Mass Fraction")
    plt.legend()
    plt.title(f"Initial Mass: {initialmass} Solar Masses")
    plt.savefig(f"Metallicity plots/{initialmass}.png")


v_wind = compute_vwind_vink(vars.mass, vars.radius, vars.gamma, vars.temp)[1]

