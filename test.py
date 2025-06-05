import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from windvelocity import StellarModel
from windvelocity import compute_vwind_vink
filename = "BoOST Data files//f020-500.izw18.dat"
data = np.genfromtxt(filename, autostrip=True)

vars = StellarModel(data, filename)

# ==== INPUT DATA ====
time_myr = vars.time_myr  # Time in Myr

mdot_msunyr = 10**vars.mdot     # Mass-loss rate in Msun/yr
#constantmdot = np.full(np.shape(mdot_msunyr), 5)

v_wind = compute_vwind_vink(vars.mass, vars.radius, vars.gamma, vars.temp)[1]  # Wind velocity in km/s
#constantvwind = np.full(np.shape(mdot_msunyr), 5000)

radius = vars.radius     # Stellar radius in R_sun
#constantrad = np.zeros(np.shape(mdot_msunyr))


# Constants
Msun = 1.989e33  # g
Rsun = 6.957e10  # cm
year = 3.154e7   # s

#Unit Conversions

time_yr = time_myr * 1e6
vwind = v_wind[:-1] * 1e5       # cm/s

# Final time (snapshot)
t_final = time_yr[-1]

#initializing list
shell_radii = []
shell_masses = []

# Loop over each time interval
for i in range(len(time_yr) - 1):
    dt = 100 #time_yr[i+1] - time_yr[i]                # years
    if dt <= 0:
        continue  # skip bad or duplicate time steps

    v = vwind[i]                                       # cm/s
    mdot = mdot_msunyr[i]* Msun / year                 # g/s
    t_release = time_yr[i]
    travel_time = t_final - t_release                  # years
    
    if travel_time <= 0:
        continue  # shell hasn't had time to propagate

    r = v * travel_time * year                         # cm
    m = mdot * dt * year                               # g

    shell_radii.append(r)
    shell_masses.append(m)

shell_radii = np.array(shell_radii)
shell_masses = np.array(shell_masses)

# Bin the shell masses by radius
r_min = 1e20
r_max = 1e29
n_bins = 1000

r_bins = np.logspace(np.log10(r_min), np.log10(r_max), n_bins + 1)
r_centers = np.sqrt(r_bins[:-1] * r_bins[1:])  # geometric mean of bin edges

hist, _ = np.histogram(shell_radii, bins=r_bins, weights=shell_masses)
shell_volumes = (4/3) * np.pi * (r_bins[1:]**3 - r_bins[:-1]**3)
density = hist / shell_volumes


# Plot
valid = density > 0
plt.figure(figsize=(8, 5))
plt.plot(r_centers[valid], density[valid], label='Shell-Tracked Density')
plt.xscale('log')
plt.yscale('log')
plt.xlim(r_min, r_max)
plt.ylim(np.min(density[valid]) * 0.5, np.max(density[valid]) * 2)
plt.xlabel("Radius [cm]")
plt.ylabel("Density [g/cm³]")
plt.title(f"Wind Density Profile at {t_final / 1e6:.2f} Myr")
plt.grid(True, which='both', ls='--', alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("testdensityplot") 




