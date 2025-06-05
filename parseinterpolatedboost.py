import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

filename = "# m =  9 MSun =  1.78963e+34 g.dat"

with open(filename) as f:
    header = f.readline().strip()
    header_names = [col.strip() for col in header.split()]

wanted = ["mZAMS", "time(02)", "mass(03)", "dotM(04)", "vwnd(05)"]
usecols = [i for i, name in enumerate(header_names) if name in wanted]
names = [header_names[i] for i in usecols]

print("Using names:", repr(names))  # Helpful debug print

data = np.genfromtxt(
    filename,
    skip_header=1,
    usecols=usecols,
    names=names,
    autostrip=True
)

# Extract the data
time_myr = data["time02"] / 1e6  # Convert time from seconds to Myr
mass = data["mass03"]
Mdot = data["dotM04"]
vwind = data["vwnd05"]


fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(8, 6), sharex=True)

# Plot vwind
ax1.plot(time_myr, vwind, color="steelblue")
ax1.set_ylabel("Wind velocity (km/s)")
ax1.set_title("Wind velocity and Mass Loss Rate over Time")
ax1.grid(True)

# Plot mdot
ax2.plot(time_myr, Mdot, color="darkred")
ax2.set_xlabel("Time (Myr)")
ax2.set_ylabel("Mass loss rate (M☉/yr)")
ax2.grid(True)

plt.tight_layout()
plt.show()