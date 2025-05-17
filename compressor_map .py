
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Defining given constants

T01 = 293.0           # Stagnation temperature at station 1 (in K)
R = 287.0             # Gas constant for air (J/(kg*K))
gamma = 1.4
mu = 1.83e-5          # Dynamic viscosity (Pa.s)
Re_target = 3.0e6     # Target Reynolds number
M_target = 0.55       # Mach number at station 4
l_prototype = 0.2     # Prototype characteristic length (m)
L = 1.5               # Total length of the diffuser (m)
D3 = 0.6               # Diameter of diffues at station 3

# Loading compressor _operation data

file_path = "Compressor_operation_map.txt"
df = pd.read_csv(file_path, sep="\t")
df.columns = ['mdot_ref', 'T0_ratio', 'p0_ratio']

'''
#If colab is being used
file_path = "/content/Compressor_operation_map.txt"
df = pd.read_csv(file_path, sep="\t")
df.columns = ['mdot_ref', 'T0_ratio', 'p0_ratio']
'''
# Defining Variables

scale_result = -1.0         # Will store the maximum (feasible) scale
p01_result = None           # Inlet stagnation pressure corresponding to the best scale
p0_ratio_result = None      # Store p0_ratio (from the compressor map) for the best point
mdot_result = None          # mdot_ref for the best point

# Iterating through each row

for index, row in df.iterrows():
    mdot_ref = row['mdot_ref']
    T0_ratio = row['T0_ratio']
    p0_ratio = row['p0_ratio']

    # Stagnation temperature at station 2:

    T02 = T0_ratio * T01  # [K]

    # Static temperature at station 4 using isentropic relation:
    T4 = T02 / (1.0 + 0.5 * (gamma - 1.0) * (M_target ** 2))

    # Computing local speed of sound and flow velocity at station 4:
    a4 = np.sqrt(gamma * R * T4)    # [in m/s]
    V4 = M_target * a4              # [in m/s]

    # Compute stagnation pressure at station 4 using the isentropic relation:
    # p04 = p03 and p03 = 0.985 * p02
    # Static pressure p4 = p04 / (1.0 + 0.5 * (gamma - 1.0) * (M_target ** 2)) ** (gamma / (gamma - 1.0)):

    factor = (1.0 + 0.5 * (gamma - 1.0) * (M_target ** 2)) ** (gamma / (gamma - 1.0))

    #  Calculating l_model: l_model = sqrt((mdot * R * T4) / (0.985 * p4 * V4 * pi  * 101325)) (using p_ref = 101325 Pa):
    # mdot = mdot_ref * p01/101325(in Pa) *(sqrt(T01/293))

    l_model = np.sqrt((mdot_ref * R * T4) / (0.985 * p0_ratio / factor * V4 * np.pi  * 101325))

    # Finding scale

    scale = l_model/ l_prototype  # [in m]

    # Calculating the density at station 4 using the Reynolds number relation:
    # rho4 = (Re_target * mu) / (V4 * l_model)

    rho4 = (Re_target * mu) / (V4 * l_model)  # [in kg/m^3]

    # Calculate static pressure at station 4 via ideal gas law:

    p4_pa = rho4 * R * T4   # [in Pa]
    p4 = p4_pa / 1000.0     # [in kPa]

    p04 = p4 * factor       # [in kPa]

    # Compressor exit pressure: p02 = p04 / 0.985:

    p02 = p04 / 0.985       # [in kPa]

    # Inlet stagnation pressure: p01 = p02 / p0_ratio:

    p01 = p02 / p0_ratio    # [in kPa]

    # Select scale only if scale <= 1 and p02 < 250 kPa:

    if scale <= 1 and p02 <= 250:
        if scale > scale_result:
            scale_result = scale
            p01_result = p01
            p02_result = p02
            l_model_result = l_model
            p0_ratio_result = p0_ratio
            mdot_result = mdot_ref

# Calculating the value of theta = arctan((D4 - D3) /2 / L)

theta_rad = np.arctan(((2* l_model_result - D3) / 2) / L)

theta_deg = np.degrees(theta_rad)

# Optimal Operating Point Results:

if scale_result > 0:
    print("\n--- Optimal Operating Point ---")
    print("Maximum Scale: {:.4f}".format(scale_result))
    print("Inlet Stagnation Pressure, p01: {:.2f} kPa".format(p01_result))
    print("p0_ratio (p02/p01) from Compressor Map: {:.4f}".format(p0_ratio_result))
    print("Compressor exit Stagnation Pressure, p02: {:.2f} kPa".format(p02_result))
    print("m_dot_ref: {:.4f} kg/s".format(mdot_result))
    print("Theta: {:.4f} degrees".format(theta_deg))
else:
    print("No feasible operating point found that meets the constraints (scale <= 1 and p02 < 250 kPa).")


# Plot p0_ratio vs. m_dot_ref and Mark the Optimal Point:

plt.figure(figsize=(10,6))
plt.plot(df['mdot_ref'], df['p0_ratio'], '-o', label='Compressor Map')
if scale_result > 0:
    plt.scatter(mdot_result, p0_ratio_result, color='red', marker='D', s=120,
                label='Optimal Operating Point', zorder=10)
plt.xlabel('m_dot_ref (kg/s)')
plt.ylabel('p0_ratio (p02/p01)')
plt.title('Compressor Map & Optimal Operating Point')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
