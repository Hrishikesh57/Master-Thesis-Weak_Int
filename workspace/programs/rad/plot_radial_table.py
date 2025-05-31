

import numpy as np
import matplotlib.pyplot as plt

# Read data from the file, skipping the header line
try:
    
    data = np.loadtxt('92_polynomial_table_2s.txt', skiprows=2)
    #data = np.loadtxt('radial_function.dat', skiprows=1)
    print("Data loaded successfully.")
except OSError as e:
    print(f"Error: {e}")
    exit(1)

# Verify data shape
print('Data shape:', data.shape)

# Ensure data has at least 3 columns
if data.shape[1] < 3:
    print("Error: Data does not have enough columns.")
    exit(1)

x = data[:, 0]
RadialG = data[:, 1]
RadialF = data[:, 2]

# Check the ranges of RadialG and RadialF
print('RadialG(x): min =', np.min(RadialG), ', max =', np.max(RadialG))
print('RadialF(x): min =', np.min(RadialF), ', max =', np.max(RadialF))

# Plot RadialG(x)
plt.figure(figsize=(10, 6))
plt.plot(x, RadialG, label='RadialG(x)', color='blue')
plt.xlabel('r (a.u.)')
plt.ylabel('RadialG(x)')
plt.title('Plot of RadialG(x) on Linear Scale')
plt.legend()
plt.grid(True)
plt.show()


# Plot RadialF(x)
plt.figure(figsize=(10, 6))
plt.plot(x, RadialF, label='RadialF(x)', color='green')
plt.xlabel('r (relativisitic units)')
plt.ylabel('RadialF(x)')
plt.title('Plot of RadialF(x) on Linear Scale')
plt.legend()
plt.grid(True)
plt.show()


# Plot both RadialG(x) and RadialF(x) on linear scale
plt.figure(figsize=(10, 6))
plt.plot(x, RadialG, label='RadialG(x)', color='blue')
plt.plot(x, RadialF, label='RadialF(x)', color='green')
plt.xlabel('r (a.u.)')
plt.ylabel('Radial Functions')
plt.title('Radial Functions RadialG(x) and RadialF(x) on Linear Scale')
plt.legend()
plt.grid(True)
plt.show()

