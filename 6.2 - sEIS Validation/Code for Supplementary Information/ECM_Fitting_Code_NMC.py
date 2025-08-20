# %%
import time as timer
import matplotlib.pyplot as plt
import numpy as np
import pybamm
import pandas as pd
import scipy as sci
import impedance
from impedance import preprocessing
from impedance.models.circuits import CustomCircuit

# %%
# Set SoC To Analyse
SoC = 100  # Set the state of charge to analyze



# %%
# Import Data

# MATLAB sEIS Data
frequencies_MATLAB, Z_MATLAB = preprocessing.readCSV('Matlab_NMC_EIS_'+str(SoC)+'_SoC.csv')

# PyBaMM sEIS Data
frequencies_PyBaMM, Z_PyBaMM = preprocessing.readCSV('PyBaMM_NMC_EIS_'+str(SoC)+'_SoC.csv')


# %%
# Fit ECM to sEIS Data
# Fit for MATLAB Data
circuit = 'R0-p(R1,C1)-p(R2,CPE1)-p(R3-Wo1,C2)'
initial_guess = [.0002, .004, 2.48, .001,  25000,0.75,0.001,.05, 100,25000]
circuit = CustomCircuit(circuit, initial_guess=initial_guess)
circuit.fit(frequencies_MATLAB, Z_MATLAB)
Z_fit_ECM = circuit.predict(frequencies_MATLAB)


# %%
# Calaculate Initial Fit Error
Z_error = (Z_MATLAB-Z_fit_ECM)/Z_MATLAB


Z_error_real = np.mean(np.abs(np.real(Z_error)))
Z_error_imag = np.mean(np.abs(np.imag(Z_error)))

Z_error_abs = np.mean([Z_error_real, Z_error_imag])

print(f"Initial ECM Fit Error for MATLAB sEIS Data at SoC={SoC}%: {Z_error_abs:.4f}")

# %%
# Limit ECM Fit to frequency range 1mHz<f< 100Hz

k = 10
# Fit for MATLAB Data
circuit = 'R0-p(R1,C1)-p(R2,CPE1)-p(R3-Wo1,C2)'
initial_guess = [.0002, .004, 2.48, .001,  25000,0.75,0.001,.05, 100,25000]
circuit = CustomCircuit(circuit, initial_guess=initial_guess)
circuit.fit(frequencies_MATLAB[:-k], Z_MATLAB[:-k])
Z_fit_ECM_2 = circuit.predict(frequencies_MATLAB[:-k])

# Calaculate Initial Fit Error
Z_error = (Z_MATLAB[:-k]-Z_fit_ECM_2)/Z_MATLAB[:-k]


Z_error_real = np.mean(np.abs(np.real(Z_error)))
Z_error_imag = np.mean(np.abs(np.imag(Z_error)))

Z_error_abs = np.mean([Z_error_real, Z_error_imag])

print(f"Improved ECM Fit Error for MATLAB sEIS Data at SoC={SoC}%: {Z_error_abs:.4f}")


# %%
# Plot Fit Comparison
if SoC == 100:
    plt.figure(1)
    plt.plot(Z_MATLAB.real, -Z_MATLAB.imag, marker="o", linestyle="None", label='MATLAB SoC=100%')
    plt.plot(Z_fit_ECM.real, -Z_fit_ECM.imag, "r--", label='ECM Fit for MATLAB SoC=100%')
    plt.xlabel(r"$Z_\mathrm{Re}$ [Ohm]")
    plt.ylabel(r"$-Z_\mathrm{Im}$ [Ohm]")
    plt.grid()
    plt.legend()
    plt.savefig('ECM_Fit_Example_Plot_Error_NMC.png', dpi=300, bbox_inches='tight')
    plt.show()
    

    plt.figure(2)
    plt.plot(Z_MATLAB[:-k].real, -Z_MATLAB[:-k].imag, marker="o", linestyle="None", label='MATLAB SoC=100%')
    plt.plot(Z_fit_ECM_2.real, -Z_fit_ECM_2.imag, "r--", label='ECM Fit for MATLAB SoC=100%')
    plt.xlabel(r"$Z_\mathrm{Re}$ [Ohm]")
    plt.ylabel(r"$-Z_\mathrm{Im}$ [Ohm]")
    plt.grid()
    plt.legend()
    plt.savefig('ECM_Fit_Example_Plot_NMC.png', dpi=300, bbox_inches='tight')
    plt.show()




# %%
# Fit for PyBaMM Data
circuit = 'R0-p(R1,C1)-p(R2,CPE1)-p(R3-Wo1,C2)'
initial_guess = [.0002, .004, 2.48, .001,  25000,0.75,0.001,.05, 100,25000]
circuit = CustomCircuit(circuit, initial_guess=initial_guess)
circuit.fit(frequencies_PyBaMM[:-k], Z_PyBaMM[:-k])
Z_fit_ECM_PyBaMM = circuit.predict(frequencies_PyBaMM[:-k])

# %%
# Calculate ECM Fit Error Between MATLAB and PyBaMM Data
Z_error = (Z_fit_ECM_PyBaMM-Z_fit_ECM_2)/Z_fit_ECM_PyBaMM

Z_error_real = np.mean(np.abs(np.real(Z_error)))
Z_error_imag = np.mean(np.abs(np.imag(Z_error)))

Z_error_abs = np.mean([Z_error_real, Z_error_imag])

print(f"ECM Fit Error between MATLAB and PyBaMM sEIS Data at SoC={SoC}%: {Z_error_abs:.4f}")


# %%
# Plot ECM Fit Comparison Between MATLAB and PyBaMM Data
plt.figure(3)
plt.plot(Z_fit_ECM_2.real, -Z_fit_ECM_2.imag, "b", label='ECM Fit MATLAB SoC=100%')
plt.plot(Z_fit_ECM_PyBaMM.real, -Z_fit_ECM_PyBaMM.imag, "r--", label='ECM Fit for PyBaMM SoC=100%')
plt.xlabel(r"$Z_\mathrm{Re}$ [Ohm]")    
plt.ylabel(r"$-Z_\mathrm{Im}$ [Ohm]")
plt.title(f'ECM Comparison for SoC = {SoC}%')   
plt.grid()
plt.legend()        
plt.savefig('ECM_Comp_SoC_'+str(SoC)+'NMC.png', dpi=300, bbox_inches='tight')
plt.show()  
# %%
