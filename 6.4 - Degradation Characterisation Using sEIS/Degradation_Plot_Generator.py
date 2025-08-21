# %%
import time as timer
import matplotlib.pyplot as plt
import numpy as np
import pybamm
import pandas as pd
import scipy as sci
import impedance
from impedance import preprocessing


# %%
# Select Degradation Condition
SEI_con = 1 # SEI_con = 1 for KL SEI Growth, SEI_con = 2 for DL SEI Growth
# %%
# Import Data
if SEI_con == 1:
    frequencies = pd.read_csv('Matlab_EIS_Frequencies_KL.csv')
    Z_real = pd.read_csv('Matlab_EIS_Real_Impedance_KL.csv')
    Z_imag = pd.read_csv('Matlab_EIS_Imag_Impedance_KL.csv')

if SEI_con == 2:
    frequencies = pd.read_csv('Matlab_EIS_Frequencies_DL.csv')
    Z_real = pd.read_csv('Matlab_EIS_Real_Impedance_DL.csv')
    Z_imag = pd.read_csv('Matlab_EIS_Imag_Impedance_DL.csv')


# %%
# Preprocess Data
m = 7   
N_spectra = int(frequencies.shape[1])
N_frequencies = int(frequencies.shape[0])

frequency_vectors = frequencies.to_numpy()
Z_vectors_real = Z_real.to_numpy()
Z_vectors_imag = Z_imag.to_numpy()
Z_complex = Z_vectors_real - Z_vectors_imag*1j

frequency_vectors = frequency_vectors[:-m,:]
Z_complex = Z_complex[:-m,:]
Z_fit_vectors = np.zeros(Z_complex.shape)
Z_fit_vectors = Z_fit_vectors.astype(np.complex128)
# %%
# Define the impedance model
# Fit Circuit to  impedance data
from impedance.models.circuits import CustomCircuit

#ECM_Parameters = np.zeros([9,N_spectra])
ECM_Parameters = np.zeros([9,N_spectra])

for k in range(0,N_spectra):
    circuit = 'R0-p(R1,C1)-p(R2,C2)-p(R3-Wo1,C3)'
    initial_guess = [.0002, .004, 2.48, .001,  25000,0.001,.05, 100,25000]
    circuit = CustomCircuit(circuit, initial_guess=initial_guess)
    circuit.fit(frequency_vectors[:,k], Z_complex[:,k])
    ECM_Parameters[:,k] = circuit.parameters_
    Z_fit_vectors[:,k] = circuit.predict(frequency_vectors[:,k])
    
# %%
# Import Cap lost and SoH data from matlab
if SEI_con == 1:
    SoH = pd.read_csv('50_cycle_KL_ageing_SoH_lost.csv')
    SoH = SoH.to_numpy()
    Cap_lost = pd.read_csv('50_cycle_KL_ageing_cap_lost.csv')
    Cap_lost = Cap_lost.to_numpy()  
if SEI_con == 2:
    SoH = pd.read_csv('50_cycle_DL_ageing_SoH_lost.csv')
    SoH = SoH.to_numpy()    
    Cap_lost = pd.read_csv('50_cycle_DL_ageing_cap_lost.csv')
    Cap_lost = Cap_lost.to_numpy()  


# %%
# Plotting the results
if SEI_con == 1:
    SEI_text = "KL"
if SEI_con == 2:
    SEI_text = "DL"

# Plot sEIS Evolution
plt.figure(figsize=(10, 6))
cycles_to_plot = [1, int(np.round(N_spectra/2)) - 1, N_spectra - 1]
colors = ['r', 'b', 'k']

for k, color in zip(cycles_to_plot, colors):
    plt.plot(Z_fit_vectors[:,k].real, 
             -Z_fit_vectors[:,k].imag, 
             color, 
             label=f'ECM Spectra for Cycle {k}')

plt.xlabel(r"$Z_\mathrm{Re}$ [Ohm]")
plt.ylabel(r"$-Z_\mathrm{Im}$ [Ohm]")
plt.grid(True)
plt.legend()
plt.title("Fitted Impedance Spectra")
plt.savefig('50_cycles_sEIS_plots_'+SEI_text+'.png', dpi=300, bbox_inches='tight')
plt.show()

# Plot R_0 vs Capacity Loss

Q_lost_percent = Cap_lost[:,1]

fig, ax1 = plt.subplots()
ax1.set_xlabel('Cycles [-]')
ax1.set_ylabel(r"$R_\mathrm{0}$ [Ohm]", color='tab:red')
ax1.plot(np.linspace(1,ECM_Parameters.shape[1],ECM_Parameters.shape[1]), ECM_Parameters[0,:], 'r',label='ECM Parameter $R_0$')
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.grid()
plt.legend(loc=0)
ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis
ax2.set_ylabel('Lost Capacity [%]', color='tab:blue')  # we already handled the x-label with ax1
ax2.plot(np.linspace(1,ECM_Parameters.shape[1],ECM_Parameters.shape[1]), Q_lost_percent, 'b--',label='Simulations Lost Capacity')
ax2.tick_params(axis='y', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
ax2.grid()
plt.legend(loc=1)
plt.savefig(SEI_text+'_R_0_Fit_Cap.png', dpi=300, bbox_inches='tight')
plt.show()

# Plot R_1 vs Capacity Loss
fig, ax1 = plt.subplots()
ax1.set_xlabel('Cycles [-]')
ax1.set_ylabel(r"$R_\mathrm{1}$ [Ohm]", color='tab:red')
ax1.plot(np.linspace(1,ECM_Parameters.shape[1],ECM_Parameters.shape[1]), ECM_Parameters[1,:], 'r',label='ECM Parameter $R_1$')
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.grid()
plt.legend(loc=0)
ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis
ax2.set_ylabel('Lost Capacity [%]', color='tab:blue')  # we already handled the x-label with ax1
#ax2.plot(np.linspace(1,ECM_Parameters.shape[1],ECM_Parameters.shape[1]), Cap_lost[:,1], 'b--',label='Simulations Lost Capacity')
ax2.plot(np.linspace(1,ECM_Parameters.shape[1],ECM_Parameters.shape[1]), Q_lost_percent, 'b--',label='Simulations Lost Capacity')
ax2.tick_params(axis='y', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
ax2.grid()
plt.legend(loc=1)
plt.savefig(SEI_text+'_R_1_Fit_Cap.png', dpi=300, bbox_inches='tight')
plt.show()

# Plot C_1 vs SoH
fig, ax1 = plt.subplots()
ax1.set_xlabel('Cycles [-]')
ax1.set_ylabel(r"$C_\mathrm{1}$ [F]", color='tab:red')
ax1.plot(np.linspace(1,ECM_Parameters.shape[1],ECM_Parameters.shape[1]), ECM_Parameters[2,:], 'r',label='ECM Parameter $C_1$')
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.grid()
plt.legend(loc=2)
ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis
ax2.set_ylabel('SoH [%]', color='tab:blue')  # we already handled the x-label with ax1
ax2.plot(np.linspace(1,ECM_Parameters.shape[1],ECM_Parameters.shape[1]), SoH[:,1], 'b--',label='Simulations State of Health')
ax2.tick_params(axis='y', labelcolor='tab:blue')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
ax2.grid()
plt.legend(loc=0)
plt.savefig(SEI_text+'_C_1_Fit_SoH.png', dpi=300, bbox_inches='tight')
plt.show()
# %%
