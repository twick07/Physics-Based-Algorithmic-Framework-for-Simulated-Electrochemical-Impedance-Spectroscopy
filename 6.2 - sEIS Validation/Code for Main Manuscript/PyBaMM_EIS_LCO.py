# %%
# Script to conduct sEIS using PyBaMM
import time as timer
import matplotlib.pyplot as plt
import numpy as np
import pybamm
import pandas as pd
from scipy.fft import fft
from pybamm import exp, constants, sqrt, max

# Import Battery Parameters
param =1
n_parameters = 31
F_const = 96487
V_cut_off_low = 2.8
V_cut_off_high = 4.2
SoC = 100
import Battery_Parameters_MDN as param_file

Model_Parameters = param_file.Battery_Parameters_P2D(param,n_parameters,SoC)

# Set Functions
# OCP Functions
def ocp_neg(sto):
    #import numpy as np
    
    OCP_ref_neg = (
        0.194
        + 1.5 * np.exp(-120.0 * sto)
        + 0.0351 * np.tanh((sto - 0.286) / 0.083)
        - 0.0045 * np.tanh((sto - 0.849) / 0.119)
        - 0.035 * np.tanh((sto - 0.9233) / 0.05)
        - 0.0147 * np.tanh((sto - 0.5) / 0.034)
        - 0.102 * np.tanh((sto - 0.194) / 0.142)
        - 0.022 * np.tanh((sto - 0.9) / 0.0164)
        - 0.011 * np.tanh((sto - 0.124) / 0.0226)
        + 0.0155 * np.tanh((sto - 0.105) / 0.029))
    
    return OCP_ref_neg

def ocp_pos(sto):
    
    OCP_ref_pos = (
        2.16216
        + 0.07645 * np.tanh(30.834 - 54.4806 * sto)
        + 2.1581 * np.tanh(52.294 - 50.294 * sto)
        - 0.14169 * np.tanh(11.0923 - 19.8543 * sto)
        + 0.2051 * np.tanh(1.4684 - 5.4888 * sto)
        + 0.2531 * np.tanh((-sto + 0.56478) / 0.1316)
        - 0.02167 * np.tanh((sto - 0.525) / 0.006)
    )
    return OCP_ref_pos

#Exchange Current Density Functions
def negative_electrode_exchange_current_density(c_e, c_s_surf, c_s_max,T):
    

    #k_ref = Model_Parameters[15]

    # multiply by Faraday's constant to get correct units
    m_ref = (
        2 * 10 ** (-5)
    )  # (A/m2)(m3/mol)**1.5 - includes ref concentrations

    return m_ref *  c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5

def positive_electrode_exchange_current_density(c_e, c_s_surf, c_s_max,T):
    

   

    # multiply by Faraday's constant to get correct units
    m_ref = ( 6 * 10 ** (-7))  # (A/m2)(m3/mol)**1.5 - includes ref concentrations

    return m_ref *  c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5

#Conductivity and Diffusivity Functions
def kappa_function(c_e,T):
    
    kappa = (
        0.0911
        + 1.9101 * (c_e / 1000)
        - 1.052 * (c_e / 1000) ** 2
        + 0.1554 * (c_e / 1000) ** 3
    )

    return kappa



#Set Parameter Values
Q_pos_initial =  (1/3600)*F_const*Model_Parameters[16]* Model_Parameters[19]*Model_Parameters[17]*Model_Parameters[20]*np.abs(1.0 - Model_Parameters[17])
Q_neg_initial = (1/3600)*F_const*Model_Parameters[7]*Model_Parameters[10]*Model_Parameters[8]*Model_Parameters[11]*np.abs(Model_Parameters[8] - 0)
max_cap = np.max([Q_pos_initial,Q_neg_initial])

param_dict = {'Thermodynamic factor': Model_Parameters[28],
              'Ambient temperature [K]': 298.15,
              'Current function [A]': max_cap,
              "Electrode height [m]": 1,
             "Electrode width [m]": 1,
             "Negative electrode OCP [V]": ocp_neg, 
             "Positive electrode OCP [V]": ocp_pos,
             'Positive electrode Bruggeman coefficient (electrode)': Model_Parameters[1],
             'Positive electrode Bruggeman coefficient (electrolyte)': Model_Parameters[1],
             'Negative electrode Bruggeman coefficient (electrode)': Model_Parameters[1],
             'Negative electrode Bruggeman coefficient (electrolyte)': Model_Parameters[1],
             'Separator Bruggeman coefficient (electrolyte)': Model_Parameters[1],
             'Separator Bruggeman coefficient (electrode)': Model_Parameters[1],
              'Separator thickness [m]': Model_Parameters[2],
              'Initial concentration in electrolyte [mol.m-3]': Model_Parameters[3],
              'Electrolyte diffusivity [m2.s-1]': Model_Parameters[4],
              'Cation transference number': Model_Parameters[5],
              'Separator porosity': Model_Parameters[6],
              'Negative electrode thickness [m]': Model_Parameters[7],
              'Initial concentration in negative electrode [mol.m-3]': Model_Parameters[8]*Model_Parameters[11],
              'Negative electrode porosity': Model_Parameters[9],
              'Negative electrode active material volume fraction': Model_Parameters[10],
              'Maximum concentration in negative electrode [mol.m-3]': Model_Parameters[11],
              'Negative particle radius [m]': Model_Parameters[12],
              'Negative electrode diffusivity [m2.s-1]': Model_Parameters[13],
              'Negative electrode conductivity [S.m-1]': Model_Parameters[14],
              'Negative electrode charge transfer coefficient': 0.5,
              'Negative electrode exchange-current density [A.m-2]': negative_electrode_exchange_current_density, 
              'Negative electrode double-layer capacity [F.m-2]': Model_Parameters[26],
              'Positive electrode thickness [m]': Model_Parameters[16],
              'Initial concentration in positive electrode [mol.m-3]': Model_Parameters[17]*Model_Parameters[20],
              'Positive electrode porosity': Model_Parameters[18],
              'Positive electrode active material volume fraction': Model_Parameters[19],
              'Maximum concentration in positive electrode [mol.m-3]': Model_Parameters[20],
              'Positive particle radius [m]': Model_Parameters[21],
              'Positive electrode diffusivity [m2.s-1]': Model_Parameters[22],
              'Positive electrode conductivity [S.m-1]': Model_Parameters[23],
              'Electrolyte conductivity [S.m-1]': kappa_function,
              #'Electrolyte conductivity [S.m-1]':Model_Parameters[25],
              'Positive electrode double-layer capacity [F.m-2]': Model_Parameters[27],
              'Positive electrode charge transfer coefficient': 0.5,
              'Positive electrode exchange-current density [A.m-2]': positive_electrode_exchange_current_density, 
              'Initial temperature [K]': 298.15,
              'Lower voltage cut-off [V]': V_cut_off_low,
              'Reference temperature [K]': 298.15,
              'Upper voltage cut-off [V]': V_cut_off_high,
              'Maximum voltage [V]': V_cut_off_high,
              'Typical electrolyte concentration [mol.m-3]': Model_Parameters[3],
              'Negative electrode electrons in reaction': 1.0,
              'Positive electrode electrons in reaction': 1.0,
              'Negative electrode OCP entropic change [V.K-1]': 0.0,
              'Positive electrode OCP entropic change [V.K-1]': 0.0,
              'Nominal cell capacity [A.h]': 5.0,
              'Number of cells connected in series to make a battery': 1.0,
              'Number of electrodes connected in parallel to make a cell': 1.0,
              'Typical current [A]': max_cap,
              'Negative particle distribution in x': 1.0, 
              'Positive particle distribution in x': 1.0,
              #'Cell capacity [A.h]': 1.0 
              'Cell capacity [A.h]': max_cap*3600, 
              'Open-circuit voltage at 100% SOC [V]': ocp_pos(Model_Parameters[17]) - ocp_neg(Model_Parameters[8]),
              'Open-circuit voltage at 0% SOC [V]': V_cut_off_low          
              }

# %%
# Simulate a 1C Discharge
# Set Model Type
var_pts = {
"x_n": 20,  # negative electrode
"x_s": 20,  # separator
"x_p": 20,  # positive electrode
"r_n": 10,  # negative particle
"r_p": 10,  # positive particle
}

model =pybamm.lithium_ion.DFN({"surface form": "differential","thermal": "isothermal","SEI":"none","lithium plating":"none"},var_pts)
parameter_values = pybamm.ParameterValues(param_dict)


# %%
#Simulate and EIS Curve
#I_app = 5e-3
safe_solver = pybamm.CasadiSolver(atol=1e-12, rtol=1e-12, mode="safe")
I_app = 0.2
number_of_periods = 10
samples_per_period = 56
No_of_frequencies = 40
No_of_output_cycles = 10
time_vector = np.zeros([No_of_frequencies,samples_per_period*No_of_output_cycles])
current_vector = np.zeros([No_of_frequencies,samples_per_period*No_of_output_cycles])
voltage_vector = np.zeros([No_of_frequencies,samples_per_period*No_of_output_cycles])
frequencies = np.logspace(-4, 3, No_of_frequencies)

def current_function(t):
    return I_app * pybamm.sin(2 * np.pi * pybamm.InputParameter("Frequency [Hz]") * t)
    #return I_app * pybamm.cos(2 * np.pi * pybamm.InputParameter("Frequency [Hz]") * t)


parameter_values["Current function [A]"] = current_function




sim = pybamm.Simulation(model,parameter_values=parameter_values,solver=safe_solver)

impedances_time = []
idx_1 = 0
for frequency in frequencies:
    # Solve
    period = 1 / frequency
    dt = period / samples_per_period
    t_eval = np.array(range(0, samples_per_period * number_of_periods)) * dt
    sol = sim.solve(t_eval, inputs={"Frequency [Hz]": frequency})
    # Extract final two periods of the solution
    #time = sol["Time [s]"].entries[-5 * samples_per_period :]
    #current = sol["Current [A]"].entries[-5 * samples_per_period :]
    #voltage = sol["Voltage [V]"].entries[-5 * samples_per_period :]
    time = sol["Time [s]"].entries[-No_of_output_cycles * samples_per_period :]
    current = sol["Current [A]"].entries[-No_of_output_cycles * samples_per_period :]
    voltage = sol["Voltage [V]"].entries[-No_of_output_cycles * samples_per_period :]

    # Store the current
    time_vector[idx_1,:] = time
    current_vector[idx_1,:] = current
    voltage_vector[idx_1,:] = voltage


    # FFT
    current_fft = fft(current)
    voltage_fft = fft(voltage)
    # Get index of first harmonic
    idx = np.argmax(np.abs(current_fft))
    impedance = -voltage_fft[idx] / current_fft[idx]
    impedances_time.append(impedance)
    idx_1 += 1



# %%

# Plotting
data = np.array(impedances_time)
plt.figure(1)
plt.plot(data.real, -data.imag, marker="o", linestyle="None")
plt.xlabel(r"$Z_\mathrm{Re}$ [Ohm]")
plt.ylabel(r"$-Z_\mathrm{Im}$ [Ohm]")
plt.show()



# %%
# Export the data
output_data_current = {}
output_data_voltage = {}
output_data_time = {}



output_data_impedance = {
    "Frequencies":frequencies,
    "Real Values":np.real(data),
    "Imag Values":np.imag(data)
}
# Export the Impedences
for idx in range(len(frequencies)):
    output_data_current["F = "+str(frequencies[idx])] = current_vector[idx,:]
    output_data_voltage["F = "+str(frequencies[idx])] = voltage_vector[idx,:]
    output_data_time["F = "+str(frequencies[idx])] = time_vector[idx,:]

# %%
# Export Data

SoC_text = str(SoC)


#load output data into a DataFrame object:
df = pd.DataFrame(output_data_time)

# writing output data frame to a CSV file
df.to_csv('LCO_EIS_time_vectors_'+SoC_text+'.csv')

#load output data into a DataFrame object:
df = pd.DataFrame(output_data_current)

# writing output data frame to a CSV file
df.to_csv('LCO_EIS_current_vectors_'+SoC_text+'.csv')

#load output data into a DataFrame object:
df = pd.DataFrame(output_data_voltage)

# writing output data frame to a CSV file
df.to_csv('LCO_EIS_voltage_vectors_'+SoC_text+'.csv')

#load output data into a DataFrame object:
df = pd.DataFrame(output_data_impedance)

# writing output data frame to a CSV file
df.to_csv('LCO_EIS_Impedance_'+SoC_text+'.csv')


# %%
