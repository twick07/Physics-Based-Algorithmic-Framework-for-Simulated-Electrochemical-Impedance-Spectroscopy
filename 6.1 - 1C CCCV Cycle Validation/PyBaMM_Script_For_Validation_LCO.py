# %%
# Script to Validate with a PyBaMM drive cycle
import time as timer
import matplotlib.pyplot as plt
import numpy as np
import pybamm
import pandas as pd
import scipy 
from scipy import interpolate
from scipy.fft import fft
from pybamm import exp, constants, sqrt, max

# %%
# Import MATLAB Drive Cycle Data
data = scipy.io.loadmat("Cycler_Output_LCO.mat")
time_drive_cycle = data.get("time_vector_total")
current_drive_cycle = data.get('I_app_total')
voltage_drive_cycle = data.get('V_cell_total')
# %%
# Start With PyBaMM Built In Model
# Define All Required Functions
def Applied_Current_Function(t):
   
    I = pybamm.Interpolant(x=time_drive_cycle.flatten(),y=current_drive_cycle.flatten(),children=t,interpolator="cubic",extrapolate=False)

   
    return(I)


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
def kappa_function_2(c_e,T):
    
    kappa = (
        0.0911
        + 1.9101 * (c_e / 1000)
        - 1.052 * (c_e / 1000) ** 2
        + 0.1554 * (c_e / 1000) ** 3
    )

    return kappa
#Set Parameter Values
param =1 # param = 1 for LCO/C6, param =2 for NMC532/C6

n_parameters = 50

import Battery_Parameters_MDN_2 as param_file

Model_Parameters = param_file.Battery_Parameters_P2D(param,n_parameters)
F_const = 96487  # Faraday Constant [C/mol]
Q_pos_initial =  (1/3600)*F_const*Model_Parameters[16]* Model_Parameters[19]*Model_Parameters[17]*Model_Parameters[20]*np.abs(1.0 - Model_Parameters[17])
Q_neg_initial = (1/3600)*F_const*Model_Parameters[7]*Model_Parameters[10]*Model_Parameters[8]*Model_Parameters[11]*np.abs(Model_Parameters[8] - 0)
max_cap = np.max([Q_pos_initial,Q_neg_initial])
V_cut_off_low = 3.3
V_cut_off_high = 4.2

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
              'Electrolyte conductivity [S.m-1]': kappa_function_2,
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

# Simulate a 1C Discharge
# Set Model Type
var_pts = {
"x_n": 20,  # negative electrode
"x_s": 10,  # separator
"x_p": 20,  # positive electrode
"r_n": 10,  # negative particle
"r_p": 10,  # positive particle
}


model =pybamm.lithium_ion.DFN({"surface form": "differential","thermal": "isothermal","SEI":"none","lithium plating":"none"},var_pts)

drive_cycle_current= np.column_stack([time_drive_cycle.T, current_drive_cycle.T])
Experiment = pybamm.Experiment([pybamm.step.current(drive_cycle_current)])

safe_solver = pybamm.CasadiSolver(atol=1e-6, rtol=1e-6, mode="safe",return_solution_if_failed_early=True,dt_max=None)

parameter_values = pybamm.ParameterValues(param_dict)

experiment =pybamm.Experiment([pybamm.step.current(drive_cycle_current)])


sim = pybamm.Simulation(model, experiment=experiment,parameter_values=parameter_values,solver=safe_solver)
solution_pybamm = sim.solve()

# %%
# Post Processing
t_pybamm = solution_pybamm["Time [s]"]
time_vector_pybamm = t_pybamm.entries
V_pybamm = solution_pybamm["Voltage [V]"]
V_cell_vector_pybamm = V_pybamm.entries
I_app_pybamm = solution_pybamm["Current [A]"]
I_app_pybamm_vector = I_app_pybamm.entries

# %%
# Build Custom P2D-DFN Model
# Define All Required Functions

# Define All Required Functions
def Applied_Current_Function(t):
   
    I = pybamm.Interpolant(x=time_drive_cycle.flatten(),y=current_drive_cycle.flatten(),children=t,interpolator="cubic",extrapolate=False)

   
    return(I)


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


#Conductivity and Diffusivity Functions
def kappa_function(c_e):
    
    kappa = (
        0.0911
        + 1.9101 * (c_e / 1000)
        - 1.052 * (c_e / 1000) ** 2
        + 0.1554 * (c_e / 1000) ** 3
    )

    return kappa

# Define Model

P2D_Model = pybamm.BaseModel()

# Define Model Variable
SEI_tol = -1e-3
SEI_con = 0

# Define Spatial Variables
x_neg = pybamm.SpatialVariable("x_neg", domain = ['Negative Electrode'], coord_sys="cartesian")
x_sep = pybamm.SpatialVariable("x_sep", domain = ['Separator'], coord_sys="cartesian")
x_pos = pybamm.SpatialVariable("x_pos", domain = ['Positive Electrode'], coord_sys="cartesian")
r_neg = pybamm.SpatialVariable("r_neg", domain = ["Negative Particle"],auxiliary_domains={"secondary": "Negative Electrode"}, coord_sys="spherical polar")
r_pos = pybamm.SpatialVariable("r_pos", domain = ["Positive Particle"],auxiliary_domains={"secondary": "Positive Electrode"}, coord_sys="spherical polar")

# Define Dependent Variables
# Potential Variables
psi_s_n = pybamm.Variable("Negative Solid Potential [V]", domain = "Negative Electrode")
psi_s_p = pybamm.Variable("Positive Solid Potential [V]", domain = "Positive Electrode")
psi_e_n = pybamm.Variable("Negative Electrolyte Potential [V]",domain="Negative Electrode")
psi_e_s = pybamm.Variable("Separator Electrolyte Potential [V]",domain="Separator")
psi_e_p = pybamm.Variable("Positive Electrolyte Potential [V]",domain="Positive Electrode")
psi_e = pybamm.concatenation(psi_e_n, psi_e_s, psi_e_p)

#Flux Variables
j_flux_n = pybamm.Variable("Negative Intercalation Flux [mols/m^2s]",domain="Negative Electrode")
j_flux_s = pybamm.Variable("Separator Intercalation Flux [mols/m^2s]",domain="Separator")
j_flux_p = pybamm.Variable("Positive Intercalation Flux [mols/m^2s]",domain="Positive Electrode")
j_flux = pybamm.concatenation(j_flux_n, j_flux_s, j_flux_p)
j_tot_n = pybamm.Variable("Negative Total Flux [mols/m^2s]",domain="Negative Electrode")
j_tot_s = pybamm.Variable("Separator Total Flux [mols/m^2s]",domain="Separator")
j_tot_p = pybamm.Variable("Positive Total Flux [mols/m^2s]",domain="Positive Electrode")
#j_tot = pybamm.concatenation(j_tot_n, j_tot_s, j_tot_p)

# SEI Variables
j_SEI_n = pybamm.Variable("Negative Intercalation Flux [mols/m^2s]",domain="Negative Electrode")
c_SEI_n = pybamm.Variable("Negative SEI Concentration [mols/m^3]",domain="Negative Electrode")
L_SEI_n = pybamm.Variable("Negative SEI Thickness [m]",domain="Negative Electrode")
R_SEI_n = pybamm.Variable("Negative SEI Resistance [Ohm m^2]",domain="Negative Electrode")

#Concentration Variables
c_n = pybamm.Variable("Negative Particle Concentration [mols/m^3]", domain = "Negative Particle",auxiliary_domains={"secondary": "Negative Electrode"})
c_p = pybamm.Variable("Positive Particle Concentration [mols/m^3]", domain = "Positive Particle",auxiliary_domains={"secondary": "Positive Electrode"})
c_e_n = pybamm.Variable("Negative Electrolyte Concentration [mol.m-3]",domain="Negative Electrode")
c_e_s = pybamm.Variable("Separator Electrolyte Concentration [mol.m-3]",domain="Separator")
c_e_p = pybamm.Variable("Positive Electrolyte Concentration [mol.m-3]",domain="Positive Electrode")
c_e = pybamm.concatenation(c_e_n, c_e_s, c_e_p)

# Define Model Parameters
R_n = pybamm.Parameter("Negative Particle Radius [m]")
R_p = pybamm.Parameter("Positive Particle Radius [m]")
eps_s_n = pybamm.Parameter("Negative Particle Solid Volume Fraction [-]")
eps_s_p = pybamm.Parameter("Positive Particle Solid Volume Fraction [-]")
T_amb = pybamm.Parameter("Ambient Temperature [K]")
SoC_init_n = pybamm.Parameter("Negative Particle Initial State of Charge [-]")
SoC_init_p = pybamm.Parameter("Positive Particle Initial State of Charge [-]")
C_max_neg = pybamm.Parameter("Negative Particle Maximum Concentration [mols/m^3]")
C_max_pos = pybamm.Parameter("Positive Particle Maximum Concentration [mols/m^3]")
I_app = pybamm.FunctionParameter("Current function [A]", {"Time [s]": pybamm.t})
D_s_n = pybamm.Parameter("Negative Solid Diffusion Coefficient [m^2/s]")
D_s_p = pybamm.Parameter("Positive Solid Diffusion Coefficient [m^2/s]")
D_e = pybamm.Parameter("Elecrolyte Ionic Diffusion Coefficient [m^2/s]")
L_n = pybamm.Parameter("Length of Negative Electrode [m]")
L_s = pybamm.Parameter("Length of Separator [m]")
L_p = pybamm.Parameter("Length of Positive Electrode [m]")
L_tot = pybamm.Parameter("Total Cell Length [m]")
k_n = pybamm.Parameter("Negative Rate Constant [mols/m^2s]")
k_p = pybamm.Parameter("Positive Rate Constant [mols/m^2s]")
sigma_n = pybamm.Parameter("Negative Electrode Electronic Conductivity [S/m]")
sigma_p = pybamm.Parameter("Positive Electrode Electronic Conductivity [S/m]")
F = pybamm.Parameter("Faraday Constant [cols/mol]")
R_const = pybamm.Parameter("Ideal Gas Constant [J/(K*mol)]")
Brugg = pybamm.Parameter("Bruggeman Coefficient [-]")
t_li = pybamm.Parameter("Transport Coefficient [-]")
c_e_init = pybamm.Parameter("Initial Electrolyte Concentration [mols/m^3]")
V_cell_low_cut_off = pybamm.Parameter("Low Voltage Cutt-off [V]")
OCP_n_init = pybamm.FunctionParameter("Initial Negative Electrode OCP [V]", {"Negative Particle Initial State of Charge [-]": SoC_init_n})
OCP_p_init = pybamm.FunctionParameter("Initial Positive Electrode OCP [V]", {"Positive Particle Initial State of Charge [-]": SoC_init_p})

# SEI Parameters
## Anode SEI Specific Parameters
kappa_SEI = pybamm.Parameter("SEI Conductivity [S/m]")
M_SEI = pybamm.Parameter("SEI Molar Mass [kg/mol]")
rho_SEI = pybamm.Parameter("SEI Layer Density [kg/m^3]")
i_0_SEI = pybamm.Parameter("SEI Exchange Current [A/m^2]")
OCP_SEI = pybamm.Parameter("SEI OCP [V]")
alpha_SEI = pybamm.Parameter("Charge transfer coefficient SEI")
c_SEI_init = pybamm.Parameter("Initial Negative SEI Concentration [mols/m^3]")
c_sol_0 = pybamm.Parameter("Concentration of Solvent in bulk electrolyte [mol/m^3]")
D_SEI = pybamm.Parameter("EC diffusion coefficient [m^2/s]")
k_SEI = pybamm.Parameter("SEI Reaction Rate constant [m/s]")

# Porosity
# Primary broadcasts are used to broadcast scalar quantities across a domain
# into a vector of the right shape, for multiplying with other vectors
#eps_e_n = pybamm.PrimaryBroadcast(pybamm.Parameter("Negative Electrode Porosity"), "Negative Electrode")
eps_e_n =pybamm.Variable("Negative Electrode Porosity Variable [-]",domain="Negative Electrode")
eps_e_n_init = pybamm.Parameter("Negative Electrode Porosity")
deps_e_n_dt = pybamm.Variable("Negative Electrode Porosity Rate of Change [s-1]",domain="Negative Electrode")

eps_e_s = pybamm.PrimaryBroadcast(pybamm.Parameter("Separator Porosity"), "Separator")
eps_e_p = pybamm.PrimaryBroadcast(pybamm.Parameter("Positive Electrode Porosity"), "Positive Electrode")
eps_e = pybamm.concatenation(eps_e_n, eps_e_s, eps_e_p)
# Electrolyte Conductivity
kappa_n = pybamm.FunctionParameter("Negative Electrolyte Conductivity [S/m]", {"Negative Electrolyte Concentration [mol.m-3]": c_e_n})
kappa_s = pybamm.FunctionParameter("Separator Electrolyte Conductivity [S/m]", {"Separator Electrolyte Concentration [mol.m-3]": c_e_s})
kappa_p = pybamm.FunctionParameter("Positive Electrolyte Conductivity [S/m]", {"Positive Electrolyte Concentration [mol.m-3]": c_e_p})
kappa = pybamm.concatenation(kappa_n,kappa_s,kappa_p)
# Define Intermediate Parameters and Corrected Parameters

stoic_n = pybamm.surf(c_n)/C_max_neg
stoic_p = pybamm.surf(c_p)/C_max_pos
OCP_p = pybamm.FunctionParameter("Positive Electrode OCP [V]", {"Positive Stoichiometry": stoic_p})
OCP_n = pybamm.FunctionParameter("Negative Electrode OCP [V]", {"Negative Stoichiometry": stoic_n})

# %%
# Define Model Equations

########################
# Interfacial reactions
########################

As_neg = 3*eps_s_n/R_n
As_pos = 3*eps_s_p/R_p
ex_current_n = k_n*pybamm.sqrt(c_e_n)*pybamm.sqrt(pybamm.surf(c_n))*pybamm.sqrt(C_max_neg-pybamm.surf(c_n))
eta_n = psi_s_n - psi_e_n - OCP_n 
j_flux_n = 2*ex_current_n*pybamm.sinh((F/(2*R_const*T_amb))*(eta_n))
eta_s = pybamm.PrimaryBroadcast(0, "Separator")
j_flux_s = pybamm.PrimaryBroadcast(0, "Separator")
ex_current_p = k_p*pybamm.sqrt(c_e_p)*pybamm.sqrt(pybamm.surf(c_p))*pybamm.sqrt(C_max_pos-pybamm.surf(c_p))
eta_p = psi_s_p - psi_e_p - OCP_p
j_flux_p = 2*ex_current_p*pybamm.sinh((F/(2*R_const*T_amb))*(eta_p))


eta = pybamm.concatenation(eta_n,eta_s,eta_p)
j_flux = pybamm.concatenation(j_flux_n,j_flux_s,j_flux_p)

j_ave_neg = I_app/(As_neg*F*L_n)
# %%
# Define SEI Flux

###########################
# Intercalation reactions
###########################


eta_SEI_n = psi_s_n - psi_e_n - OCP_SEI 

if SEI_con == 0:
    j_SEI_n = 0


if SEI_con == 1:
    j_SEI_n = -(I_app < SEI_tol)*(i_0_SEI/F) * pybamm.exp(  -(alpha_SEI*F/(R_const*T_amb))*eta_SEI_n)

if SEI_con == 2:
    
    L_SEI_var = (1/As_neg)*M_SEI*c_SEI_init*(1/rho_SEI)
    j_SEI_n = (I_app < SEI_tol)*(-c_sol_0/( (1/(k_SEI*pybamm.exp(-(alpha_SEI*F/(R_const*T_amb))*eta_SEI_n))) + (L_SEI_n/D_SEI) ) )
    


# %%
# Define Total Flux

j_tot_n = j_flux_n  + j_SEI_n
j_tot_s = j_flux_s 
j_tot_p = j_flux_p 
j_tot = pybamm.concatenation(j_tot_n,j_tot_s,j_tot_p)


deps_e_n_dt = As_neg*(M_SEI/rho_SEI)*j_SEI_n
P2D_Model.rhs[eps_e_n] = deps_e_n_dt
P2D_Model.initial_conditions[eps_e_n] = eps_e_n_init

# %%
# Define Electrolyte Potential Equation

################################
# Electrolyte Potential Equation
################################

kappa_eff_n = kappa_n* (eps_e_n**Brugg)
kappa_eff_s = kappa_s* (eps_e_s**Brugg)
kappa_eff_p = kappa_p* (eps_e_p**Brugg)
kappa_eff = pybamm.concatenation(kappa_eff_n,kappa_eff_s,kappa_eff_p)

kappa_d_n = (2*kappa_eff_n*R_const*T_amb*(1-t_li)/F)
kappa_d_s = (2*kappa_eff_s*R_const*T_amb*(1-t_li)/F)
kappa_d_p = (2*kappa_eff_p*R_const*T_amb*(1-t_li)/F)
kappa_d = pybamm.concatenation(kappa_d_n,kappa_d_s,kappa_d_p)

P2D_Model.algebraic[psi_e_n] = (L_n**2)*(pybamm.div(kappa_eff_n*pybamm.grad(psi_e_n)) - pybamm.div(kappa_d_n*(1/c_e_n)*pybamm.grad((c_e_n))) + As_neg*F*j_tot_n)
P2D_Model.initial_conditions[psi_e_n] = pybamm.Scalar(0)
P2D_Model.algebraic[psi_e_s] = (L_s**2)*(pybamm.div(kappa_eff_s*pybamm.grad(psi_e_s)) - pybamm.div(kappa_d_s*(1/c_e_s)*pybamm.grad((c_e_s))))

P2D_Model.initial_conditions[psi_e_s] = pybamm.Scalar(0)
P2D_Model.algebraic[psi_e_p] = (L_p**2)*(pybamm.div(kappa_eff_p*pybamm.grad(psi_e_p)) - pybamm.div(kappa_d_p*(1/c_e_p)*pybamm.grad((c_e_p))) + As_pos*F*j_tot_p)
P2D_Model.initial_conditions[psi_e_p] = pybamm.Scalar(0)

psi_e = pybamm.concatenation(psi_e_n,psi_e_s,psi_e_p)



P2D_Model.boundary_conditions[psi_e] = {"left": (pybamm.Scalar(0), "Dirichlet"),"right": (pybamm.Scalar(0), "Neumann")}

# %%
# Define Solid Potential Equation

################################
# Solid Potential Equation
################################

sigma_eff_n = sigma_n*(eps_s_n**Brugg)

P2D_Model.algebraic[psi_s_n] = (L_tot**2)*(pybamm.div(sigma_eff_n*pybamm.grad(psi_s_n)) -  F*As_neg*j_flux_n)
P2D_Model.boundary_conditions[psi_s_n] = {"left": (-I_app/sigma_eff_n, "Neumann"),"right": (pybamm.Scalar(0), "Neumann")}
P2D_Model.initial_conditions[psi_s_n] = OCP_n_init


psi_s_s = pybamm.PrimaryBroadcast(0, "Separator")

sigma_eff_p = sigma_p*(eps_s_p**Brugg)
P2D_Model.algebraic[psi_s_p] = (L_tot**2)*(pybamm.div(sigma_eff_p*pybamm.grad(psi_s_p)) -  F*As_pos*j_flux_p)
P2D_Model.boundary_conditions[psi_s_p] = {"left": (pybamm.Scalar(0), "Neumann"),"right": (-I_app/sigma_eff_p, "Neumann")}
P2D_Model.initial_conditions[psi_s_p] = OCP_p_init
psi_s = pybamm.concatenation(psi_s_n,psi_s_s,psi_s_p)

V_cell = pybamm.boundary_value(psi_s, "right") - pybamm.boundary_value(psi_s, "left")

psi_delta = psi_s - psi_e



# %%
# Calculate i_e


i_e_n = 0 + pybamm.IndefiniteIntegral(As_neg*F*j_tot_n,x_neg)
i_e_s = pybamm.PrimaryBroadcast(I_app,"Separator")
i_e_p = I_app + pybamm.IndefiniteIntegral(As_pos*F*j_tot_p,x_pos)

i_e = pybamm.concatenation(i_e_n,i_e_s,i_e_p)





# %%
# Define Electrolyte Concentration Equation

################################
# Electrolyte Concentration Equation
################################

N_e_n = -(eps_e_n**Brugg) * D_e * pybamm.grad(c_e_n) 
dcen_dt = (1 / eps_e_n) * (-pybamm.div(N_e_n) + (1 - t_li) * As_neg*j_tot_n)

N_e_s = -(eps_e_s**Brugg) * D_e * pybamm.grad(c_e_s)
dces_dt = (1 / eps_e_s) * (-pybamm.div(N_e_s)) 

N_e_p = -(eps_e_p**Brugg) * D_e * pybamm.grad(c_e_p)
dcep_dt = (1 / eps_e_p) * (-pybamm.div(N_e_p) + (1 - t_li) * As_pos*j_tot_p)

dce_dt = pybamm.concatenation(dcen_dt, dces_dt, dcep_dt)
P2D_Model.rhs[c_e] = dce_dt
P2D_Model.boundary_conditions[c_e] = {"left": (pybamm.Scalar(0), "Neumann"),"right": (pybamm.Scalar(0), "Neumann"),}
P2D_Model.initial_conditions[c_e] = c_e_init




# %%
# Define Solid Concentration Equation
######################
# Solid Phase Diffusion
######################

# Negative Electrode Diffusion Equation
N_n = -D_s_n*pybamm.grad(c_n)
P2D_Model.rhs[c_n] = -pybamm.div(N_n)
P2D_Model.initial_conditions[c_n] = SoC_init_n*C_max_neg
P2D_Model.boundary_conditions[c_n] = {"left": (pybamm.Scalar(0), "Neumann"), "right": (-j_flux_n/D_s_n, "Neumann")}

# Positive Electrode Diffusion Equation
N_p = -D_s_p*pybamm.grad(c_p)
P2D_Model.rhs[c_p] = -pybamm.div(N_p)
P2D_Model.initial_conditions[c_p] = SoC_init_p*C_max_pos
P2D_Model.boundary_conditions[c_p] = {"left": (pybamm.Scalar(0), "Neumann"), "right": (-j_flux_p/D_s_p, "Neumann")}

# %%
# Define SEI Growth Model
P2D_Model.rhs[c_SEI_n] = -As_neg*j_SEI_n
P2D_Model.initial_conditions[c_SEI_n] = c_SEI_init

#L_SEI_n = (1/As_neg)*M_SEI*c_SEI_n*(1/rho_SEI)

d_L_SEI_dt = -M_SEI*j_SEI_n*(1/rho_SEI)

R_SEI_n = L_SEI_n/kappa_SEI 

d_R_SEI_n_dt =  -M_SEI*j_SEI_n*(1/rho_SEI)*(1/kappa_SEI )

eta_n += - F*R_SEI_n*j_ave_neg
eta_SEI_n += - F*R_SEI_n*j_ave_neg 

P2D_Model.rhs[L_SEI_n] = -M_SEI*j_SEI_n*(1/rho_SEI)
P2D_Model.initial_conditions[L_SEI_n] = (1/As_neg)*M_SEI*c_SEI_init*(1/rho_SEI)


if SEI_con == 2:
    #L_SEI_var += L_SEI_n
    L_SEI_var += d_L_SEI_dt
    #R_SEI_n_var += d_R_SEI_n_dt
    #L_SEI_var = 0
    #L_SEI_var = L_SEI_n

# %%
# Define Model Variables

c_e_0 = pybamm.boundary_value(c_e_n,"left")
c_e_n_s_1 = pybamm.boundary_value(c_e_n,"right")
c_e_n_s_2 = pybamm.boundary_value(c_e_s,"left")
c_e_s_p_1 = pybamm.boundary_value(c_e_s,"right")
c_e_s_p_2 = pybamm.boundary_value(c_e_p,"left")
c_e_p_end = pybamm.boundary_value(c_e_p,"right")
c_s_n_surf = pybamm.surf(c_n)
c_s_p_surf = pybamm.surf(c_p)
c_s_p_surf = pybamm.surf(c_p)



P2D_Model.variables = {"Current function [A]":I_app,"Negative Particle Concentration [mols/m^3]":c_n,"Positive Particle Concentration [mols/m^3]":c_p,"Negative Particle Length":r_neg,"Positive Particle Length":r_pos ,"Electrolyte Concentration [mols/m^3]":c_e,"Differential Voltage [V]": psi_delta,"Cell Voltage [V]": V_cell,"Intercalation Flux [mols/m^3]": j_flux,"Total Flux [mols/m^3]": j_tot,"Electrolyte Potential [V]":psi_e,"Length of Negative Electrode [m]": x_neg,"Length of Separator [m]": x_sep,"Length of Positive Electrode [m]": x_pos,
                       "Electorlyte Concentration x = 0":c_e_0,"Electorlyte Concentration x = L_n-":c_e_n_s_1,"Electorlyte Concentration x = L_n+":c_e_n_s_2,"Electorlyte Concentration x = L_n+L_s-":c_e_s_p_1,"Electorlyte Concentration x = L_n+L_s+":c_e_s_p_2,"Electorlyte Concentration x = L_tot":c_e_p_end,
                       "Solid Surface Concentration Negative":c_s_n_surf,"Solid Surface Concentration Positive":c_s_p_surf,"Negative Intercalation Flux [mols/m^2s]":j_SEI_n,"Negative SEI Concentration [mols/m^3]":c_SEI_n,"Negative SEI Thickness [m]":L_SEI_n,"Negative SEI Resistance [Ohm m^2]":R_SEI_n
                       ,"Negative Electrode Porosity Variable [-]":eps_e_n}

Summary_Variables_List = ["Current function [A]", "Negative Particle Concentration [mols/m^3]","Positive Particle Concentration [mols/m^3]","Negative Particle Length","Positive Particle Length","Electrolyte Concentration [mols/m^3]","Differential Voltage [V]","Cell Voltage [V]","Intercalation Flux [mols/m^3]","Total Flux [mols/m^3]","Electrolyte Potential [V]"]

P2D_Model.summary_variables = {"Current function [A]":I_app,"Negative Particle Concentration [mols/m^3]":c_n,"Positive Particle Concentration [mols/m^3]":c_p,"Negative Particle Length":r_neg,"Positive Particle Length":r_pos ,"Electrolyte Concentration [mols/m^3]":c_e,"Differential Voltage [V]": psi_delta,"Cell Voltage [V]": V_cell,"Intercalation Flux [mols/m^3]": j_flux,"Total Flux [mols/m^3]": j_tot,"Electrolyte Potential [V]":psi_e}
# %%

# Import Battery Parameters
param =1
#n_parameters = 42
n_parameters = 50

import Battery_Parameters_MDN_2 as param_file

Model_Parameters = param_file.Battery_Parameters_P2D(param,n_parameters)
F_const = 96487  # Faraday Constant [C/mol]

#Set Parameter Values
Q_pos_initial =  (1/3600)*F_const*Model_Parameters[16]* Model_Parameters[19]*Model_Parameters[17]*Model_Parameters[20]*np.abs(1.0 - Model_Parameters[17])
Q_neg_initial = (1/3600)*F_const*Model_Parameters[7]*Model_Parameters[10]*Model_Parameters[8]*Model_Parameters[11]*np.abs(Model_Parameters[8] - 0)
max_cap = np.max([Q_pos_initial,Q_neg_initial])

param_values = pybamm.ParameterValues({
    "Bruggeman Coefficient [-]": Model_Parameters[1],
    "Length of Separator [m]": Model_Parameters[2],
    "Initial Electrolyte Concentration [mols/m^3]": Model_Parameters[3],
    "Elecrolyte Ionic Diffusion Coefficient [m^2/s]": Model_Parameters[4],
    "Transport Coefficient [-]": Model_Parameters[5],
    "Separator Porosity": Model_Parameters[6],
    "Length of Negative Electrode [m]": Model_Parameters[7],
    "Negative Particle Initial State of Charge [-]": Model_Parameters[8],
    "Negative Electrode Porosity": Model_Parameters[9],
    "Negative Particle Solid Volume Fraction [-]": Model_Parameters[10],
    "Negative Particle Maximum Concentration [mols/m^3]": Model_Parameters[11],
    "Negative Particle Radius [m]": Model_Parameters[12],
    "Negative Solid Diffusion Coefficient [m^2/s]": Model_Parameters[13],
    "Negative Electrode Electronic Conductivity [S/m]": Model_Parameters[14],
    "Negative Rate Constant [mols/m^2s]": Model_Parameters[15],
    "Length of Positive Electrode [m]": Model_Parameters[16],
    "Negative Double Layer Capacitance [F/m^2]": Model_Parameters[26],
    "Positive Particle Initial State of Charge [-]": Model_Parameters[17],
    "Positive Electrode Porosity": Model_Parameters[18],
    "Positive Particle Solid Volume Fraction [-]": Model_Parameters[19],
    "Positive Particle Maximum Concentration [mols/m^3]": Model_Parameters[20],
    "Positive Particle Radius [m]": Model_Parameters[21],
    "Positive Solid Diffusion Coefficient [m^2/s]": Model_Parameters[22],
    "Positive Electrode Electronic Conductivity [S/m]": Model_Parameters[23],
    "Positive Rate Constant [mols/m^2s]": Model_Parameters[24],
    "Positive Double Layer Capacitance [F/m^2]": Model_Parameters[27],
    "Total Cell Length [m]": Model_Parameters[7] + Model_Parameters[2] + Model_Parameters[16],
    "Negative Electrolyte Conductivity [S/m]": kappa_function,
    "Positive Electrolyte Conductivity [S/m]": kappa_function,
    "Separator Electrolyte Conductivity [S/m]": kappa_function,
    "Positive Electrode OCP [V]": ocp_pos,
    "Negative Electrode OCP [V]": ocp_neg,
    "Initial Negative Electrode OCP [V]": ocp_neg,
    "Initial Positive Electrode OCP [V]": ocp_pos,
    "Current function [A]": Applied_Current_Function,
    "Faraday Constant [cols/mol]": 96487,
    "Low Voltage Cutt-off [V]": 2.8,
    "Ambient Temperature [K]": 298.15,
    "Ideal Gas Constant [J/(K*mol)]": 8.314,
    "SEI Conductivity [S/m]": Model_Parameters[32],
    "SEI Molar Mass [kg/mol]": Model_Parameters[33],
    "SEI Layer Density [kg/m^3]": Model_Parameters[34],
    "SEI Exchange Current [A/m^2]": Model_Parameters[35],
    "SEI OCP [V]": Model_Parameters[36],
    "Charge transfer coefficient SEI": Model_Parameters[37],
    "Initial Negative SEI Concentration [mols/m^3]": Model_Parameters[31],
    "SEI Reaction Rate constant [m/s]": Model_Parameters[39],
    "Concentration of Solvent in bulk electrolyte [mol/m^3]": Model_Parameters[40],
    "EC diffusion coefficient [m^2/s]": Model_Parameters[41]

})

# Define Geometry
geometry = {"Negative Electrode": {x_neg: {"min": pybamm.Scalar(0), "max": L_n}},"Negative Particle": {r_neg: {"min": pybamm.Scalar(0), "max": R_n}},"Separator": {x_sep: {"min": L_n, "max": L_n+L_s}},"Positive Electrode": {x_pos: {"min": L_n+L_s, "max": L_n+L_s+L_p}},"Positive Particle": {r_pos: {"min": pybamm.Scalar(0), "max": R_p}}}

# Set Termination Events
# Define Model Events
P2D_Model.events = [
    
    pybamm.Event("Low Cut-off Voltage Condition", V_cell - V_cell_low_cut_off),

]

# Process Model
param_values.process_model(P2D_Model)
param_values.process_geometry(geometry)
# %%
# Disctretise Model
submesh_types = {
    "Negative Electrode": pybamm.Uniform1DSubMesh,
    "Negative Particle": pybamm.Uniform1DSubMesh,
    "Separator": pybamm.Uniform1DSubMesh,
    "Positive Electrode": pybamm.Uniform1DSubMesh,
    "Positive Particle": pybamm.Uniform1DSubMesh,
}
var_pts = {x_neg: 20,r_neg:10,x_sep:20, x_pos: 20, r_pos: 10}

mesh = pybamm.Mesh(geometry, submesh_types, var_pts)

spatial_methods = {
    "Negative Electrode": pybamm.FiniteVolume(),
    "Negative Particle": pybamm.FiniteVolume(),
    "Separator": pybamm.FiniteVolume(),
    "Positive Electrode": pybamm.FiniteVolume(),
    "Positive Particle": pybamm.FiniteVolume(),
}
disc = pybamm.Discretisation(mesh, spatial_methods)
disc.process_model(P2D_Model,inplace=False)



# %%
# Solve Model
# Solve the model

drive_cycle_current= np.column_stack([time_drive_cycle.T, current_drive_cycle.T])
Experiment = pybamm.Experiment([pybamm.step.current(drive_cycle_current)])



safe_solver = pybamm.CasadiSolver(atol=1e-6, rtol=1e-6, mode="safe",return_solution_if_failed_early=True,dt_max=None)



sim = pybamm.Simulation(P2D_Model, experiment=Experiment,solver = safe_solver,parameter_values=param_values,submesh_types={
    "Negative Electrode": pybamm.Uniform1DSubMesh,
    "Negative Particle": pybamm.Uniform1DSubMesh,
    "Separator": pybamm.Uniform1DSubMesh,
    "Positive Electrode": pybamm.Uniform1DSubMesh,
    "Positive Particle": pybamm.Uniform1DSubMesh,
},geometry=
{"Negative Electrode": {x_neg: {"min": pybamm.Scalar(0), "max": L_n}},
 "Negative Particle": {r_neg: {"min": pybamm.Scalar(0), "max": R_n}},
 "Separator": {x_sep: {"min": L_n, "max": L_n+L_s}},
 "Positive Electrode": {x_pos: {"min": L_n+L_s, "max": L_n+L_s+L_p}},
 "Positive Particle": {r_pos: {"min": pybamm.Scalar(0), "max": R_p}}
 },var_pts=
 {x_neg: 18,
  r_neg:5,
  x_sep:4,
    x_pos: 18,
      r_pos: 5
      },
        spatial_methods = {
    "Negative Electrode": pybamm.FiniteVolume(),
    "Negative Particle": pybamm.FiniteVolume(),
    "Separator": pybamm.FiniteVolume(),
    "Positive Electrode": pybamm.FiniteVolume(),
    "Positive Particle": pybamm.FiniteVolume(),
})

sim.parameter_values.update({"Current function [A]":current_drive_cycle})

solution_1 = sim.solve(calc_esoh=False,showprogress=True)


print(solution_1.termination)

# %%
# Post Processing

c_n = solution_1["Negative Particle Concentration [mols/m^3]"]
r_neg = solution_1["Negative Particle Length"]
r_pos = solution_1["Positive Particle Length"]
x_neg = solution_1["Length of Negative Electrode [m]"]
x_sep = solution_1["Length of Separator [m]"]
x_pos = solution_1["Length of Positive Electrode [m]"]
c_p = solution_1["Positive Particle Concentration [mols/m^3]"]
c_e = solution_1["Electrolyte Concentration [mols/m^3]"]
psi_e = solution_1["Electrolyte Potential [V]"]
psi_delta = solution_1["Differential Voltage [V]"]
j_flux = solution_1["Intercalation Flux [mols/m^3]"]
j_tot = solution_1["Total Flux [mols/m^3]"]
V_cell = solution_1["Cell Voltage [V]"]
I_app = solution_1["Current function [A]"]
c_e_0 = solution_1["Electorlyte Concentration x = 0"]
c_e_n_s_1 = solution_1["Electorlyte Concentration x = L_n-"]
c_e_n_s_2 = solution_1["Electorlyte Concentration x = L_n+"]
c_e_s_p_1 = solution_1["Electorlyte Concentration x = L_n+L_s-"]
c_e_s_p_2 = solution_1["Electorlyte Concentration x = L_n+L_s+"]
c_e_p_end = solution_1["Electorlyte Concentration x = L_tot"]
c_s_n_surf = solution_1["Solid Surface Concentration Negative"]
c_s_p_surf = solution_1["Solid Surface Concentration Positive"]
j_SEI_n = solution_1["Negative Intercalation Flux [mols/m^2s]"]
c_SEI_n = solution_1["Negative SEI Concentration [mols/m^3]"]
L_SEI_n = solution_1["Negative SEI Thickness [m]"]
R_SEI_n = solution_1["Negative SEI Resistance [Ohm m^2]"]
eps_e_n = solution_1["Negative Electrode Porosity Variable [-]"]
#"Negative Electrode Porosity Variable [-]":eps_e_n


# %%
# Get Data into Numpy Arrays
Negative_conc = c_n.data
Positive_conc = c_p.data
r_n_vector = r_neg.data
r_p_vector = r_pos.data
x_neg_vector = x_neg.data
x_sep_vector = x_sep.data
x_pos_vector = x_pos.data
c_e_conc_model_1 = c_e.data
c_e_0_vector = c_e_0.data
c_e_n_s_1_vector = c_e_n_s_1.data
c_e_n_s_2_vector = c_e_n_s_2.data
c_e_s_p_1_vector = c_e_s_p_1.data
c_e_s_p_2_vector = c_e_s_p_2.data
c_e_p_end_vector = c_e_p_end.data
c_s_n_surf_vector_model_1 = c_s_n_surf.data
c_s_p_surf_vector_model_1 = c_s_p_surf.data
c_s_n_ave_vector_model_1 = np.mean(Negative_conc,axis=0)
c_s_p_ave_vector_model_1 = np.mean(Positive_conc,axis=0)
j_flux_data_model_1 = j_flux.data
j_tot_data_model_1 = j_tot.data
j_dl_data_model_1 = j_tot_data_model_1 - j_flux_data_model_1
psi_delta_data = psi_delta.data
psi_e_data = psi_e.data
j_SEI_n_data_model_1 = j_SEI_n.data
c_SEI_n_data_model_1 = c_SEI_n.data
L_SEI_n_data_model_1 = L_SEI_n.data
R_SEI_n_data_model_1 = R_SEI_n.data
eps_e_n_model_1 = eps_e_n.data
#V_cell_vector = V_cell.data
time_vector_model_1 = solution_1.t
V_cell_vector_model_1 = V_cell.data
I_app_vector_model_1 = I_app.data
x_tot_vector = np.concatenate([x_neg_vector[:,0],x_sep_vector[:,0],x_pos_vector[:,0]])

x_elements = c_e_conc_model_1[:,0].size

x_elements_vector = np.linspace(1,x_elements,x_elements)
x_tot_vector = np.concatenate([x_neg_vector[:,0],x_sep_vector[:,0],x_pos_vector[:,0]])

# %%

plt.figure(1)
plt.grid
plt.plot(time_vector_model_1,V_cell_vector_model_1,'b--',label = 'PyBaMM Custom Model')
plt.plot(time_drive_cycle.flatten(),voltage_drive_cycle.flatten(),'k',label = 'MATLAB Custom Model')
plt.plot(time_vector_pybamm,V_cell_vector_pybamm,'r-.',label = 'PyBaMM Built-in Model')
plt.xlabel("Time [s]")
plt.ylabel("Cell Voltage [V]")
plt.legend()
plt.show()



plt.figure(2)
plt.grid
plt.plot(time_vector_model_1,I_app_vector_model_1,'b--',label = 'PyBaMM Custom Model')
plt.plot(time_drive_cycle.flatten(),current_drive_cycle.flatten(),'k',label = 'MATLAB Custom Model')
plt.plot(time_vector_pybamm,I_app_pybamm_vector,'r-.',label = 'PyBaMM Built-in Model')
plt.xlabel("Time [s]")
plt.ylabel("Applied Current [A/m^2]")
plt.legend()
plt.show()




# %%
# Export Data For Analysis
output_cell_voltage = {
    "time (s)": time_vector_model_1,
  "voltage (v)": V_cell_vector_model_1
}
output_cell_current = {
    "time (s)": time_vector_model_1,
  "current (A/m^2)": I_app_vector_model_1
}
output_c_e = {
    "time (s)": time_vector_model_1,
  "c_e_0 (mols/m^3)": c_e_0_vector,
  "c_e_n_s (mols/m^3)":c_e_n_s_1_vector,
  "c_e_s_p (mols/m^3)":c_e_s_p_1_vector,
  "c_e_p_end (mols/m^3)":c_e_p_end_vector
}
output_c_se = {"time (s)": time_vector_model_1,
  "c_se_0 (mols/m^3)": c_s_n_surf_vector_model_1[0,:],
  "c_se_n_s (mols/m^3)":c_s_n_surf_vector_model_1[-1,:],
  "c_se_s_p (mols/m^3)":c_s_p_surf_vector_model_1[0,:],
  "c_se_p_end (mols/m^3)":c_s_p_surf_vector_model_1[-1,:]}

output_c_s_ave = {"time (s)": time_vector_model_1,
  "c_s_ave_0 (mols/m^3)": c_s_n_ave_vector_model_1[0,:],
  "c_s_ave_n_s (mols/m^3)":c_s_n_ave_vector_model_1[-1,:],
  "c_s_ave_s_p (mols/m^3)":c_s_p_ave_vector_model_1[0,:],
  "c_s_ave_end (mols/m^3)":c_s_p_ave_vector_model_1[-1,:]}





# %%
# Export Data

#load output data into a DataFrame object:
df = pd.DataFrame(output_cell_voltage)

# writing output data frame to a CSV file
df.to_csv('Pybamm_custom_CCCV_voltage_ageing_LCO.csv')

#load output data into a DataFrame object:
df = pd.DataFrame(output_cell_current)

# writing output data frame to a CSV file
df.to_csv('Pybamm_custom_CCCV_current_ageing_LCO.csv')

#load output data into a DataFrame object:
df = pd.DataFrame(output_c_e)

# writing output data frame to a CSV file
df.to_csv('Pybamm_custom_CCCV_c_e_ageing_LCO.csv')

#load output data into a DataFrame object:
df = pd.DataFrame(output_c_se)

# writing output data frame to a CSV file
df.to_csv('Pybamm_custom_CCCV_c_se_ageing_LCO.csv')

#load output data into a DataFrame object:
df = pd.DataFrame(output_c_s_ave)


# writing output data frame to a CSV file
df.to_csv('Pybamm_custom_CCCV_c_s_ave_ageing_LCO.csv')



# %%
