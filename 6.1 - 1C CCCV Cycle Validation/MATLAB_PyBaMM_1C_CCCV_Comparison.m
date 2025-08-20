%% Script to conduct one CCCV cycle
clear;
close all;


%% Model Selection

param  =1; %% param = 1 for LCO, param = 2 for NMC

%% P2D Modelling Parameters
Comparison_points = 1000;
time_step = 1;
order = 5;
N_tot = 40;
N_shell = 5;
V_cut_off_low = 3.35;
V_cut_off_high_interim = 4.1;
V_cut_off_high_80 = 3.9;
V_cut_off_high_60 = 3.84;
V_cut_off_high_40 = 3.81;
V_cut_off_high_end = V_cut_off_high_interim;
V_cut_off_high = V_cut_off_high_interim;
V_cut_off_test_vectors = [V_cut_off_high_60  V_cut_off_high_interim];
%I_app_cut_off = -0.1; %A/m^2
I_app_cut_off_percent = 0.05; % Percentage of the 1C current
Cut_off_capacity = 0.1;
tol_CV = 5e-4;
N_cycles = 1;
Test_cycles = [3 6 9 12 15 18];
n_parameters = 46;
C_rate_non_test = 1;
C_rate_test = 1;


%% Choose Degredation Models
SEI = 0;
li_pl = 0;
coupling = 2;
Ageing_models = [SEI li_pl coupling];



%% Simulation Setting
%sim_type = 2; %% sim_type = 1 for discharge, sim_type =2 for charge, sim_type =3 discharge/charge
C_rate = 1;
t_rest_discharge = 5*60; %%% Rest time after discharge
t_rest_charge = 5*60; %%% Rest time after CV charge

%% Import Parameters
Battery_Parameters_P2D_MATLAB_2;
Model_Parameters = zeros(1,n_parameters);

%%% Seperator Parameters
Model_Parameters(1) = A; % [m] %Cell Area
Model_Parameters(2) = brugg;
Model_Parameters(3) = L_sep ; % Length of separator [m] 
Model_Parameters(4)= C_e_initial; %Initial Electrolyte Lithium Concentration [mol/[m^3]] 
Model_Parameters(5) = D_elec; % Diffusion Coefficient for Lithium in electrolyte [[m^2]/s] %Defined in model params
Model_Parameters(6) = t_li; % Lithium Ion Transference Number [unitless]
Model_Parameters(7) = eps_sep_elec; %Electrolyte Volume Fraction
Model_Parameters(29) =dlnfdlnce_sep;

%%%  Anode Parameters
Model_Parameters(8) = L_neg;
Model_Parameters(9) = soc_initial_neg;
Model_Parameters(10) = eps_neg_elec; % Electrolyte Volume Fraction (porosity)
Model_Parameters(11)= eps_neg; % Solid Volume Fraction
Model_Parameters(12) = C_max_neg ; % Maximum Solid Lithium Concentration [mol/[m^3]]
Model_Parameters(13)= Rs_neg; % Electrode Particle Size [m] %changed
Model_Parameters(14)= D_neg_s; % Solid Diffusion Coefficient [m^2/s]
Model_Parameters(15)= sigma_neg ; % Electrical Conductivity [S/m]
Model_Parameters(16) = K_0_neg; % Negative Rate Constant %changed
Model_Parameters(27) = Cdl_neg;


%%% Cathode Parameters
Model_Parameters(17)= L_pos;
Model_Parameters(18)= soc_initial_pos;
Model_Parameters(19) = eps_pos_elec; %% Electrolyte Volume Fraction 
Model_Parameters(20) = eps_pos; % Solid Volume Fraction
Model_Parameters(21) =C_max_pos ; % Maximum Solid Lithium Concentration [mol/[m^3]]
Model_Parameters(22) = Rs_pos; % Electrode Particle Size [m] %changed
Model_Parameters(23) = D_pos_s; % Solid Diffusion Coefficient [m^2/s]
Model_Parameters(24) = sigma_pos ; % Electrical Conductivity [S/m] From COMSOL
Model_Parameters(25) = K_0_pos; % Negative Rate Constant %changed
Model_Parameters(26) =kappa_const;
Model_Parameters(28) = Cdl_pos;

%%% SEI Parameters
Model_Parameters(30) = L_SEI; 
Model_Parameters(31) = kappa_SEI;
Model_Parameters(32) = M_SEI;
Model_Parameters(33) = rho_SEI;
Model_Parameters(34) = i_0_SEI;
Model_Parameters(35) = OCP_SEI;
Model_Parameters(36) = alpha_SEI;
Model_Parameters(43) = c_sol_0;
Model_Parameters(44) = k_SEI;
Model_Parameters(45) = D_SEI;

%%% Plating Parameters
Model_Parameters(37) = k_li;
Model_Parameters(38) = alpha_pl_neg;
Model_Parameters(39) = alpha_pl_pos;
Model_Parameters(40) = kappa_li;
Model_Parameters(41) = M_Li;
Model_Parameters(42) = rho_Li;
Model_Parameters(46) = OCP_li;

%% EIS Stored Data
N_tot_store = 29;
soc_points = length(V_cut_off_test_vectors)*(length(Test_cycles));
C_e_EIS = zeros(soc_points,N_tot_store);
C_s_full_EIS = zeros(soc_points,N_shell,N_tot_store);
C_se_EIS = zeros(soc_points,N_tot_store);
c_SEI_EIS = zeros(soc_points,N_tot_store);
c_li_EIS = zeros(soc_points,N_tot_store);
psi_s_EIS = zeros(soc_points,N_tot_store);
psi_e_EIS = zeros(soc_points,N_tot_store);
j_tot_EIS = zeros(soc_points,N_tot_store);

%% Capacity Data
q_lost = zeros(1,N_cycles);


%% Simulations
Test_cycle_count = 0;

%% Import Drive Cycle

if param == 1
    Input_Data = xlsread("Pybamm_custom_CCCV_current_ageing_LCO.csv");
end

if param == 2
    Input_Data = xlsread("Pybamm_custom_CCCV_current_ageing_NMC.csv");
end

time_vector_input = Input_Data(:,2);
Current_Vector = Input_Data(:,3);


%% Forward Copute simulation
sim_type = 0;
C_e_init = zeros(1,N_tot);
C_s_full_init = zeros(1,N_shell,N_tot);
C_se_init = zeros(1,N_tot);
I_app_init = Current_Vector';
c_SEI_init = zeros(1,N_tot);
c_li_init = zeros(1,N_tot);
psi_s_init = zeros(1,N_tot);
psi_e_init = zeros(1,N_tot);
j_tot_init = zeros(1,N_tot);
eps_sol_init = zeros(1,N_tot);


[V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,c_SEI,c_li,I_app,psi_s,psi_e,j_tot,eps_sol,L_SEI,R_SEI,j_SEI,eps_sol_store] = P2D_function_general(Model_Parameters,N_tot,N_shell,order,param,time_step,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);

%%% Update local variables
   
%%% Update Global Data
time_vector_total = time_vector;
I_app_total =I_app;
V_cell_total = V_cell;
C_e_total = C_e;
C_se_total = C_se;
C_s_ave_total = squeeze(mean(C_s_full,2));
c_SEI_tot = c_SEI;
L_SEI_tot = L_SEI;
R_SEI_tot = R_SEI;
j_SEI_tot = j_SEI;
eps_sol_store_tot = eps_sol_store;
j_tot_total = j_tot;


%% Plot Results

close all;

figure(1);
plot(time_vector_total,I_app_total,'r');
grid on;
xlabel('Time (s)');
ylabel('Cell Current (A/m^2)');

figure(2)
plot(time_vector_total,V_cell_total,'r');
grid on;
xlabel('Time (s)');
ylabel('Cell Voltage (V)');



%%% Find Points of boundaries
idx_disc = find(C_se_total(1,:)==0);
N_neg = idx_disc(1)-1;
N_sep = length(idx_disc);


figure(3);
hold on;
plot(time_vector_total,C_e_total(:,1)','r');
C_e_n_s_matlab = mean([C_e_total(:,N_neg)'; C_e_total(:,N_neg+1)'],1);
plot(time_vector_total,C_e_n_s_matlab,'b');
C_e_s_p_matlab = mean([C_e_total(:,N_neg+N_sep)'; C_e_total(:,N_neg+N_sep+1)'],1);
plot(time_vector_total,C_e_s_p_matlab,'g');
plot(time_vector_total,C_e_total(:,end)','k');
grid on;
xlabel('Time (s)');
ylabel('Electrolyte Concentration (mols/m^3)');

figure(4);
hold on;
plot(time_vector_total,C_se_total(:,1)','r');
plot(time_vector_total,C_se_total(:,N_neg)','b');
plot(time_vector_total,C_se_total(:,N_neg+N_sep+1)','g');
plot(time_vector_total,C_se_total(:,end)','k');
grid on;
xlabel('Time (s)');
ylabel('Electrode Surface Concentration (mols/m^3)');

figure(5);
hold on;
plot(time_vector_total,C_s_ave_total(:,1)','r');
plot(time_vector_total,C_s_ave_total(:,N_neg)','b');
plot(time_vector_total,C_s_ave_total(:,N_neg+N_sep+1)','g');
plot(time_vector_total,C_s_ave_total(:,end)','k');
grid on;
xlabel('Time (s)');
ylabel('Electrode Average Concentration (mols/m^3)');


%% Import Drive Cycle Data


if param == 1
    Voltage_data_pybamm = xlsread("Pybamm_custom_CCCV_voltage_ageing_LCO.csv");
    Current_data_pybamm = xlsread("Pybamm_custom_CCCV_current_ageing_LCO.csv");
    C_e_data_pybamm = xlsread("Pybamm_custom_CCCV_c_e_ageing_LCO.csv");
    C_se_data_pybamm = xlsread("Pybamm_custom_CCCV_c_se_ageing_LCO.csv");
    C_s_ave_data_pybamm = xlsread("Pybamm_custom_CCCV_c_s_ave_ageing_LCO.csv");

end

if param == 2
    Voltage_data_pybamm = xlsread("Pybamm_custom_CCCV_voltage_ageing_NMC.csv");
    Current_data_pybamm = xlsread("Pybamm_custom_CCCV_current_ageing_NMC.csv");
    C_e_data_pybamm = xlsread("Pybamm_custom_CCCV_c_e_ageing_NMC.csv");
    C_se_data_pybamm = xlsread("Pybamm_custom_CCCV_c_se_ageing_NMC.csv");
    C_s_ave_data_pybamm = xlsread("Pybamm_custom_CCCV_c_s_ave_ageing_NMC.csv");

end

%}
%% Plot Pybamm data
figure(1);
hold on;
plot(Current_data_pybamm(:,2)',Current_data_pybamm(:,3)','b--');
grid on;
xlabel('Time (s)');
ylabel('Cell Current (A/m^2)');
legend('MATLAB Solver','PyBaMM Solver',Location='northeast');
fontsize(figure(1),'increase');


figure(2)
hold on;
plot(Voltage_data_pybamm(:,2)',Voltage_data_pybamm(:,3)','b--');
grid on;
xlabel('Time (s)');
ylabel('Cell Voltage (V)');
legend('MATLAB Solver','PyBaMM Solver',Location='northwest');
fontsize(figure(2),'increase');


figure(3);
hold on;
plot(C_e_data_pybamm(:,2)',C_e_data_pybamm(:,3)','r--');
plot(C_e_data_pybamm(:,2)',C_e_data_pybamm(:,4)','b--');
plot(C_e_data_pybamm(:,2)',C_e_data_pybamm(:,5)','g--');
plot(C_e_data_pybamm(:,2)',C_e_data_pybamm(:,6)','k--');
grid on;
xlabel('Time (s)');
ylabel('Electrolyte Concentration (mols/m^3)');
fontsize(figure(3),'increase');


figure(4);
hold on;
plot(C_se_data_pybamm(:,2)',C_se_data_pybamm(:,3)','r--');
plot(C_se_data_pybamm(:,2)',C_se_data_pybamm(:,4)','b--');
plot(C_se_data_pybamm(:,2)',C_se_data_pybamm(:,5)','g--');
plot(C_se_data_pybamm(:,2)',C_se_data_pybamm(:,6)','k--');
grid on;
xlabel('Time (s)');
ylabel('Surface Lithium Concentration (mols/m^3)');
fontsize(figure(4),'increase');


figure(5);
hold on;
plot(C_s_ave_data_pybamm(:,2)',C_s_ave_data_pybamm(:,3)','r--');
plot(C_s_ave_data_pybamm(:,2)',C_s_ave_data_pybamm(:,4)','b--');
plot(C_s_ave_data_pybamm(:,2)',C_s_ave_data_pybamm(:,5)','g--');
plot(C_s_ave_data_pybamm(:,2)',C_s_ave_data_pybamm(:,6)','k--');
grid on;
xlabel('Time (s)');
ylabel('Average Lithium Concentration (mols/m^3)');
fontsize(figure(5),'increase');





%% Calculate Errors
RMSE_Current_percentage = mean(abs(I_app_total(1:end-1)-Current_data_pybamm(:,3)'));

RMSE_Voltage_percentage = mean(abs(V_cell_total(1:end-1)-Voltage_data_pybamm(:,3)')./Voltage_data_pybamm(:,3)');
RMSE_Voltage_abs = mean(abs(V_cell_total(1:end-1)-Voltage_data_pybamm(:,3)'));

fprintf('The Precentage Voltage RMSE = %.5f\n', RMSE_Voltage_percentage);

RMSE_C_e_cc_neg = mean(abs(C_e_total(1:end-1,1)'-C_e_data_pybamm(:,3)')./C_e_data_pybamm(:,3)');
RMSE_C_e_neg_sep = mean(abs(C_e_n_s_matlab(1:end-1)-C_e_data_pybamm(:,4)')./C_e_data_pybamm(:,4)');
RMSE_C_e_sep_pos = mean(abs(C_e_s_p_matlab(1:end-1)-C_e_data_pybamm(:,5)')./C_e_data_pybamm(:,5)');
RMSE_C_e_pos_cc = mean(abs(C_e_total(1:end-1,end)'-C_e_data_pybamm(:,6)')./C_e_data_pybamm(:,6)');

RMSE_C_e = mean([RMSE_C_e_cc_neg RMSE_C_e_neg_sep RMSE_C_e_sep_pos RMSE_C_e_pos_cc]);

fprintf('The Precentage C_e RMSE = %.5f\n', RMSE_C_e);

RMSE_C_se_cc_neg = mean(abs(C_se_total(1:end-1,1)'-C_se_data_pybamm(:,3)')./C_se_data_pybamm(:,3)');
RMSE_C_se_neg_sep = mean(abs(C_se_total(1:end-1,N_neg)'-C_se_data_pybamm(:,4)')./C_se_data_pybamm(:,4)');
RMSE_C_se_sep_pos = mean(abs(C_se_total(1:end-1,N_neg+N_sep+1)'-C_se_data_pybamm(:,5)')./C_se_data_pybamm(:,5)');
RMSE_C_se_pos_cc = mean(abs(C_se_total(1:end-1,end)'-C_se_data_pybamm(:,6)')./C_se_data_pybamm(:,6)');

RMSE_C_se = mean([RMSE_C_se_cc_neg RMSE_C_se_neg_sep RMSE_C_se_sep_pos RMSE_C_se_pos_cc]);

fprintf('The Precentage C_se RMSE = %.5f\n', RMSE_C_se);

RMSE_C_s_ave_cc_neg = mean(abs(C_s_ave_total(1:end-1,1)'-C_s_ave_data_pybamm(:,3)')./C_s_ave_data_pybamm(:,3)');
RMSE_C_s_ave_neg_sep = mean(abs(C_s_ave_total(1:end-1,N_neg)'-C_s_ave_data_pybamm(:,4)')./C_s_ave_data_pybamm(:,4)');
RMSE_C_s_ave_sep_pos = mean(abs(C_s_ave_total(1:end-1,N_neg+N_sep+1)'-C_s_ave_data_pybamm(:,5)')./C_s_ave_data_pybamm(:,5)');
RMSE_C_s_ave_pos_cc = mean(abs(C_s_ave_total(1:end-1,end)'-C_s_ave_data_pybamm(:,6)')./C_s_ave_data_pybamm(:,6)');

RMSE_C_s_ave = mean([RMSE_C_s_ave_cc_neg RMSE_C_s_ave_neg_sep RMSE_C_s_ave_sep_pos RMSE_C_s_ave_pos_cc]);

fprintf('The Precentage C_se RMSE = %.5f\n', RMSE_C_se);



%}
%% Export Figures
dpi = 300;
if param == 1
    figure(1);
    print('1C CCCV Current LCO', '-dpng', ['-r', num2str(dpi)]);
    figure(2);
    print('1C CCCV Voltage LCO', '-dpng', ['-r', num2str(dpi)]);
    figure(3);
    print('C_e 1C CCCV LCO', '-dpng', ['-r', num2str(dpi)]);
    figure(4);
    print('C_se 1C CCCV LCO', '-dpng', ['-r', num2str(dpi)]);
    figure(5);
    print('C_ave 1C CCCV LCO', '-dpng', ['-r', num2str(dpi)]);
end 


if param == 2
    figure(1);
    print('1C CCCV Current NMC', '-dpng', ['-r', num2str(dpi)]);
    figure(2);
    print('1C CCCV Voltage NMC', '-dpng', ['-r', num2str(dpi)]);
    figure(3);
    print('C_e 1C CCCV NMC', '-dpng', ['-r', num2str(dpi)]);
    figure(4);
    print('C_se 1C CCCV NMC', '-dpng', ['-r', num2str(dpi)]);
    figure(5);
    print('C_ave 1C CCCV NMC', '-dpng', ['-r', num2str(dpi)]);
end 
