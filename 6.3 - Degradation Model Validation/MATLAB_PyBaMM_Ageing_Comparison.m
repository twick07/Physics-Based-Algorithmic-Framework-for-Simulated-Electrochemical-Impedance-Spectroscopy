%% Script to conduct one CCCV cycle
clear;
close all;


%% Model Selection

param  =1; %% param = 1 for LCO, param = 2 for NMC

%% Choose Degredation Models

SEI = 2;
li_pl = 0;
coupling = 2;
Ageing_models = [SEI li_pl coupling];

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

if SEI == 1
    
    Input_Data = xlsread("Pybamm_custom_CCCV_current_KL_ageing.csv");

end

if SEI == 2
    
    Input_Data = xlsread("Pybamm_custom_CCCV_current_DL_ageing.csv");

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



%% Import Drive Cycle Data

if SEI == 1
    
    Voltage_data_pybamm = xlsread("Pybamm_custom_CCCV_voltage_KL_ageing.csv");
    Current_data_pybamm = xlsread("Pybamm_custom_CCCV_current_KL_ageing.csv");
    j_SEI_data_pybamm = xlsread("Pybamm_custom_CCCV_j_SEI_KL_ageing.csv");
    c_SEI_data_pybamm = xlsread("Pybamm_custom_CCCV_c_SEI_KL_ageing.csv");
    L_SEI_data_pybamm = xlsread("Pybamm_custom_CCCV_L_SEI_KL_ageing.csv");
    R_SEI_data_pybamm = xlsread("Pybamm_custom_CCCV_R_SEI_KL_ageing.csv");
    eps_s_data_pybamm = xlsread("Pybamm_custom_CCCV_eps_e_n_KL_ageing.csv");


end

if SEI == 2
    
    Voltage_data_pybamm = xlsread("Pybamm_custom_CCCV_voltage_DL_ageing.csv");
    Current_data_pybamm = xlsread("Pybamm_custom_CCCV_current_DL_ageing.csv");
    j_SEI_data_pybamm = xlsread("Pybamm_custom_CCCV_j_SEI_DL_ageing.csv");
    c_SEI_data_pybamm = xlsread("Pybamm_custom_CCCV_c_SEI_DL_ageing.csv");
    L_SEI_data_pybamm = xlsread("Pybamm_custom_CCCV_L_SEI_DL_ageing.csv");
    R_SEI_data_pybamm = xlsread("Pybamm_custom_CCCV_R_SEI_DL_ageing.csv");
    eps_s_data_pybamm = xlsread("Pybamm_custom_CCCV_eps_e_n_DL_ageing.csv");

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
plot(time_vector_total,-j_SEI_tot(:,1)','r');
plot(time_vector_total,-j_SEI_tot(:,N_neg)','b');
plot(j_SEI_data_pybamm(:,2)',-j_SEI_data_pybamm(:,3)','r--');
plot(j_SEI_data_pybamm(:,2)',-j_SEI_data_pybamm(:,4)','b--');
xlabel('Time (s)');
ylabel('SEI Growth Flux (mols/m^2 s)');
legend('MATLAB SEI Growth CC/Neg','MATLAB SEI Growth Neg/Sep','PyBaMM SEI Growth CC/Neg','PyBaMM SEI Growth Neg/Sep',Location='northeast');
grid on;
fontsize(figure(3),'increase');


figure(4);
hold on;
plot(time_vector_total,c_SEI_tot(:,1)','r');
plot(time_vector_total,c_SEI_tot(:,N_neg)','b');
plot(c_SEI_data_pybamm(:,2)',c_SEI_data_pybamm(:,3)','r--');
plot(c_SEI_data_pybamm(:,2)',c_SEI_data_pybamm(:,4)','b--');
grid on;
xlabel('Time (s)');
ylabel('SEI Concentration (mols/m^3)');
legend('MATLAB SEI Concentration CC/Neg','MATLAB SEI Concentration Neg/Sep','PyBaMM SEI Concentration CC/Neg','PyBaMM SEI Concentration Neg/Sep',Location='northwest');
fontsize(figure(4),'increase');


figure(5);
hold on;
plot(time_vector_total,L_SEI_tot(:,1)','r');
plot(time_vector_total,L_SEI_tot(:,N_neg)','b');
plot(L_SEI_data_pybamm(:,2)',L_SEI_data_pybamm(:,3)','r--');
plot(L_SEI_data_pybamm(:,2)',L_SEI_data_pybamm(:,4)','b--');
grid on;
xlabel('Time (s)');
ylabel('SEI Thickness (m)');
legend('MATLAB SEI Thickness CC/Neg','MATLAB SEI Thickness Neg/Sep','PyBaMM SEI Thickness CC/Neg','PyBaMM SEI Thickness Neg/Sep',Location='northwest');
fontsize(figure(5),'increase');

figure(6);
hold on;
plot(time_vector_total,R_SEI_tot(:,1)','r');
plot(time_vector_total,R_SEI_tot(:,N_neg)','b');
plot(R_SEI_data_pybamm(:,2)',R_SEI_data_pybamm(:,3)','r--');
plot(R_SEI_data_pybamm(:,2)',R_SEI_data_pybamm(:,4)','b--');
grid on;
xlabel('Time (s)');
ylabel('SEI Resistance (Ohm m^2)');
legend('MATLAB SEI Resistance CC/Neg','MATLAB SEI Resistance Neg/Sep','PyBaMM SEI Resistance CC/Neg','PyBaMM SEI Resistance Neg/Sep',Location='northwest');
fontsize(figure(6),'increase');


figure(7);
hold on;
plot(time_vector_total,eps_sol_store_tot(:,1)','r');
plot(time_vector_total,eps_sol_store_tot(:,N_neg)','b');
plot(eps_s_data_pybamm(:,2)',eps_s_data_pybamm(:,3)','r--');
plot(eps_s_data_pybamm(:,2)',eps_s_data_pybamm(:,4)','b--');
grid on;
xlabel('Time (s)');
ylabel('Anode Porosity (-)');
legend('MATLAB SEI Porosity CC/Neg','MATLAB SEI Porosity Neg/Sep','PyBaMM SEI Porosity CC/Neg','PyBaMM SEI Porosity Neg/Sep',Location='northeast');
fontsize(figure(7),'increase');



%% Calculate Errors
RMSE_Current_percentage = mean(abs(I_app_total(1:end-1)-Current_data_pybamm(:,3)'));

RMSE_Voltage_percentage = mean(abs(V_cell_total(1:end-1)-Voltage_data_pybamm(:,3)')./Voltage_data_pybamm(:,3)');
RMSE_Voltage_abs = mean(abs(V_cell_total(1:end-1)-Voltage_data_pybamm(:,3)'));

fprintf('The Precentage Voltage RMSE = %.5f\n', RMSE_Voltage_percentage);

temp = abs((j_SEI_tot(1:end-1,1)' - j_SEI_data_pybamm(:,3)'));
idx = find(temp~=0);
RMSE_j_SEI_cc_neg = mean(temp(idx)/abs(j_SEI_data_pybamm(idx,3)'));

temp = abs((j_SEI_tot(1:end-1,N_neg)' - j_SEI_data_pybamm(:,4)'));
idx = find(temp~=0);
RMSE_j_SEI_neg_sep = mean(temp(idx)/abs(j_SEI_data_pybamm(idx,4)'));

RMSE_j_SEI = mean([RMSE_j_SEI_cc_neg RMSE_j_SEI_neg_sep]);

fprintf('The Precentage j_SEI RMSE = %.5f\n', RMSE_j_SEI);

RMSE_c_SEI_cc_neg = mean(abs((c_SEI_tot(1:end-1,1)' - c_SEI_data_pybamm(:,3)')./c_SEI_data_pybamm(:,3)'));
RMSE_c_SEI_neg_sep = mean(abs((c_SEI_tot(1:end-1,N_neg)' - c_SEI_data_pybamm(:,4)')./c_SEI_data_pybamm(:,4)'));

RMSE_c_SEI = mean([RMSE_c_SEI_cc_neg RMSE_c_SEI_neg_sep]);

fprintf('The Precentage c_SEI RMSE = %.5f\n', RMSE_c_SEI);

RMSE_L_SEI_cc_neg = mean(abs((L_SEI_tot(1:end-1,1)' - L_SEI_data_pybamm(:,3)')./L_SEI_data_pybamm(:,3)'));
RMSE_L_SEI_neg_sep = mean(abs((L_SEI_tot(1:end-1,N_neg)' - L_SEI_data_pybamm(:,4)')./L_SEI_data_pybamm(:,4)'));

RMSE_L_SEI = mean([RMSE_L_SEI_cc_neg RMSE_L_SEI_neg_sep]);

fprintf('The Precentage L_SEI RMSE = %.5f\n', RMSE_L_SEI);

RMSE_R_SEI_cc_neg = mean(abs((R_SEI_tot(1:end-1,1)' - R_SEI_data_pybamm(:,3)')./R_SEI_data_pybamm(:,3)'));
RMSE_R_SEI_neg_sep = mean(abs((R_SEI_tot(1:end-1,N_neg)' - R_SEI_data_pybamm(:,4)')./R_SEI_data_pybamm(:,4)'));

RMSE_R_SEI = mean([RMSE_R_SEI_cc_neg RMSE_R_SEI_neg_sep]);

fprintf('The Precentage R_SEI RMSE = %.5f\n', RMSE_R_SEI);


RMSE_eps_sol_store_cc_neg = mean(abs((eps_sol_store_tot(1:end-1,1)' - eps_s_data_pybamm(:,3)')./eps_s_data_pybamm(:,3)'));
RMSE_eps_sol_store_neg_sep = mean(abs((eps_sol_store_tot(1:end-1,N_neg)' - eps_s_data_pybamm(:,4)')./eps_s_data_pybamm(:,4)'));

RMSE_eps_s = mean([RMSE_eps_sol_store_cc_neg RMSE_eps_sol_store_neg_sep]);

fprintf('The Precentage eps_s RMSE = %.5f\n', RMSE_eps_s);



%% Export Figures
dpi = 300;
if SEI == 1
    figure(1);
    print('I_app_4_cycle_KL', '-dpng', ['-r', num2str(dpi)]);
    figure(2);
    print('V_cell_4_cycle_KL', '-dpng', ['-r', num2str(dpi)]);
    figure(3);
    print('SEI_Flux_KL_4_cycle', '-dpng', ['-r', num2str(dpi)]);
    figure(4);
    print('SEI_Conc_KL_4_cycle', '-dpng', ['-r', num2str(dpi)]);
    figure(5);
    print('SEI_Thick_KL_4_cycle', '-dpng', ['-r', num2str(dpi)]);
    figure(6);
    print('SEI_Res_KL_4_cycle', '-dpng', ['-r', num2str(dpi)]);
    figure(7);
    print('eps_e_KL_4_cycle', '-dpng', ['-r', num2str(dpi)]);
end 


if SEI == 2
    figure(1);
    print('I_app_4_cycle_DL', '-dpng', ['-r', num2str(dpi)]);
    figure(2);
    print('V_cell_4_cycle_DL', '-dpng', ['-r', num2str(dpi)]);
    figure(3);
    print('SEI_Flux_DL_4_cycle', '-dpng', ['-r', num2str(dpi)]);
    figure(4);
    print('SEI_Conc_DL_4_cycle', '-dpng', ['-r', num2str(dpi)]);
    figure(5);
    print('SEI_Thick_DL_4_cycle', '-dpng', ['-r', num2str(dpi)]);
    figure(6);
    print('SEI_Res_DL_4_cycle', '-dpng', ['-r', num2str(dpi)]);
    figure(7);
    print('eps_e_DL_4_cycle', '-dpng', ['-r', num2str(dpi)]);
end 
