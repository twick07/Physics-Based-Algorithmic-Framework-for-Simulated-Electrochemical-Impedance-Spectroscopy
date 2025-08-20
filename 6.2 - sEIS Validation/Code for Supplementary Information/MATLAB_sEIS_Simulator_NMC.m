clear;
close all;
%% Script to time domain model EIS measurements and compare with PyBaMM
SoC = 100;

%% Import PyBaMM Data
Impedance_data = xlsread("NMC_EIS_Impedance_"+int2str(SoC)+".csv");
Time_data = xlsread("NMC_EIS_time_vectors_"+int2str(SoC)+".csv");
Current_data = xlsread("NMC_EIS_current_vectors_"+int2str(SoC)+".csv");
Voltage_data = xlsread("NMC_EIS_voltage_vectors_"+int2str(SoC)+".csv");

%% MDN Modelling Parameters

Comparison_points = 1000;
time_step = 1;
order = 5;
N_tot = 50;
N_shell = 10;
V_cut_off_low = 2.8;
V_cut_off_high = 4.0;
n_parameters = 42;
Cut_off_capacity = 0.1;
%%% Model Selection
%param  =1; %% param = 1 for LCO, param = 2 for NMC, 
param  =2;
C_rate = 1;
sim_type = 4;

%% Import Parameters
%Battery_Parameters_P2D_5;
Battery_Parameters_P2D_EIS;
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

%%% Plating Parameters
Model_Parameters(37) = k_li;
Model_Parameters(38) = alpha_pl_neg;
Model_Parameters(39) = alpha_pl_pos;
Model_Parameters(40) = kappa_li;
Model_Parameters(41) = M_Li;
Model_Parameters(42) = rho_Li;


%% Choose Degredation Models
SEI = 0;
li_pl = 0;
coupling = 0;
Ageing_models = [SEI li_pl coupling];

%%% Initialisation for MATLAB Code
C_e_init = zeros(1,N_tot);
C_s_full_init = zeros(1,N_shell,N_tot);
C_se_init = zeros(1,N_tot);
I_app_init = zeros(1,N_tot);
c_SEI_init = zeros(1,N_tot);
c_li_init = zeros(1,N_tot);
psi_s_init = zeros(1,N_tot);
psi_e_init = zeros(1,N_tot);
j_tot_init = zeros(1,N_tot);

%% EIS Parameters
N_frequencies = length(Impedance_data(:,2));
Points_Per_cycle = 56;
No_of_samples = 56;
Cycles_per_freq = 10;
Total_cycles = 10;
Output_cycles = 5;
Out_cycle = Total_cycles - Output_cycles;
last_x_cycles = 5;
time_step = 1;
No_of_frequencies = 40;
Frequency_Vector = Impedance_data(:,2)';
Delta_I = 0.2;
Impedance_calculated = zeros(N_frequencies,1);
Impedance_pybamm = zeros(N_frequencies,1);
V_rmse = zeros(N_frequencies,1);
Impedance_rmse = zeros(N_frequencies,1);
Impedance_rmse_real = zeros(N_frequencies,1);
Impedance_rmse_imag = zeros(N_frequencies,1);
break_con = 0;
%idx_init = 1;
idx_init = 1;
Freq_turn = 100;
output_voltage_current_plots = [2 21 34];

%% EIS Simulation
for i = idx_init:length(Frequency_Vector)
    
    

    %% Make the EIS Current
    if Frequency_Vector(i) < 20e-3
        [time_vector, I_app_EIS,idx_last_cycle] = Single_EIS_Current_Generator(Delta_I,Frequency_Vector(i),time_step,Cycles_per_freq,last_x_cycles,Points_Per_cycle);
        I_app_init = I_app_EIS; 
        time_step_vector = time_vector; 

    else

        I_app_init = Current_data(:,i+1)';
        time_step_vector = Time_data(:,i+1)';
        idx_last_cycle = No_of_samples*Out_cycle+1;

    end
    %% Generate EIS Response
    
    
    if Frequency_Vector(i) < Freq_turn
        
        [V_cell_EIS,Flag_convergence,time_vector_EIS,C_e,C_s_full,C_se,I_app_EIS_singles] = P2D_function_general(Model_Parameters,N_tot,N_shell,order,param,time_step_vector,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,C_e_init,C_s_full_init,C_se_init,I_app_init,C_rate);

    else
        [V_cell_EIS,Flag_convergence,time_vector_EIS,C_e,C_s_full,C_se,c_SEI,c_li_pl,I_app_EIS_singles,psi_s,psi_e,j_tot] = P2D_MDN_function_ODE_int(Model_Parameters,N_tot,N_shell,order,param,time_step_vector,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init);

    end
  
    
    if and(sum(i == output_voltage_current_plots) == 1, SoC == 100) 
        figure(i);
        hold on;
        yyaxis left;
        plot(time_vector_EIS(1:end),I_app_EIS_singles(1:end),'r');
        plot(Time_data(:,i+1),Current_data(:,i+1),'r*');
        xlabel('Time (s)');
        ylabel('Current [A/m^2]');
        yyaxis right;
        plot(time_vector_EIS(1:end),V_cell_EIS(1:end),'b');
        plot(Time_data(:,i+1),Voltage_data(:,i+1),'b*');
        ylabel('Voltage [V]');
        title('Frequency = ',num2str(Frequency_Vector(i)));
        grid on;
        legend('Applied Current MATLAB','Applied Current Pybamm','Cell Voltage MATLAB','Cell Voltage Pybamm');
        fontsize(figure(i),'increase');

    end
    if Flag_convergence ~= 0

        output = ['Break at Freq = ', num2str(Frequency_Vector(i))];
        disp(output)
        break_con = 1;

        break;



    end

    %% Calculate EIS From MDN Solver 

    time_vector_SS = time_vector_EIS(idx_last_cycle:end);
    I_app_EIS_SS = I_app_EIS_singles(idx_last_cycle:end);
    V_cell_EIS_SS = V_cell_EIS(idx_last_cycle:end);

    % FFT
    current_fft = fft(I_app_EIS_SS);
    voltage_fft = fft(V_cell_EIS_SS);

    % Get index of first harmonic
    [~, idx] = max(abs(current_fft));

    Impedance_calculated(i,:) = -voltage_fft(idx) / current_fft(idx);

    %% Caclulate EIS from pybamm
    
    time_ss_pybamm = Time_data(No_of_samples*Out_cycle+1:end,i+1);

    current_ss_pybamm = Current_data(No_of_samples*Out_cycle+1:end,i+1);
     
    voltage_ss_pybamm = Voltage_data(No_of_samples*Out_cycle+1:end,i+1);

    

    %% Caclulate fft from pybamm normalised to same points
    V_cell_normalised = pchip(time_ss_pybamm,voltage_ss_pybamm,time_vector_SS);
    I_app_normalised = pchip(time_ss_pybamm,current_ss_pybamm,time_vector_SS);

    current_fft_2 = fft(I_app_normalised);
    voltage_fft_2 = fft(V_cell_normalised);

    % Get index of first harmonic
    [~, idx_2] = max(abs(current_fft_2));

    Impedance_pybamm(i,:) = -voltage_fft_2(idx_2) / current_fft_2(idx_2);

   
    %% Calculate RMSE
    
    Impedance_rmse_real(i,:) = abs((real(Impedance_calculated(i,:))-real(Impedance_pybamm(i,:))) ./real(Impedance_pybamm(i,:)));
    Impedance_rmse_imag(i,:) = abs((imag(Impedance_calculated(i,:))-imag(Impedance_pybamm(i,:))) ./imag(Impedance_pybamm(i,:)));
    Impedance_rmse(i,:) = mean([Impedance_rmse_real(i,:) Impedance_rmse_imag(i,:)]);
    V_rmse(i,:) = mean(abs(V_cell_EIS_SS - V_cell_normalised)./V_cell_normalised);
    
     
    %% Output Info
    if sum(i == output_voltage_current_plots) == 1
        output_Reigon_Errors = ['Frequency = ',num2str(Frequency_Vector(i)),', Voltage Error = ',num2str(V_rmse(i,:)), ', Impedance Error = ',num2str(Impedance_rmse(i,:)),...
            ', Real Impedance Error = ',num2str(Impedance_rmse_real(i,:)),', Imag Impedance Error = ',num2str(Impedance_rmse_imag(i,:))];
        disp(output_Reigon_Errors);
    end
end

%% Post Processing
if break_con == 1
    idx_plot = i-1;

else
    idx_plot = i;
    
end


%% Plot Impedance Data

multi = 1e3;
figure(i+1);
grid on;
hold on
plot(multi*abs(real((Impedance_calculated(idx_init:idx_plot,:)))), multi*abs(imag((Impedance_calculated(idx_init:idx_plot,:)))),'b')
plot(multi*abs(real(Impedance_pybamm(idx_init:idx_plot,:))),multi*abs(imag(Impedance_pybamm(idx_init:idx_plot,:))),'r');

xlabel('Z_{Re} (mOhms)');
ylabel('-Z_{Im} (mOhms)');

legend('MATLAB Solver','PyBaMM Solver',location='northwest')
title("sEIS Results for SoC = "+int2str(SoC)+"%");
fontsize(figure(i+1),'increase');


%% Calculate Error

Mean_Impedance_Error = mean(100*Impedance_rmse(idx_init:idx_plot,:));

fprintf('The Precentage Impedance RMSE = %.5f\n', Mean_Impedance_Error);

%% Export Data
dpi = 300;
if SoC == 100

    figure(output_voltage_current_plots(1));
    print('sEIS_Voltage_Current_f_1mHz_NMC', '-dpng', ['-r', num2str(dpi)]);
    figure(output_voltage_current_plots(2));
    print('sEIS_Voltage_Current_f_1Hz_NMC', '-dpng', ['-r', num2str(dpi)]);
    figure(output_voltage_current_plots(3));
    print('sEIS_Voltage_Current_f_100Hz_NMC', '-dpng', ['-r', num2str(dpi)]);

end

figure(i+1);
print("sEIS SoC "+int2str(SoC)+" NMC", '-dpng', ['-r', num2str(dpi)]);

Export_matrix = [Frequency_Vector(idx_init:idx_plot)'  real((Impedance_calculated(idx_init:idx_plot,:)))  imag((Impedance_calculated(idx_init:idx_plot,:)))];
Export_matrix = [Frequency_Vector(idx_init:idx_plot)'  real((Impedance_calculated(idx_init:idx_plot,:)))  imag((Impedance_calculated(idx_init:idx_plot,:)))];
Export_matrix_pybamm = [Frequency_Vector(idx_init:idx_plot)' real(Impedance_pybamm(idx_init:idx_plot,:)) imag(Impedance_pybamm(idx_init:idx_plot,:))];
csvwrite("Matlab_NMC_EIS_"+int2str(SoC)+"_SoC.csv",Export_matrix);
csvwrite("PyBaMM_NMC_EIS_"+int2str(SoC)+"_SoC.csv",Export_matrix_pybamm);

