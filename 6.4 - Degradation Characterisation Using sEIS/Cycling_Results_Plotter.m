%% Script to Generate Cycling Results Plots
clear;
close all;

%% Import Data
%%% Import KL SEI Growth Cycling Results
Input_Data = xlsread("Cycling_Data_KL_50_Cyles.csv");
time_vector_KL = Input_Data(:,1);
I_app_KL = Input_Data(:,2);
V_cell_total_KL = Input_Data(:,3);

Input_Data = xlsread("50_cycle_KL_ageing_cap_lost.csv");

Q_lost_KL = Input_Data(2:end,2);

Input_Data = xlsread("50_cycle_KL_ageing_SoH_lost.csv");

SoH_vector_KL = Input_Data(2:end,2);

%%% Import DL SEI Growth Cycling Results
Input_Data = xlsread("Cycling_Data_DL_50_Cyles.csv");
time_vector_DL = Input_Data(:,1);
I_app_DL = Input_Data(:,2);
V_cell_total_DL = Input_Data(:,3);

Input_Data = xlsread("50_cycle_DL_ageing_cap_lost.csv");

Q_lost_DL = Input_Data(2:end,2);

Input_Data = xlsread("50_cycle_DL_ageing_SoH_lost.csv");

SoH_vector_DL = Input_Data(2:end,2);


%% Plotting the Results
%% Plot Data
N_cycles_KL = length(Q_lost_KL);
N_cycles_DL = length(Q_lost_DL);

figure(1);
hold on;
plot(time_vector_DL,I_app_DL,'b-');
plot(time_vector_KL,I_app_KL,'r--');
xlabel('Time (s)');
ylabel('Cell Current (A/m^2)');
legend('Diffusion-Limited SEI Growth','Kinetic-Limited SEI Growth');
grid on;
fontsize(figure(1),'increase')
fontsize(figure(1),'increase')


figure(2);
hold on;
plot(time_vector_DL,V_cell_total_DL,'b-');
plot(time_vector_KL,V_cell_total_KL,'r--');
xlabel('Time (s)');
ylabel('Cell Voltage (V)');
legend('Diffusion-Limited SEI Growth','Kinetic-Limited SEI Growth');
grid on;
fontsize(figure(2),'increase')
fontsize(figure(2),'increase')


figure(3);
hold on;
plot(linspace(1,N_cycles_DL,N_cycles_DL),Q_lost_DL,'r-');
plot(linspace(1,N_cycles_KL,N_cycles_KL),Q_lost_KL,'b-');
xlabel('Cycles [-]');
ylabel('Capacity Lost [%]');
legend('Diffusion-Limited SEI Growth','Kinetic-Limited SEI Growth');
grid on;
fontsize(figure(3),'increase')
fontsize(figure(3),'increase')




figure(4);
hold on;
plot(linspace(1,N_cycles_DL,N_cycles_DL),SoH_vector_DL,'r-');
plot(linspace(1,N_cycles_KL,N_cycles_KL),SoH_vector_KL,'b-');
xlabel('Cycles [-]');
ylabel('State of Health [%]');
legend('Diffusion-Limited SEI Growth','Kinetic-Limited SEI Growth');
grid on;
fontsize(figure(4),'increase')
fontsize(figure(4),'increase')



%% Print Figures
dpi = 300;

figure(1);
print('50_cycle_I_app', '-dpng', ['-r', num2str(dpi)]);
figure(2);
print('50_cycle_V_cell', '-dpng', ['-r', num2str(dpi)]);
figure(3);
print('Cap_Lost_Comp', '-dpng', ['-r', num2str(dpi)]);
figure(4);
print('SoH_Lost_Comp', '-dpng', ['-r', num2str(dpi)]);




