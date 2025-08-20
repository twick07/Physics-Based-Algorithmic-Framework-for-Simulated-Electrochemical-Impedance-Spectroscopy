%% Script to conduct one CCCV cycle
clear;
close all;


%% Model Selection

param  =1; %% param = 1 for LCO/C6, param = 2 for NMC/C6
N_cycles = 4; 

%% Choose Degredation Models
SEI = 1;
li_pl = 0;
coupling = 2;
Ageing_models = [SEI li_pl coupling];

%% P2D Modelling Parameters
Comparison_points = 1000;
time_step = 1;
order = 5;
%N_tot = 30;
N_tot = 40;
N_shell = 5;



% Cut Off setting for LCO
if param == 1
    V_cut_off_low = 3.4;
    V_cut_off_high_interim = 3.97;
    V_cut_off_high_80 = 3.9;
    V_cut_off_high_60 = 3.84;
    V_cut_off_high_40 = 3.81;
    V_cut_off_high_end = V_cut_off_high_interim;
    V_cut_off_high = V_cut_off_high_interim;
end

if param == 2
    % Cut off setting for NMC
    V_cut_off_low = 3.4;
    V_cut_off_high_100 = 3.97;
    V_cut_off_high_80 = 3.848;
    V_cut_off_high_60 = 3.776;
    V_cut_off_high_40 = 3.712;
    V_cut_off_high_interim = V_cut_off_high_100;
    V_cut_off_high_end = V_cut_off_high_100;
    V_cut_off_high = V_cut_off_high_100;

end


%V_cut_off_test_vectors = [V_cut_off_high_60  V_cut_off_high_interim];
V_cut_off_test_vectors = [V_cut_off_high_interim];
%I_app_cut_off = -0.1; %A/m^2
I_app_cut_off_percent = 0.05; % Percentage of the 1C current
Cut_off_capacity = 0.1;
tol_CV = 5e-4;
test_freq = 1;
Test_cycles = [51]; 
n_parameters = 46;
C_rate_non_test = 1;
C_rate_test = 1;





%% Simulation Setting
%sim_type = 2; %% sim_type = 1 for discharge, sim_type =2 for charge, sim_type =3 discharge/charge
C_rate = 1;
t_rest_discharge = 5*60; %%% Rest time after discharge
t_rest_charge = 5*60; %%% Rest time after CV charge

%% Import Parameters
Battery_Parameters_P2D_MATLAB_2
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
N_tot_store = 40;
soc_points = length(V_cut_off_test_vectors)*(length(Test_cycles));
C_e_EIS = zeros(soc_points,N_tot_store);
C_s_full_EIS = zeros(soc_points,N_shell,N_tot_store);
C_se_EIS = zeros(soc_points,N_tot_store);
c_SEI_EIS = zeros(soc_points,N_tot_store);
c_li_EIS = zeros(soc_points,N_tot_store);
psi_s_EIS = zeros(soc_points,N_tot_store);
psi_e_EIS = zeros(soc_points,N_tot_store);
j_tot_EIS = zeros(soc_points,N_tot_store);
eps_sol_EIS = zeros(soc_points,N_tot_store);

%% Capacity Data
q_lost = zeros(1,N_cycles);


%% Simulations
Test_cycle_count = 0;
cycle_end_points = zeros(1,N_cycles);

%% CCCV Cycling Simulation
for k =1:N_cycles

    %% CC Discharge Simulation
    tic;
    sim_type = 1;
    if k ==1
        C_e_init = zeros(1,N_tot);
        C_s_full_init = zeros(1,N_shell,N_tot);
        C_se_init = zeros(1,N_tot);
        I_app_init = zeros(1,N_tot);
        c_SEI_init = zeros(1,N_tot);
        c_li_init = zeros(1,N_tot);
        psi_s_init = zeros(1,N_tot);
        psi_e_init = zeros(1,N_tot);
        j_tot_init = zeros(1,N_tot);
        eps_sol_init = zeros(1,N_tot);
        

    else
        C_e_init = C_e_end;
        C_s_full_init = C_s_full_end;
        C_se_init = C_se_end;
        I_app_init = zeros(1,N_tot);
        c_SEI_init = c_SEI_end;    
        c_li_init = c_li_end;
        j_tot_init = j_tot_end;
        psi_s_init = psi_s_end;
        psi_e_init = psi_e_end;
        eps_sol_init = eps_sol_end;


    end
    
   
    
    
    
    [V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,c_SEI,c_li,I_app,psi_s,psi_e,j_tot,eps_sol,L_SEI,R_SEI,j_SEI,eps_sol_store] = P2D_function_general(Model_Parameters,N_tot,N_shell,order,param,time_step,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);

    
    %%% Store Data
    time_vector_CC_discharge = time_vector;
    V_cell_CC_discharge = V_cell;
    I_app_CC_discharge = I_app;
    
    C_e_end = C_e(end,:);
    C_s_full_end = C_s_full(end,:,:);
    C_se_end = C_se(end,:);
    c_SEI_end = c_SEI(end,:);
    c_li_end = c_li(end,:);
    j_tot_end = j_tot(end,:);
    psi_s_end = psi_s(end-1:end,:);
    psi_e_end = psi_e(end-1:end,:);
    eps_sol_end = eps_sol;
    
    %%% Create Data Containers for ageing
    if k == 1
        c_SEI_store = zeros(N_cycles,length(c_SEI_end));
        c_li_store = zeros(N_cycles,length(c_li_end));
    end

    %%% Update Global Data
    if k == 1
        time_vector_total = [time_vector];
        I_app_total =[I_app];
        V_cell_total = [V_cell];
        C_e_total = [C_e];
        C_se_total = [C_se];
        C_s_ave_total = [squeeze(mean(C_s_full,2))];
        %L_SEI_total = [L_SEI];
    else
        time_vector_total = [time_vector_total time_vector_total(end)+time_vector+time_step];
        I_app_total =[I_app_total I_app];
        V_cell_total = [V_cell_total V_cell];
        C_e_total = [C_e_total; C_e];
        C_se_total = [C_se_total; C_se];
        C_s_ave_total = [C_s_ave_total;squeeze(mean(C_s_full,2))];
        %L_SEI_total =[L_SEI_total L_SEI];
    end


    %% Post Discharge Rest
    sim_type = 0;
    I_app_rest = zeros(1,round(t_rest_discharge./time_step));
    C_e_init = C_e_end;
    C_s_full_init = C_s_full_end;
    C_se_init = C_se_end;
    I_app_init = I_app_rest;
    c_SEI_init = c_SEI_end;
    c_li_init = c_li_end;
    j_tot_init = j_tot_end;
    psi_s_init = psi_s_end;
    psi_e_init = psi_e_end;
    eps_sol_init = eps_sol_end;
    
    
    
  
    
    [V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,c_SEI,c_li,I_app,psi_s,psi_e,j_tot,eps_sol,L_SEI,R_SEI,j_SEI,eps_sol_store] = P2D_function_general(Model_Parameters,N_tot,N_shell,order,param,time_step,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);

    
    
    %%% Update local variables
    C_e_end = C_e(end,:);
    C_s_full_end = C_s_full(end,:,:);
    C_se_end = C_se(end,:);
    c_SEI_end = c_SEI(end,:);
    c_li_end = c_li(end,:);
    j_tot_end = j_tot(end,:);
    psi_s_end = psi_s(end,:);
    psi_e_end = psi_e(end,:);
    eps_sol_end = eps_sol;
    %%% Update Global Data
    time_vector_total = [time_vector_total time_vector_total(end)+time_vector+time_step];
    I_app_total =[I_app_total I_app];
    V_cell_total = [V_cell_total V_cell];
    C_e_total = [C_e_total; C_e];
    C_se_total = [C_se_total; C_se];
    C_s_ave_total = [C_s_ave_total;squeeze(mean(C_s_full,2))];
    


    
    %% CC & CV Charge Simulation
    if sum(k == Test_cycles) == 0
        %% CC Charge Simulation
        C_rate = C_rate_non_test;
        V_cut_off_high = V_cut_off_high_interim;
    
        sim_type = 2;
        C_e_init = C_e_end;
        C_s_full_init = C_s_full_end;
        C_se_init = C_se_end;
        I_app_init = I_app;
        c_SEI_init = c_SEI_end;
        c_li_init = c_li_end;
        j_tot_init = j_tot_end;
        psi_s_init = psi_s_end;
        psi_e_init = psi_e_end;
        eps_sol_init = eps_sol_end;
        
       
       
       
        [V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,c_SEI,c_li,I_app,psi_s,psi_e,j_tot,eps_sol,L_SEI,R_SEI,j_SEI,eps_sol_store] = P2D_function_general(Model_Parameters,N_tot,N_shell,order,param,time_step,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);

        
        %%% Store Data
        time_vector_CC_charge = time_vector;
        V_cell_CC_charge = V_cell;
        I_app_CC_charge = I_app;
        
        %%% Update local variables
        C_e_end = C_e(end,:);
        C_s_full_end = C_s_full(end,:,:);
        C_se_end = C_se(end,:);
        c_SEI_end = c_SEI(end,:);
        c_li_end = c_li(end,:);
        j_tot_end = j_tot(end,:);
        psi_s_end = psi_s(end,:);
        psi_e_end = psi_e(end,:);
        eps_sol_end = eps_sol;

        %%% Output the current psi_delta at end of charge
        output_Reigon_Errors = ['Psi_Delta at end of charge = ',num2str(psi_s_end - psi_e_end)];
        disp(output_Reigon_Errors);
        
        %%% Update Global Data
        time_vector_total = [time_vector_total time_vector_total(end)+time_vector+time_step];
        I_app_total =[I_app_total I_app];
        V_cell_total = [V_cell_total V_cell];
        C_e_total = [C_e_total; C_e];
        C_se_total = [C_se_total; C_se];
        C_s_ave_total = [C_s_ave_total;squeeze(mean(C_s_full,2))];
        
            
        %% CV charge simulation
        %%% Calculate CV Charge Current   
        I_app_cut_off = I_app_cut_off_percent*I_app(end);
        V_cell_CV = V_cell_total(end);
        I_app_CV = I_app(end);
        j=1;
            
        while abs(I_app_CV(end))>abs(I_app_cut_off)
            
            %%% Update local variables
            C_e_init = C_e(end,:);
            C_s_full_init = C_s_full(end,:,:);
            C_se_init = C_se(end,:);
            c_SEI_init = c_SEI(end,:);
            c_li_init = c_li(end,:);
            j_tot_init = j_tot(end,:);
            psi_s_init = psi_s(end,:);
            psi_e_init = psi_e(end,:);
            eps_sol_init = eps_sol;

            
            
            time_end = 1*time_step; 
            err = 1;
            err = ones(3,1);
            tol_con = tol_CV;
           
            
            %%% Set I_app limits
            I_app_guesses = [1.1*I_app(end) I_app_cut_off mean([I_app(end) 0.9*I_app_cut_off])];
                
            while (abs(err(3))>tol_con)
            
            
                %%% Calculate Initial errors
                for i = 1:length(I_app_guesses)
                    I_app_init = I_app_guesses(i);
                   
                   
                    
                     [V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,c_SEI,c_li_pl,I_app,psi_s,psi_e,j_tot,eps_sol_store,eps_sol] = P2D_function_delta_t(Model_Parameters,N_tot,N_shell,order,param,V_cut_off_low,V_cut_off_high,Ageing_models,time_step,time_end,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);
   
                    err(i) = (V_cell(end) - V_cut_off_high);
                
                end
        
                if sign(err(1)*err(3)) == sign(err(2)*err(3))
            
                       
                        I_app_guesses(1) = I_app_guesses(1) - 1;%Changed from 0.1
                        I_app_guesses(2) = I_app_guesses(2) + 1;
                        I_app_guesses(3) = (I_app_guesses(1)+I_app_guesses(2))/2 ; 
        
                else
                
                    if err(1)*err(3) < 0 
                            
                        I_app_guesses(1) = I_app_guesses(1);
                    
                        I_app_guesses(2) = I_app_guesses(3);
                    
                        
                        I_app_guesses(3) = (I_app_guesses(1)+I_app_guesses(2))/2 ;   
                                    
                    end
                    
                    if err(2)*err(3) < 0
                    
                        I_app_guesses(1) = I_app_guesses(3);
                        I_app_guesses(2) = I_app_guesses(2);
                        
                        I_app_guesses(3) = (I_app_guesses(1)+I_app_guesses(2))/2 ;
                    end
                end
            
            end
                
               
            
            V_cell_CV = [V_cell_CV V_cell(2:end)];
            I_app_CV = [I_app_CV I_app(2:end)];
            j= j+1;
            
        end
            

        %end
            
        %% Forward Copute CV simulation
        sim_type = 0;
        I_app_init = I_app_CV;
        C_e_init = C_e_end;
        C_s_full_init = C_s_full_end;
        C_se_init = C_se_end;
        c_SEI_init = c_SEI_end;
        c_li_init = c_li_end;
        j_tot_init = j_tot_end;
        psi_s_init = psi_s_end;
        psi_e_init = psi_e_end;
        eps_sol_init = eps_sol_end;
    
        
        
       
        
        [V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,c_SEI,c_li,I_app,psi_s,psi_e,j_tot,eps_sol,L_SEI,R_SEI,j_SEI,eps_sol_store] = P2D_function_general(Model_Parameters,N_tot,N_shell,order,param,time_step,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);

        
        %%% Store Data
        time_vector_CC_discharge_rest = time_vector;
        V_cell_CC_discharge_rest = V_cell;
        I_app_CC_discharge_rest = I_app;
        
        %%% Update local variables
        C_e_end = C_e(end,:);
        C_s_full_end = C_s_full(end,:,:);
        C_se_end = C_se(end,:);
        c_SEI_end = c_SEI(end,:);
        c_li_end = c_li(end,:);
        j_tot_end = j_tot(end,:);
        psi_s_end = psi_s(end-1:end,:);
        psi_e_end = psi_e(end-1:end,:);
        eps_sol_end = eps_sol;
        
        
        %%% Update Global Data
        time_vector_total = [time_vector_total time_vector_total(end)+time_vector+time_step];
        I_app_total =[I_app_total I_app];
        V_cell_total = [V_cell_total V_cell];
        C_e_total = [C_e_total; C_e];
        C_se_total = [C_se_total; C_se];
        C_s_ave_total = [C_s_ave_total;squeeze(mean(C_s_full,2))];
    
        %% Post CV Rest
        sim_type = 0;
        I_app_rest = zeros(1,round(t_rest_discharge./time_step));
        C_e_init = C_e_end;
        C_s_full_init = C_s_full_end;
        C_se_init = C_se_end;
        I_app_init = I_app_rest;
        c_SEI_init = c_SEI_end;
        c_li_init = c_li_end;
        j_tot_init = j_tot_end;
        psi_s_init = psi_s_end;
        psi_e_init = psi_e_end;
        eps_sol_init = eps_sol_end;
    
    
        
        
        
       
       
        [V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,c_SEI,c_li,I_app,psi_s,psi_e,j_tot,eps_sol,L_SEI,R_SEI,j_SEI,eps_sol_store] = P2D_function_general(Model_Parameters,N_tot,N_shell,order,param,time_step,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);

        
        
        %%% Update local variables
        C_e_end = C_e(end,:);
        C_s_full_end = C_s_full(end,:,:);
        C_se_end = C_se(end,:);
        c_SEI_end = c_SEI(end,:);
        c_li_end = c_li(end,:);
        j_tot_end = j_tot(end,:);
        psi_s_end = psi_s(end,:);
        psi_e_end = psi_e(end,:);
        eps_sol_end = eps_sol;
    
        %%% Update Global Data
        time_vector_total = [time_vector_total time_vector_total(end)+time_vector+time_step];
        I_app_total =[I_app_total I_app];
        V_cell_total = [V_cell_total V_cell];
        C_e_total = [C_e_total; C_e];
        C_se_total = [C_se_total; C_se];
        C_s_ave_total = [C_s_ave_total;squeeze(mean(C_s_full,2))];
        %L_SEI_total =[L_SEI_total L_SEI];

    end
    if sum(k == Test_cycles) ~= 0
        for m = 1:length(V_cut_off_test_vectors)
            %% CC Charge Simulation
            C_rate = C_rate_test;
            V_cut_off_high = V_cut_off_test_vectors(m);
            sim_type = 2;
            C_e_init = C_e_end;
            C_s_full_init = C_s_full_end;
            C_se_init = C_se_end;
            I_app_init = I_app;
            c_SEI_init = c_SEI_end;
            c_li_init = c_li_end;
            j_tot_init = j_tot_end;
            psi_s_init = psi_s_end;
            psi_e_init = psi_e_end;
            eps_sol_init = eps_sol_end;

            
            [V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,c_SEI,c_li,I_app,psi_s,psi_e,j_tot,eps_sol,L_SEI,R_SEI,j_SEI,eps_sol_store] = P2D_function_general(Model_Parameters,N_tot,N_shell,order,param,time_step,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);

            %%% Update local variables
            C_e_end = C_e(end,:);
            C_s_full_end = C_s_full(end,:,:);
            C_se_end = C_se(end,:);
            c_SEI_end = c_SEI(end,:);
            c_li_end = c_li(end,:);
            j_tot_end = j_tot(end,:);
            psi_s_end = psi_s(end,:);
            psi_e_end = psi_e(end,:);
            eps_sol_init = eps_sol;

            %%% Output the current psi_delta at end of charge
            output_Reigon_Errors = ['Psi_Delta at end of charge = ',num2str(psi_s_end - psi_e_end)];
            disp(output_Reigon_Errors);
        
            %%% Update Global Data
            time_vector_total = [time_vector_total time_vector_total(end)+time_vector+time_step];
            I_app_total =[I_app_total I_app];
            V_cell_total = [V_cell_total V_cell];
            C_e_total = [C_e_total; C_e];
            C_se_total = [C_se_total; C_se];
            C_s_ave_total = [C_s_ave_total;squeeze(mean(C_s_full,2))];
            

            %% CV charge simulation
            %%% Calculate CV Charge Current
            %if k <= 2
            I_app_cut_off = I_app_cut_off_percent*I_app(end);
            %V_cell_CV = V_cell(end);
            V_cell_CV = V_cell_total(end);
            I_app_CV = I_app(end);
            j=1;
                
            while abs(I_app_CV(end))>abs(I_app_cut_off)
                %%% Initialise
                
                %%% Update local variables
                C_e_init = C_e(end,:);
                C_s_full_init = C_s_full(end,:,:);
                C_se_init = C_se(end,:);
                c_SEI_init = c_SEI(end,:);
                c_li_init = c_li(end,:);
                j_tot_init = j_tot(end,:);
                psi_s_init = psi_s(end,:);
                psi_e_init = psi_e(end,:);
                eps_sol_init = eps_sol;
    
                
                
                time_end = 1*time_step; 
                err = 1;
                err = ones(3,1);
                tol_con = tol_CV;
               
                
                %%% Set I_app limits
                I_app_guesses = [1.1*I_app(end) I_app_cut_off mean([I_app(end) 0.9*I_app_cut_off])];
                    
                while (abs(err(3))>tol_con)
                
                
                    %%% Calculate Initial errors
                    for i = 1:length(I_app_guesses)
                        I_app_init = I_app_guesses(i);
                       
                       
                        
                        [V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,c_SEI,c_li_pl,I_app,psi_s,psi_e,j_tot,eps_sol_store,eps_sol] = P2D_function_delta_t(Model_Parameters,N_tot,N_shell,order,param,V_cut_off_low,V_cut_off_high,Ageing_models,time_step,time_end,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);
        
                        err(i) = (V_cell(end) - V_cut_off_high);
                    
                    end
            
                    if sign(err(1)*err(3)) == sign(err(2)*err(3))
                
                            
                            I_app_guesses(1) = I_app_guesses(1) - 1;
                            I_app_guesses(2) = I_app_guesses(2) + 1;
                            I_app_guesses(3) = (I_app_guesses(1)+I_app_guesses(2))/2 ; 
            
                    else
                    
                        if err(1)*err(3) < 0 
                                
                            I_app_guesses(1) = I_app_guesses(1);
                        
                            I_app_guesses(2) = I_app_guesses(3);
                        
                            I_app_guesses(3) = (err(2)*I_app_guesses(1)-err(1)*I_app_guesses(2))/(err(2)-err(1)) ;
                                
                                        
                        end
                        
                        if err(2)*err(3) < 0
                        
                            I_app_guesses(1) = I_app_guesses(3);
                            I_app_guesses(2) = I_app_guesses(2);
                            I_app_guesses(3) = (err(2)*I_app_guesses(1)-err(1)*I_app_guesses(2))/(err(2)-err(1)) ;
                        
                        end
                    end
                
                end
                    
                   
                
                V_cell_CV = [V_cell_CV V_cell(2:end)];
                I_app_CV = [I_app_CV I_app(2:end)];
                j= j+1;
                
            end
                        
            
                    %end            
            %% Forward Copute CV simulation
            sim_type = 0;
            I_app_init = I_app_CV;
            C_e_init = C_e_end;
            C_s_full_init = C_s_full_end;
            C_se_init = C_se_end;
            c_SEI_init = c_SEI_end;
            c_li_init = c_li_end;
            j_tot_init = j_tot_end;
            psi_s_init = psi_s_end;
            psi_e_init = psi_e_end;
            eps_sol_init = eps_sol_end;
        
            
            
           
            
            [V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,c_SEI,c_li,I_app,psi_s,psi_e,j_tot,eps_sol,L_SEI,R_SEI,j_SEI,eps_sol_store] = P2D_function_general(Model_Parameters,N_tot,N_shell,order,param,time_step,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);
      
            
            %%% Store Data
            time_vector_CC_discharge_rest = time_vector;
            V_cell_CC_discharge_rest = V_cell;
            I_app_CC_discharge_rest = I_app;
            
            %%% Update local variables
            C_e_end = C_e(end,:);
            C_s_full_end = C_s_full(end,:,:);
            C_se_end = C_se(end,:);
            c_SEI_end = c_SEI(end,:);
            c_li_end = c_li(end,:);
            j_tot_end = j_tot(end,:);
            psi_s_end = psi_s(end-1:end,:);
            psi_e_end = psi_e(end-1:end,:);
            eps_sol_end = eps_sol;
            
            
            %%% Update Global Data
            time_vector_total = [time_vector_total time_vector_total(end)+time_vector+time_step];
            I_app_total =[I_app_total I_app];
            V_cell_total = [V_cell_total V_cell];
            C_e_total = [C_e_total; C_e];
            C_se_total = [C_se_total; C_se];
            C_s_ave_total = [C_s_ave_total;squeeze(mean(C_s_full,2))];
        
            

           
        
            %% Post CV Rest
            sim_type = 0;
            I_app_rest = zeros(1,round(t_rest_discharge./time_step));
            C_e_init = C_e_end;
            C_s_full_init = C_s_full_end;
            C_se_init = C_se_end;
            I_app_init = I_app_rest;
            c_SEI_init = c_SEI_end;
            c_li_init = c_li_end;
            j_tot_init = j_tot_end;
            psi_s_init = psi_s_end;
            psi_e_init = psi_e_end;
            eps_sol_init = eps_sol_end;
        
        
            
            
            
           
           
            [V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,c_SEI,c_li,I_app,psi_s,psi_e,j_tot,eps_sol,L_SEI,R_SEI,j_SEI,eps_sol_store] = P2D_function_general(Model_Parameters,N_tot,N_shell,order,param,time_step,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);

            
            
            %%% Update local variables
            C_e_end = C_e(end,:);
            C_s_full_end = C_s_full(end,:,:);
            C_se_end = C_se(end,:);
            c_SEI_end = c_SEI(end,:);
            c_li_end = c_li(end,:);
            j_tot_end = j_tot(end,:);
            psi_s_end = psi_s(end,:);
            psi_e_end = psi_e(end,:);
            eps_sol_end = eps_sol;
        
            %%% Update Global Data
            time_vector_total = [time_vector_total time_vector_total(end)+time_vector+time_step];
            I_app_total =[I_app_total I_app];
            V_cell_total = [V_cell_total V_cell];
            C_e_total = [C_e_total; C_e];
            C_se_total = [C_se_total; C_se];
            C_s_ave_total = [C_s_ave_total;squeeze(mean(C_s_full,2))];
            %L_SEI_total =[L_SEI_total L_SEI];

            %% Store Data for EIS test
            idx_store = Test_cycle_count + m;
            C_e_EIS(idx_store,:) = C_e_end;
            C_s_full_EIS(idx_store,:,:) = C_s_full_end;
            C_se_EIS(idx_store,:) = C_se_end;
            c_SEI_EIS(idx_store,:) = c_SEI_end;
            c_li_EIS(idx_store,:) = c_li_end;
            j_tot_EIS(idx_store,:) = j_tot_end;
            psi_s_EIS(idx_store,:) = psi_s_end;
            psi_e_EIS(idx_store,:) = psi_e_end;
            eps_sol_EIS(idx_store,:) = eps_sol_end;

        end
        Test_cycle_count = Test_cycle_count + m;

    end
    toc;
    output_Reigon_Errors = [num2str(k),' Cycles Complete'];
    disp(output_Reigon_Errors);
    q_lost(k) = (1./As_neg).*(F).*sum(c_SEI_end + c_li_end);
    cycle_end_points(k) = length(time_vector_total);
    c_SEI_store(k,:) = c_SEI_end;
    c_li_store(k,:) = c_li_end;
end

%% Plot Results

close all;

figure(1);
plot(time_vector_total,I_app_total,'r');
grid on;
xlabel('Time (s)');
ylabel('Cell Current (A/m^2)');

figure(2)
plot(time_vector_total,V_cell_total,'b');
grid on;
xlabel('Time (s)');
ylabel('Cell Voltage (V)');





%% Conduct EIS Test

if Test_cycle_count > 0

    %% EIS Parameters
    
    Comparison_points = 1000;
    time_step = 1.0;
    order = 5;
    V_cut_off_high = 4.2;
    Cut_off_capacity = 0.1;
    N_frequencies = 40;
    freq_low = 1e-3;
    freq_high = 1e3;
    Points_Per_cycle = 56;
    Cycles_per_freq = 10;
    last_x_cycles = 5;
    Output_cycles = 1;
    time_step = 1;
    Delta_I = 0.5;
    C_rate = 1;
    sim_type = 4;
    
    %% Store EIS Data
    Impedance_plots = zeros(soc_points,N_frequencies);
    idx_end = zeros(soc_points,1);
    Real_Impedance = zeros(soc_points,N_frequencies);
    Imag_Impedance = zeros(soc_points,N_frequencies);
    Freq_Impedance = zeros(soc_points,N_frequencies);
    
    
    %% Pass Parameters
    for idx_store = 1:soc_points
        %idx_store = 1;
        C_e_init = C_e_EIS(idx_store,:);
        C_s_full_init = C_s_full_EIS(idx_store,:,:);
        C_se_init = C_se_EIS(idx_store,:);
        c_SEI_init = c_SEI_EIS(idx_store,:);
        c_li_init = c_li_EIS(idx_store,:);
        j_tot_init = j_tot_EIS(idx_store,:);
        psi_s_init = psi_s_EIS(idx_store,:);
        psi_e_init = psi_e_EIS(idx_store,:);
        eps_sol_init = eps_sol_EIS(idx_store,:);
        
        
        
        %% Conduct EIS Test
        [Impedance_calculated,Frequency_Vector] = EIS_Calculator_Cycling_3(Model_Parameters,N_tot,N_shell,order,param,time_step,freq_low,freq_high,Delta_I, N_frequencies,Points_Per_cycle,Cycles_per_freq,last_x_cycles,C_rate,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);
        
        idx_end(idx_store) = length(Impedance_calculated);
        Impedance_plots(idx_store,1:idx_end(idx_store)) = Impedance_calculated';
         %% Store Impedance Data
        Real_Impedance(idx_store,:) = real(Impedance_plots(idx_store,:));
        Imag_Impedance(idx_store,:) = imag(Impedance_plots(idx_store,:));
        Freq_Impedance(idx_store,:) = Frequency_Vector;
        
        %%% Output the progress of EIS sims
        output_Reigon_Errors = ['EIS Completed Cycle = ',num2str(idx_store)];
        disp(output_Reigon_Errors);
    end
    
    

end

%% Export Drive Cycle Data
if SEI == 1
    
    exampleObject = matfile('Cycler_Output_KL_ageing.mat','Writable',true);
    save Cycler_Output_KL_ageing.mat time_vector_total I_app_total V_cell_total


end

if SEI == 2
    
    exampleObject = matfile('Cycler_Output_DL_ageing.mat','Writable',true);
    save Cycler_Output_DL_ageing.mat time_vector_total I_app_total V_cell_total


end



