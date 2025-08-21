function [Impedance_calculated,Frequency_Vector] = EIS_Calculator_Cycling(Model_Parameters,N_tot,N_shell,order,param,time_step,freq_low,freq_high,Delta_I, N_frequencies,Points_Per_cycle,Cycles_per_freq,last_x_cycles,C_rate,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,psi_s_init,psi_e_init,j_tot_init,eps_sol_init)

    %% EIS Parameters
    
    Frequency_Vector = logspace(log10(freq_low), log10(freq_high), N_frequencies);
    Impedance_calculated = zeros(N_frequencies,1);
    break_con = 0;
    idx_init = 1;
    I_app_init = zeros(1,N_tot);
    Freq_turn = 400;

    for i = idx_init:length(Frequency_Vector)
    
        
    
        %% Make the EIS Current
        [time_vector, I_app_EIS,idx_last_cycle,idx_penultimate] = Single_EIS_Current_Generator(Delta_I,Frequency_Vector(i),time_step,Cycles_per_freq,last_x_cycles,Points_Per_cycle);
        
        %% Generate EIS Response
    
        I_app_init = I_app_EIS;
        time_step_vector = time_vector;
        
        if Frequency_Vector(i) < Freq_turn
        
        [V_cell_EIS,Flag_convergence,time_vector_EIS,C_e,C_s_full,C_se,c_SEI,c_li,I_app_EIS_singles,psi_s,psi_e,j_tot,eps_sol,L_SEI,R_SEI,j_SEI,eps_sol_store] = P2D_function_general(Model_Parameters,N_tot,N_shell,order,param,time_step_vector,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);
        else
        
        [V_cell_EIS,Flag_convergence,time_vector_EIS,C_e,C_s_full,C_se,c_SEI,c_li_pl,I_app_EIS_singles,psi_s,psi_e,j_tot] = P2D_MDN_function_ODE_int(Model_Parameters,N_tot,N_shell,order,param,time_step_vector,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init,eps_sol_init);
        end
        
    
        if Flag_convergence ~= 0
    
            output = ['Break at Freq = ', num2str(Frequency_Vector(i))];
            disp(output)
            break_con = 1;
            
    
            break;
    
    
    
        end
    
        %% Calculate EIS From MDN Solver 
        
        time_vector_SS = time_vector_EIS(idx_last_cycle:idx_penultimate);
        I_app_EIS_SS = I_app_EIS_singles(idx_last_cycle:idx_penultimate);
        V_cell_EIS_SS = V_cell_EIS(idx_last_cycle:idx_penultimate);
        %}
        % FFT
        current_fft = fft(I_app_EIS_SS);
        voltage_fft = fft(V_cell_EIS_SS);
    
        % Get index of first harmonic
        [~, idx] = max(abs(current_fft));
    
        Impedance_calculated(i,:) = -voltage_fft(idx) / current_fft(idx);
    
    
    
         
        
    end


%% Post Processing
    if break_con == 1
        idx_plot = i-1;
    
    else
        idx_plot = i;
    
    end

    Impedance_calculated = Impedance_calculated(idx_init:idx_plot);

end