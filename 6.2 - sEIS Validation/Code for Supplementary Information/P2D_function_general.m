function [V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,I_app] = P2D_function_general(Model_Parameters,N_tot,N_shell,order,param,time_step,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,C_e_init,C_s_full_init,C_se_init,I_app_init,C_rate)

    %% Get Model Parameters

    %Constants
    R_const = 8.314; % Ideal gas constant [J/[K*mol]]
    F = 96487; % Faraday Constant [Coulombs / mol]
    T_amb = 298.15; %Temperature in [K]
    alpha = 0.5;


    A = Model_Parameters(1); % [m] %Cell Area
    brugg = Model_Parameters(2);
    brugg_neg = brugg; %Bruggeman coefficient for tortuosity in negative electrode
    brugg_sep = brugg; %Bruggeman coefficient for tortuosity in separator
    brugg_pos = brugg; %Bruggeman coefficient for tortuosity in positive electrode


    %% Separator/Electrolyte Parameters
    
    L_sep = Model_Parameters(3); % Length of separator [m] 
    C_e_initial = Model_Parameters(4); %Initial Electrolyte Lithium Concentration [mol/[m^3]] 
    D_elec = Model_Parameters(5); % Diffusion Coefficient for Lithium in electrolyte [[m^2]/s] %Defined in model params
    D_sep_e_eff = D_elec;
    t_li = Model_Parameters(6); % Lithium Ion Transference Number [unitless]
    eps_sep_elec = Model_Parameters(7); %Electrolyte Volume Fraction
    kappa_const = Model_Parameters(26);
    dlnfdlnce = Model_Parameters(29);

    %% Anode Parameters
    L_neg = Model_Parameters(8);
    soc_initial_neg = Model_Parameters(9);
    eps_neg_elec = Model_Parameters(10); % Electrolyte Volume Fraction (porosity)
    eps_neg = Model_Parameters(11); % Solid Volume Fraction
    C_max_neg = Model_Parameters(12) ; % Maximum Solid Lithium Concentration [mol/[m^3]]
    C_neg_initial = soc_initial_neg*C_max_neg ;% Initial Lithium Concentration [mol/[m^3]]
    Rs_neg = Model_Parameters(13); % Electrode Particle Size [m] %changed
    D_neg_s = Model_Parameters(14); % Solid Diffusion Coefficient [m^2/s]
    D_neg_eff_s = D_neg_s;
    D_neg_e_eff= D_elec;
    sigma_neg = Model_Parameters(15) ; % Electrical Conductivity [S/m]
    K_0_neg = Model_Parameters(16); % Negative Rate Constant %changed
    As_neg = 3*(eps_neg/Rs_neg) ; % Surface area to volume ration [m-1]
    Cdl_neg = Model_Parameters(27); % Negative Rate Constant %changed
    R_film_neg = Model_Parameters(30);

    %% Cathode Parameters
    L_pos = Model_Parameters(17);
    soc_initial_pos = Model_Parameters(18);
    eps_pos_elec = Model_Parameters(19); %% Electrolyte Volume Fraction 
    eps_pos = Model_Parameters(20); % Solid Volume Fraction
    C_max_pos = Model_Parameters(21) ; % Maximum Solid Lithium Concentration [mol/[m^3]]
    C_pos_initial = soc_initial_pos*C_max_pos ;% Initial Lithium Concentration [mol/[m^3]]
    Rs_pos = Model_Parameters(22); % Electrode Particle Size [m] %changed
    D_pos_s = Model_Parameters(23); % Solid Diffusion Coefficient [m^2/s]
    D_pos_eff_s = D_pos_s;
    D_pos_e_eff = D_elec;
    sigma_pos = Model_Parameters(24) ; % Electrical Conductivity [S/m] From COMSOL
    K_0_pos = Model_Parameters(25); % Negative Rate Constant %changed
    As_pos = 3*(eps_pos/Rs_pos) ; % Surface area to volume ratio [m-1]
    Cdl_pos = Model_Parameters(28); % Negative Rate Constant %changed
    R_film_pos = Model_Parameters(31);
    
    L_tot = L_neg+L_sep+L_pos;

    %% Determin I_1C

    if soc_initial_pos > soc_initial_neg
        Q_pos_initial =  (1/3600)*A*F*L_pos*eps_pos*C_pos_initial*abs(soc_initial_pos - 0); % Units (Ah)
        Q_neg_initial = (1/3600)*A*F*L_neg*eps_neg*C_neg_initial*abs(1.0 - soc_initial_neg ); % Units (Ah/m^2)
    end

    if soc_initial_pos < soc_initial_neg
        Q_pos_initial =  (1/3600)*A*F*L_pos*eps_pos*C_pos_initial*abs(1.0 - soc_initial_pos); % Units (Ah/m^2)
        Q_neg_initial = (1/3600)*A*F*L_neg*eps_neg*C_neg_initial*abs(soc_initial_neg - 0); % Units (Ah/m^2)
    end

    I_1c = max(Q_neg_initial,Q_pos_initial)/A; % A/m^2 1c discharge rate;

    %% Create I_app
    
    % Discharge only
    if sim_type == 0
        I_app = I_app_init;
        time_max = length(I_app_init);
    end

    if sim_type == 1
        time_stop = (1/C_rate)*3600;
        time_max = round(time_stop/time_step);
        time_max = round(2*time_max);
    
        I_const =  C_rate*I_1c; 

        I_app = I_const*ones(1,time_max);

    end

    % Charge only
    if sim_type == 2

        time_stop = (1/C_rate)*3600;
        time_max = round(time_stop/time_step);
        time_max = round(2*time_max);
    
        I_const =  -C_rate*I_1c; 

        I_app = I_const*ones(1,time_max);


    end

     % Charge only
    if sim_type == 2.5

        time_stop = (1/C_rate)*3600;
        time_max = round(time_stop/time_step);
        time_max = round(2*time_max);
    
        I_const =  -C_rate*I_1c; 

        %I_app = I_const*ones(1,time_max);

        I_app_CV = 0;
    end


    % Discharge/charge only
    if sim_type == 3
        
        time_stop = (1/C_rate)*7200;
        time_max = round(time_stop/time_step);
        time_max = round(2*time_max);
    
        I_const =  C_rate*I_1c; 

        I_app = I_const*ones(1,time_max);
    
    end

    % EIS Simulation

    if sim_type == 4

        I_app = I_app_init;
        time_max = length(I_app_init);
        time_step_vector = diff(time_step,1);
        time_step_vector = [time_step_vector time_step_vector(end)];

    end

     %% Run Meshing Algorithm

    del_tot = L_tot/N_tot;
    
    N_neg = round(L_neg/del_tot);                                                                                                                                                                                                                              0;
    del_neg = L_neg/N_neg ;

    N_sep = round(L_sep/del_tot);
    del_sep = L_sep/N_sep;

    N_pos = round(L_pos/del_tot);
    del_pos = L_pos/N_pos ;

    
    N_tot = N_neg + N_sep +N_pos;
    N_tot_active = N_tot;

    %% Apply parameter corrections
    sigma_neg_eff = sigma_neg*(eps_neg); % Effective Electrical Conductivity in negative electrode 
    sigma_pos_eff = sigma_pos*(eps_pos); % Effective Electrical Conductivity in positive electrode 
    
    D_sep_e_eff = (eps_neg_elec^brugg).*D_elec;
    D_neg_e_eff = (eps_neg_elec^brugg).*D_elec;
    D_pos_e_eff = (eps_neg_elec^brugg).*D_elec;


    kappa_const = (eps_neg_elec^brugg).*kappa_const;


    %% Create data containers
    Cdl = [Cdl_neg.*ones(1,N_neg) zeros(1,N_sep) Cdl_pos.*ones(1,N_pos)];
    C_max = [C_max_neg*ones(1,N_neg) zeros(1,N_sep) C_max_pos*ones(1,N_pos)];
    K_0 = [K_0_neg*ones(1,N_neg) zeros(1,N_sep) K_0_pos*ones(1,N_pos)];
    dR_pos = Rs_pos/N_shell; % width of each "shell"
    dR_neg = Rs_neg/N_shell; % width of each "shell"
    Sa_pos = 4*pi*((Rs_pos*(1:N_shell)/N_shell).^2); % outer surface area of each shell
    Sa_neg = 4*pi*((Rs_neg*(1:N_shell)/N_shell).^2); % outer surface area of each shell
    dV_pos = (4/3)*pi*( ((Rs_pos*(1:N_shell)/N_shell).^3) - ((Rs_pos*(0:N_shell-1)/N_shell).^3) ); % vol. of ea. shell
    dV_neg = (4/3)*pi*( ((Rs_neg*(1:N_shell)/N_shell).^3) - ((Rs_neg*(0:N_shell-1)/N_shell).^3) ); % vol. of ea. shell
    diff_psi_s_m = zeros(time_max,N_tot_active+1);
    diff_psi_e_m = zeros(time_max,N_tot_active+1);
    del = [del_neg*ones(1,N_neg) del_sep*ones(1,N_sep) del_pos*ones(1,N_pos)];
    brugg = [brugg_neg*ones(1,N_neg) brugg_sep*ones(1,N_sep) brugg_pos*ones(1,N_pos)];
    
    K_0_init = [K_0_neg*ones(1,N_neg) zeros(1,N_sep) K_0_pos*ones(1,N_pos)];
    kappa = kappa_const*ones(time_max,N_tot_active);
    kd = zeros(time_max,N_tot_active);
    D_e_eff = [D_neg_e_eff*ones(1,N_neg) D_sep_e_eff*ones(1,N_sep) D_pos_e_eff*ones(1,N_pos)];
    D_e = D_e_eff;
    D_s_eff = [D_neg_eff_s*ones(1,N_neg) zeros(1,N_sep) D_pos_eff_s*ones(1,N_pos)];
    D_s = D_s_eff;
    eps_elec = [eps_neg_elec*ones(1,N_neg) eps_sep_elec*ones(1,N_sep) eps_pos_elec*ones(1,N_pos)];
    eps_sol = [eps_neg*ones(1,N_neg) zeros(1,N_sep) eps_pos*ones(1,N_pos)];
    As = [As_neg*ones(1,N_neg) zeros(1,N_sep) As_pos*ones(1,N_pos)];
    Rs = [Rs_neg*ones(1,N_neg) zeros(1,N_sep) Rs_pos*ones(1,N_pos)];
    T = T_amb*ones(time_max,N_tot);
    sigma = [sigma_neg_eff*ones(1,N_neg) zeros(1,N_sep) sigma_pos_eff*ones(1,N_pos)];
    R_film = [R_film_neg*ones(1,N_neg) zeros(1,N_sep) R_film_pos*ones(1,N_pos)];

    %% Dependent Variables
    V_cell = zeros(1,time_max);
    j_flux = zeros(time_max,N_tot_active);
    j_tot = zeros(time_max,N_tot_active);
    j_dl  = zeros(time_max,N_tot_active);
    psi_e = zeros(time_max,N_tot_active);
    psi_s = zeros(time_max,N_tot_active);
    psi_delta = zeros(time_max,N_tot_active);
    eta = zeros(time_max,N_tot_active);
    i_e = zeros(time_max,N_tot_active+1);
    i_s = zeros(time_max,N_tot_active+1);
    ex_current = zeros(time_max,N_tot_active);
    OCP = zeros(time_max,N_tot_active);
    stoic = zeros(time_max,N_tot_active);
    if sum(C_e_init) ==0
        C_e = C_e_initial*ones(time_max+1,N_tot_active);
    else
        C_e = zeros(time_max+1,N_tot_active);
        C_e(1,:) = C_e_init;
    end
    if sum(C_se_init) ==0
        C_se = [C_neg_initial*ones(time_max+1,N_neg) zeros(time_max+1,N_sep) C_pos_initial*ones(time_max+1,N_pos)]; 
    else
        C_se = zeros(time_max+1,N_tot_active);
        C_se(1,:) = C_se_init;
    end
    if sum(sum(squeeze(C_s_full_init))) ==0
        C_s_full = zeros(time_max,N_shell,N_tot_active);
        C_s_full(1,:,1:N_neg) = C_neg_initial*ones(1,N_shell,N_neg);
        C_s_full(1,:,N_neg+N_sep+1:N_tot_active) = C_pos_initial*ones(1,N_shell,N_pos);
    else    
        C_s_full = zeros(time_max,N_shell,N_tot_active);
        C_s_full(1,:,:) = C_s_full_init;
    end
    C_s_ave = zeros(time_max+1,N_tot_active);
    C_s_ave(1,:) = squeeze(mean(C_s_full(1,:,:)))' ;
    diff_psi_s_m = zeros(time_max,N_tot_active+1);
    diff_psi_e_m = zeros(time_max,N_tot_active+1);
    diff_t_psi_s_m = zeros(time_max,N_tot_active);
    diff_t_psi_e_m = zeros(time_max,N_tot_active);
    diff_psi_delta_m = zeros(time_max,N_tot_active);
    time_vector = zeros(1,time_max);

    %% Solver Settings
    tol = 1e-6; % within 0.001 accuracy
    tol_change = 10;

    %tol_I_app = 1e-2;
    tol_I_app = 1e-4;

    if sim_type == 4
         %tol_I_app = 1e-2;
         %tol_I_app = 1e-9;
         %tol_I_app = (1e-3)*max(I_app_init);
         tol_I_app = (tol_I_app)*max(I_app_init);
    end

    %% Create a Flag for breaking iterations
    Flag_convergence = 0;
    

    %% Start Main Loop
    time = 0;
    polarity_change = 0;
    %damp_f = 0.001;
    damp_f = 1;
    case_flux = 3;
    Convergence_max = 5000;


    for i = 1:time_max
        
        %% Get Time_Step if EIS Simulated
        if sim_type == 4
            time_step = time_step_vector(i);
        end
        
        %% Get C_e for current timestep
        C_e_reduced = C_e(i,:);

        %% Calculate kappa vector
    
        [kappa_const] = kappa_calculator(C_e_reduced,param,kappa(i,:),eps_elec,brugg);

        kappa(i,:) = kappa_const;


        %% Get Open Circuit potential and temperature for current time step

        stoic(i,:) = [C_se(i,1:N_neg)./C_max(1:N_neg) zeros(1,N_sep) C_se(i,N_neg+N_sep+1:N_tot_active)./C_max(N_neg+N_sep+1:N_tot_active) ];

        x = stoic(i,1:N_neg);
        y =  stoic(i,N_neg+N_sep+1:N_tot_active);


        [OCP_ref_pos,OCP_ref_neg] = OCP_calculator(x,y,param);

        OCP(i,:) = [OCP_ref_neg zeros(1,N_sep) OCP_ref_pos];

        %Calculate exchange current
        ex_current(i,: ) = [ (K_0(1:N_neg).*sqrt(C_e_reduced(1:N_neg).*(C_max_neg*ones(1,N_neg) - C_se(i,1:N_neg)).*C_se(i,1:N_neg)))  zeros(1,N_sep) ...
              (K_0(N_neg+N_sep+1:N_tot_active).*sqrt(C_e_reduced(N_neg+N_sep+1:N_tot_active).*(C_max_pos*ones(1,N_pos) - C_se(i,N_neg+N_sep+1:N_tot_active)).* C_se(i,N_neg+N_sep+1:N_tot_active)))] ;
        
        
        % Calculate jdl for this time step
        %j_dl(i,:) = damp_f*(Cdl./F).*(diff_t_psi_s_m(i,:) - diff_t_psi_e_m(i,:)); %Changed
        %% Iterative Solver

       if abs(I_app(i)) >= tol_I_app
       

            %% Iterative Solver For Negative Electrode
            j_ave_neg = I_app(i)/(As_neg*F*L_neg)  ;
    
    
            %Set inital guesses
            beta1 = 0.1;
            beta2 = 1.5;
            beta3 = (beta1 + beta2)/2 ;
    
    
            err = ones(3,1);
            Convergence_count = 0;
    
            while norm(err(3),2) > tol
    
                % Create Guesses
                j1 = beta1*j_ave_neg;
                j2 = beta2*j_ave_neg;
            
                j3 = beta3*j_ave_neg;
    
                j_attempts = [j1 j2 j3];
    
                j_flux_temp = zeros(3,N_tot_active);
                j_tot_temp = zeros(3,N_tot_active);
                j_dl_temp = zeros(3,N_tot_active);
                psi_e_temp = zeros(3,N_tot_active);
                psi_s_temp = zeros(3,N_tot_active);
                psi_delta_temp = zeros(3,N_tot_active);
                eta_temp = zeros(3,N_tot_active);
                diff_psi_s_m_temp = zeros(3,N_tot_active+1);
                diff_psi_e_m_temp = zeros(3,N_tot_active+1);
    
                if i>1 
                    j_SEI = j_flux(i-1,:);

                else
                    j_SEI = j_flux(i,:);
                end
                j_flux_temp(:,1) = j_attempts';
    
                %Set Reference Potential
                psi_e_temp(:,1) = 0;
                
                %Caclulate eta for j_flux guess 
                eta_temp(:,1) = (R_const*T(i,1)/(alpha*F)).*asinh((1/2).*(j_flux_temp(:,1)./ex_current(i,1))) ;
    
                %Calculate psi_s
                psi_s_temp(:,1) = eta_temp(:,1) + psi_e_temp(:,1) + OCP(i,1) + F*R_film(1)*j_SEI(1);
    
                %Calculate psi_delta
                psi_delta_temp(:,1) = psi_s_temp(:,1) - psi_e_temp(:,1);
                
                % Calculate double layer flux
                if i == 1
                    j_dl_temp(:,1) = (Cdl(1)/F)*((psi_delta_temp(:,1))./time_step);
                else
                    j_dl_temp(:,1) = (Cdl(1)/F)*((psi_delta_temp(:,1) - psi_delta(i-1,1))./time_step);
                end
                
                % Calculate total flux
                j_tot_temp(:,1) = j_flux_temp(:,1) + j_dl_temp(:,1);
                
                % Set Boundary Conditions
                diff_psi_s_m_temp(:,1) = -I_app(i)/(sigma_neg_eff) ;
            
                diff_psi_e_m_temp(:,1) = 0;
    
                A = 0.5*eye(N_tot_active+1,N_tot_active);
            
                B = A + circshift(A,[1 0]);
                
                del_mid = B*del';
                
                del_mid = del_mid';
                
                diffx_C_e_log = diff(log(C_e_reduced),1,2) ;
        
                diffx_C_e_log = [0 diffx_C_e_log 0]./del_mid;
                
    
                for j = 1:N_neg-1
    
                    
                    
                    diff_psi_s_m_temp(:,j+1) = diff_psi_s_m_temp(:,j) + (F*As(j)*del(j)/sigma(j)).*(j_tot_temp(:,j)); %% Changed


                    psi_s_temp(:,j+1) = psi_s_temp(:,j) + del(j)*diff_psi_s_m_temp(:,j+1);
    
    
                    %Calculate psi_e for j+1 element
                    kappa_k_plus_half = (kappa(i,j+1) + kappa(i,j))/2 ;
                    if j == 1
                        kappa_k_minus_half = kappa_k_plus_half;
                    else
                        kappa_k_minus_half = (kappa(i,j) + kappa(i,j-1))/2 ;
                    end
    
                    T_temp = T_amb;
    
                    kd_k_plus_half = ((2*kappa_k_plus_half*R_const*T_temp)/(F))*(1-t_li) ;
    
                    if j == 1
                        kd_k_minus_half = kd_k_plus_half;
                    else
                        kd_k_minus_half = ((2*kappa_k_minus_half*R_const*T_temp)/(F))*(1-t_li) ;
                    end
    
        
                    diff_psi_e_m_temp(:,j+1) = (1/kappa_k_plus_half).*( (kappa_k_minus_half).*diff_psi_e_m_temp(:,j) + (1+dlnfdlnce)*( kd_k_plus_half.*diffx_C_e_log(j+1) - kd_k_minus_half.*diffx_C_e_log(j) ) - (F.*As(j).*del(j).*j_tot_temp(:,j)) );

                            

                    psi_e_temp(:,j+1) = psi_e_temp(:,j) + del(j)*diff_psi_e_m_temp(:,j+1);
    
                    % Determine j_fux for j+1 element
                    eta_temp(:,j+1) = psi_s_temp(:,j+1) - psi_e_temp(:,j+1) - OCP(i,j+1) - F*R_film(j+1)*j_SEI(j+1);
    
                    j_flux_temp(:,j+1) = 2.*ex_current(i,j+1).*sinh((alpha*F/(R_const*T(i,j+1)))*eta_temp(:,j+1)) ;

                    %Calculate psi_delta
                    psi_delta_temp(:,j+1) = psi_s_temp(:,j+1) - psi_e_temp(:,j+1);

                    % Calculate double layer flux
                    if i == 1
                        j_dl_temp(:,j+1) = (Cdl(j+1)/F)*((psi_delta_temp(:,j+1))./time_step);
                    else
                        j_dl_temp(:,j+1) = (Cdl(j+1)/F)*((psi_delta_temp(:,j+1) - psi_delta(i-1,j+1))./time_step);
                    end

                    % Calculate total flux
                    j_tot_temp(:,j+1) = j_flux_temp(:,j+1) + j_dl_temp(:,j+1);


                end
    
                err =  1  - mean(j_tot_temp(:,1:N_neg),2)./j_ave_neg;
                
    
                if sign(err(1)*err(3)) == sign(err(2)*err(3))
        
                    %error("Doesnt Converge")
                    %beta1 = beta1/beta1;
                    beta1 = beta1 - 3;%Changed from 0.1
                    beta2 = beta2 + 3;
                    beta3 = (beta1+beta2)/2 ; 
        
                Convergence_count = Convergence_count +1;
                else
        
                    %%% Hybeid Bisection + Regula-falsi Root Finding
                    if norm(err(1:3),2) > tol_change
                        if err(1)*err(3) < 0 
                
                            beta1 = beta1;
                
                            beta2 = beta3;
                        
                            beta3 = (beta1 +beta2)/2;
                
                        
                        end
    
                        if err(2)*err(3) < 0
                            beta1 = beta3;
                            beta2 = beta2;
                            beta3 = (beta1 +beta2)/2;
                
                        end
    
                    else
    
                        if err(1)*err(3) < 0 
                    
                            beta_temp = (beta1*err(2) - beta2*err(1))/(err(2) - err(1)) ;
                            beta1 = beta1;
                    
                            beta2 = beta3;
                    
                            beta3 = beta_temp;
                    
                        end
            
                        if err(2)*err(3) < 0
                            beta_temp = (beta1*err(2) - beta2*err(1))/(err(2) - err(1)) ;
                            beta1 = beta3;
                            beta2 = beta2;
                            beta3 = beta_temp;
                        end
    
                    end
                Convergence_count = Convergence_count +1;
        
                end
    
                if Convergence_count>Convergence_max
                    Flag_convergence = 1;
                    break;
                end
    
            end
            % Update Fluxes
            j_tot(i,:) = j_tot_temp(3,:);
            j_flux(i,:) = j_flux_temp(3,:);
            j_dl(i,:) = j_dl_temp(3,:);
    
            
            if Convergence_count>Convergence_max
                Flag_convergence = 1;
                break;
            end
    
            if sum(isinf(j_flux(i,1:N_neg))) > 0
                Flag_convergence = 1;
                break;
    
            end
    
            if sum(isnan(j_flux(i,1:N_neg))) > 0
                Flag_convergence = 1;
                break;
    
            end
    
            
            if isreal(j_flux(i,1:N_neg)) == 0
                Flag_convergence = 1;
                break;
    
            end

           %% Caclulate psi_e from current collector to first element of cathode
           diff_psi_e_m(i,1) = 0;
           psi_e(i,1) = 0;
           for j = 1:N_neg+N_sep
    
               %Calculate psi_e for j+1 element
                kappa_k_plus_half = (kappa(i,j+1) + kappa(i,j))/2 ;
                
                if j == 1
                   
                    kappa_k_minus_half = kappa_k_plus_half;
                else
                    kappa_k_minus_half = (kappa(i,j) + kappa(i,j-1))/2 ;
                end
    
                T_temp = T_amb;
    
                kd_k_plus_half = ((2*kappa_k_plus_half*R_const*T_temp)/(F))*(1-t_li) ;
        
                if j == 1
                      kd_k_minus_half = kd_k_plus_half;
                else
                      T_temp = T_amb;
                      kd_k_minus_half = ((2*kappa_k_minus_half*R_const*T_temp)/(F))*(1-t_li) ;
                end
    
    
                diff_psi_e_m(i,j+1) = (1/kappa_k_plus_half)*((kappa_k_minus_half)*diff_psi_e_m(i,j) + (kd_k_plus_half*diffx_C_e_log(j+1) - kd_k_minus_half*diffx_C_e_log(j)) - (F*As(j)*del(j)*j_tot(i,j)));
        
                

                psi_e(i,j+1) = psi_e(i,j) + del(j)*diff_psi_e_m(i,j+1);
    
    
           end

           %% Iterative Solver for Positive Electrode
           j_ave_pos = -I_app(i)/(As_pos*F*L_pos) ;
    
           beta1 = 0.1;
           beta2 = 1.9;
           beta3 = (beta1 +beta2)/2;
            
           err = ones(3,1);
    
          
           Convergence_count = 0;
    
           while norm(err(3),2) > tol
                
               % Create Guesses
                j1 = beta1*j_ave_pos;
                j2 = beta2*j_ave_pos;
            
                j3 = beta3*j_ave_pos;
    
                j_attempts = [j1 j2 j3];
    
                j_flux_temp = ones(3,N_tot_active)*diag(j_flux(i,:));
                j_tot_temp = ones(3,N_tot_active)*diag(j_tot(i,:));
                j_dl_temp = ones(3,N_tot_active)*diag(j_dl(i,:));
                psi_e_temp = ones(3,N_tot_active)*diag(psi_e(i,:));
                %psi_s_temp = zeros(3,N_tot_active);
                eta_temp = zeros(3,N_tot_active);
                %diff_psi_s_m_temp = zeros(3,N_tot_active+1);
                diff_psi_e_m_temp = ones(3,N_tot_active+1)*diag(diff_psi_e_m(i,:));
    
                j_flux_temp(:,N_neg+N_sep+1) = j_attempts';
    
                %Caclulate eta for j_flux guess 
                eta_temp(:,N_neg+N_sep+1) = (R_const*T(i,N_neg+N_sep+1)/(alpha*F)).*asinh((1/2).*(j_flux_temp(:,N_neg+N_sep+1)./ex_current(i,N_neg+N_sep+1))) ;
    
                %Calculate psi_s
                psi_s_temp(:,N_neg+N_sep+1) = eta_temp(:,N_neg+N_sep+1) + psi_e_temp(:,N_neg+N_sep+1) + OCP(i,N_neg+N_sep+1) + F*R_film(N_neg+N_sep+1)*j_SEI(N_neg+N_sep+1);
                
                %Calculate psi_delta
                psi_delta_temp(:,N_neg+N_sep+1) = psi_s_temp(:,N_neg+N_sep+1) - psi_e_temp(:,N_neg+N_sep+1);
                
                % Calculate double layer flux
                if i == 1
                    j_dl_temp(:,N_neg+N_sep+1) = (Cdl(N_neg+N_sep+1)/F)*((psi_delta_temp(:,N_neg+N_sep+1))./time_step);
                else
                    j_dl_temp(:,N_neg+N_sep+1) = (Cdl(N_neg+N_sep+1)/F)*((psi_delta_temp(:,N_neg+N_sep+1) - psi_delta(i-1,N_neg+N_sep+1))./time_step);
                end
                
                % Calculate total flux
                j_tot_temp(:,N_neg+N_sep+1) = j_flux_temp(:,N_neg+N_sep+1) + j_dl_temp(:,N_neg+N_sep+1);

                % Set Boundary Conditions
                diff_psi_s_m(:,N_neg+N_sep+1) = 0 ;
    
                for j = N_neg+N_sep+1:N_tot_active-1
               
                    diff_psi_s_m_temp(:,j+1) = diff_psi_s_m_temp(:,j) + (F*As(j)*del(j)/sigma(j)).*(j_tot_temp(:,j)); %% Changed

                    psi_s_temp(:,j+1) = psi_s_temp(:,j) + del(j)*diff_psi_s_m_temp(:,j+1);
    
                    kappa_k_plus_half = (kappa(i,j+1) + kappa(i,j))/2 ;
    
                    kappa_k_minus_half = (kappa(i,j) + kappa(i,j-1))/2 ;
    
    
                    T_temp = T_amb;
    
                    kd_k_plus_half = ((2*kappa_k_plus_half*R_const*T_temp)/(F))*(1-t_li) ;
        
                    kd_k_minus_half = ((2*kappa_k_minus_half*R_const*T_temp)/(F))*(1-t_li) ;
        

                    diff_psi_e_m_temp(:,j+1) = (1/kappa_k_plus_half).*( (kappa_k_minus_half).*diff_psi_e_m_temp(:,j) + (1+dlnfdlnce)*( kd_k_plus_half.*diffx_C_e_log(j+1) - kd_k_minus_half.*diffx_C_e_log(j) ) - (F.*As(j).*del(j).*j_tot_temp(:,j)) );
                                 
                    psi_e_temp(:,j+1) = psi_e_temp(:,j) + del(j)*diff_psi_e_m_temp(:,j+1);
    
                    % Determine j_fux for j+1 element
                    eta_temp(:,j+1) = psi_s_temp(:,j+1) - psi_e_temp(:,j+1) - OCP(i,j+1) - F*R_film(j+1)*j_SEI(j+1);
    
                    j_flux_temp(:,j+1) = 2.*ex_current(i,j+1).*sinh((alpha*F/(R_const*T(i,j+1)))*eta_temp(:,j+1)) ;

                    %Calculate psi_delta
                    psi_delta_temp(:,j+1) = psi_s_temp(:,j+1) - psi_e_temp(:,j+1);

                    % Calculate double layer flux
                    if i == 1
                        j_dl_temp(:,j+1) = (Cdl(j+1)/F)*((psi_delta_temp(:,j+1))./time_step);
                    else
                        j_dl_temp(:,j+1) = (Cdl(j+1)/F)*((psi_delta_temp(:,j+1) - psi_delta(i-1,j+1))./time_step);
                    end

                    % Calculate total flux
                    j_tot_temp(:,j+1) = j_flux_temp(:,j+1) + j_dl_temp(:,j+1);

    
                end
    
                err =  1  - mean(j_tot_temp(:,N_neg+N_sep+1:N_tot_active),2)./j_ave_pos;
                
                
                if sign(err(1)*err(3)) == sign(err(2)*err(3))
        
                    %error("Doesnt Converge")
                    %beta1 = beta1/beta1;
                    beta1 = beta1 - 3;%Changed from 0.1
                    beta2 = beta2 + 3;
                    beta3 = (beta1+beta2)/2 ; 
        
                Convergence_count = Convergence_count +1;
                else
        
                    %%% Hybeid Bisection + Regula-falsi Root Finding
                    if norm(err(1:3),2) > tol_change
                        if err(1)*err(3) < 0 
                
                            beta1 = beta1;
                
                            beta2 = beta3;
                        
                            beta3 = (beta1 +beta2)/2;
                
                        
                        end
    
                        if err(2)*err(3) < 0
                            beta1 = beta3;
                            beta2 = beta2;
                            beta3 = (beta1 +beta2)/2;
                
                        end
    
                    else
    
                        if err(1)*err(3) < 0 
                    
                            beta_temp = (beta1*err(2) - beta2*err(1))/(err(2) - err(1)) ;
                            beta1 = beta1;
                    
                            beta2 = beta3;
                    
                            beta3 = beta_temp;
                    
                        end
            
                        if err(2)*err(3) < 0
                            beta_temp = (beta1*err(2) - beta2*err(1))/(err(2) - err(1)) ;
                            beta1 = beta3;
                            beta2 = beta2;
                            beta3 = beta_temp;
                        end
    
                    end
                Convergence_count = Convergence_count +1;
    
                
        
                end
    
                if Convergence_count>Convergence_max
                    Flag_convergence = 2;
                    break;
                end
           end
    
           % Update Fluxes
           j_tot(i,:) = j_tot_temp(3,:);
           j_flux(i,:) = j_flux_temp(3,:);
           j_dl(i,:) = j_dl_temp(3,:);
           diff_psi_e_m(i,:) = diff_psi_e_m_temp(3,:);
           diff_psi_s_m(i,:) = diff_psi_s_m_temp(3,:);

           psi_s(i,:) = psi_s_temp(3,:);
           psi_e(i,:) = psi_e_temp(3,:);
           psi_delta(i,:) = psi_s(i,:) - psi_e(i,:);

    
            
            if Convergence_count>Convergence_max
                Flag_convergence = 2;
                break;
            end
    
            if sum(isinf(j_flux(i,:))) > 0
                Flag_convergence = 2;
                break;
    
            end
    
            if sum(isnan(j_flux(i,:))) > 0
                Flag_convergence = 2;
                break;
    
            end
    
            
            if isreal(j_flux(i,:)) == 0
                Flag_convergence = 2;
                break;
    
            end
       

       
       
       %%% End of Iterative Solver


 

       end  
     
       if abs(I_app(i)) < tol_I_app
           if case_flux == 1
           % Calculate j_tot
           j_tot(i,:) = zeros(1,N_tot_active);

           
           %% Caclulate psi_e 
           A = 0.5*eye(N_tot_active+1,N_tot_active);
            
           B = A + circshift(A,[1 0]);
            
           del_mid = B*del';
            
           del_mid = del_mid';
            
           diffx_C_e_log = diff(log(C_e_reduced),1,2) ;
    
           diffx_C_e_log = [0 diffx_C_e_log 0]./del_mid;

           diff_psi_e_m(i,1) = 0;
           psi_e(i,1) = 0;

           for j = 1:N_tot_active -1
    
               %Calculate psi_e for j+1 element
                kappa_k_plus_half = (kappa(i,j+1) + kappa(i,j))/2 ;
                
                if j == 1
                   
                    kappa_k_minus_half = kappa_k_plus_half;
                else
                    kappa_k_minus_half = (kappa(i,j) + kappa(i,j-1))/2 ;
                end
    
                T_temp = T_amb;
    
                kd_k_plus_half = ((2*kappa_k_plus_half*R_const*T_temp)/(F))*(1-t_li) ;
        
                if j == 1
                      kd_k_minus_half = kd_k_plus_half;
                else
                      T_temp = T_amb;
                      kd_k_minus_half = ((2*kappa_k_minus_half*R_const*T_temp)/(F))*(1-t_li) ;
                end
    
    
                diff_psi_e_m(i,j+1) = (1/kappa_k_plus_half)*((kappa_k_minus_half)*diff_psi_e_m(i,j) + (kd_k_plus_half*diffx_C_e_log(j+1) - kd_k_minus_half*diffx_C_e_log(j)) - (F*As(j)*del(j)*j_tot(i,j)));
        

                psi_e(i,j+1) = psi_e(i,j) + del(j)*diff_psi_e_m(i,j+1);
    
    
           end

           %% Calculate psi_s (Need a better option for this)(psi_delta approach?)
           if i>1
                psi_s(i,:) = psi_s(i-1,:);
           end
           
           eta(i,1:N_neg) = psi_s(i,1:N_neg) - psi_e(i,1:N_neg) - OCP(i,1:N_neg);

           eta(i,N_neg+N_sep+1:N_tot_active) = psi_s(i,N_neg+N_sep+1:N_tot_active) - psi_e(i,N_neg+N_sep+1:N_tot_active) - OCP(i,N_neg+N_sep+1:N_tot_active);

           

          
           if i>1
                diff_t_psi_s_m(i,:) = diff(psi_s(i-1:i,:),1)./time_step;
                diff_t_psi_e_m(i,:) = diff(psi_e(i-1:i,:),1)./time_step;
                j_dl(i,:) = (Cdl./F).*(diff_t_psi_s_m(i,:) - diff_t_psi_e_m(i,:)); 
                j_flux(i,:) = -j_dl(i,:);
                
           else
               diff_t_psi_s_m(i,:) = (psi_s(i,:) - OCP(i,:))./time_step;
               diff_t_psi_e_m(i,:) = (psi_e(i,:)-zeros(1,N_tot_active))./time_step;
               j_dl(i,:) = (Cdl./F).*(diff_t_psi_s_m(i,:) - diff_t_psi_e_m(i,:));
               j_flux(i,:) = -j_dl(i,:);
           end

           end
           
           if case_flux ==2 

               j_flux(i,:) = zeros(1,N_tot_active);
    
               eta(i,:) = (R_const/(alpha*F))*T(i,:).*asinh((1/2).*(j_flux(i,:)./ex_current(i,:))) ;
                
               if i>1
                    psi_e(i,:) = psi_e(i-1,:);
               end
    
               psi_s(i,1:N_neg) = eta(i,1:N_neg) + psi_e(i,1:N_neg) + OCP(i,1:N_neg);
    
               psi_s(i,N_neg+N_sep+1:N_tot_active) = eta(i,N_neg+N_sep+1:N_tot_active) + psi_e(i,N_neg+N_sep+1:N_tot_active) + OCP(i,N_neg+N_sep+1:N_tot_active);
    
               if i>1
                    diff_t_psi_s_m(i,:) = diff(psi_s(i-1:i,:),1)./time_step;
                    diff_t_psi_e_m(i,:) = diff(psi_e(i-1:i,:),1)./time_step;
                    j_dl(i,:) = damp_f*(Cdl./F).*(diff_t_psi_s_m(i,:) - diff_t_psi_e_m(i,:)); %Changed
                    
               else
                   diff_t_psi_s_m(i,:) = (psi_s(i,:) - OCP(i,:))./time_step;
                   diff_t_psi_e_m(i,:) = (psi_e(i,:)-zeros(1,N_tot_active))./time_step;
                   j_dl(i,:) = damp_f*(Cdl./F).*(diff_t_psi_s_m(i,:) - diff_t_psi_e_m(i,:)); %Changed
               end
    
               j_tot(i,:) = j_flux(i,:) + j_dl(i,:);
           
           end

           if case_flux == 3
                % Calculate  psi_delta then calculate fluxes
                if i == 1
                   j_dl_temp = zeros(1,N_tot_active);

                   psi_delta_temp = [OCP(i,1:N_neg) zeros(1,N_sep) OCP(i,N_neg+N_sep+1:N_tot_active)];

                else
                    j_dl_temp = j_dl(i-1,:);

                    psi_delta_temp = psi_delta(i-1,:);

                end
                diff_t_psi_delta = (1./(Cdl*F)).*j_dl_temp;

                psi_delta(i,:) = psi_delta_temp + diff_t_psi_delta*time_step ;

                eta(i,:) = psi_delta(i,:) - OCP(i,:);

                eta(i,N_neg+1:N_neg+N_sep) = zeros(1,N_sep);

                j_flux(i,:) = 2.*ex_current(i,:).*sinh((alpha*F./(R_const.*T(i,:))).*eta(i,:)) ;

                j_dl(i,:) = - j_flux(i,:);

                j_tot(i,:) = j_flux(i,:) + j_dl(i,:);

               %% Caclulate psi_e 
               A = 0.5*eye(N_tot_active+1,N_tot_active);
                
               B = A + circshift(A,[1 0]);
                
               del_mid = B*del';
                
               del_mid = del_mid';
                
               diffx_C_e_log = diff(log(C_e_reduced),1,2) ;
        
               diffx_C_e_log = [0 diffx_C_e_log 0]./del_mid;
    
               diff_psi_e_m(i,1) = 0;
               psi_e(i,1) = 0;
    
               for j = 1:N_tot_active -1
        
                   %Calculate psi_e for j+1 element
                    kappa_k_plus_half = (kappa(i,j+1) + kappa(i,j))/2 ;
                    
                    if j == 1
                       
                        kappa_k_minus_half = kappa_k_plus_half;
                    else
                        kappa_k_minus_half = (kappa(i,j) + kappa(i,j-1))/2 ;
                    end
        
                    T_temp = T_amb;
        
                    kd_k_plus_half = ((2*kappa_k_plus_half*R_const*T_temp)/(F))*(1-t_li) ;
            
                    if j == 1
                          kd_k_minus_half = kd_k_plus_half;
                    else
                          T_temp = T_amb;
                          kd_k_minus_half = ((2*kappa_k_minus_half*R_const*T_temp)/(F))*(1-t_li) ;
                    end
        
        
                    diff_psi_e_m(i,j+1) = (1/kappa_k_plus_half)*((kappa_k_minus_half)*diff_psi_e_m(i,j) + (kd_k_plus_half*diffx_C_e_log(j+1) - kd_k_minus_half*diffx_C_e_log(j)) - (F*As(j)*del(j)*j_tot(i,j)));
            
    
                    psi_e(i,j+1) = psi_e(i,j) + del(j)*diff_psi_e_m(i,j+1);
        
        
               end

               psi_s(i,:) = psi_delta(i,:) + psi_e(i,:);

               psi_s(i,N_neg+1:N_neg+N_sep) = zeros(1,N_sep);

           end

           if case_flux == 4
                if i == 1
                    j_dl(i,:) = zeros(1,N_tot_active);

                else
                    j_dl(i,:) = j_dl(i-1,:);
                    
                end

                j_flux(i,:) = - j_dl(i,:);

                j_tot(i,:) = j_flux(i,:) + j_dl(i,:);


                %% Caclulate eta
                eta(i,:) = (R_const/(alpha*F))*T(i,:).*asinh((1/2).*(j_flux(i,:)./ex_current(i,:))) ;


                %% Caclulate psi_e 
               A = 0.5*eye(N_tot_active+1,N_tot_active);
                
               B = A + circshift(A,[1 0]);
                
               del_mid = B*del';
                
               del_mid = del_mid';
                
               diffx_C_e_log = diff(log(C_e_reduced),1,2) ;
        
               diffx_C_e_log = [0 diffx_C_e_log 0]./del_mid;
    
               diff_psi_e_m(i,1) = 0;
               psi_e(i,1) = 0;
    
               for j = 1:N_tot_active -1
        
                   %Calculate psi_e for j+1 element
                    kappa_k_plus_half = (kappa(i,j+1) + kappa(i,j))/2 ;
                    
                    if j == 1
                       
                        kappa_k_minus_half = kappa_k_plus_half;
                    else
                        kappa_k_minus_half = (kappa(i,j) + kappa(i,j-1))/2 ;
                    end
        
                    T_temp = T_amb;
        
                    kd_k_plus_half = ((2*kappa_k_plus_half*R_const*T_temp)/(F))*(1-t_li) ;
            
                    if j == 1
                          kd_k_minus_half = kd_k_plus_half;
                    else
                          T_temp = T_amb;
                          kd_k_minus_half = ((2*kappa_k_minus_half*R_const*T_temp)/(F))*(1-t_li) ;
                    end
        
        
                    diff_psi_e_m(i,j+1) = (1/kappa_k_plus_half)*((kappa_k_minus_half)*diff_psi_e_m(i,j) + (kd_k_plus_half*diffx_C_e_log(j+1) - kd_k_minus_half*diffx_C_e_log(j)) - (F*As(j)*del(j)*j_tot(i,j)));
            
    
                    psi_e(i,j+1) = psi_e(i,j) + del(j)*diff_psi_e_m(i,j+1);
        
        
               end

               %% Caclulate psi_s

               psi_s(i,:) = eta(i,:) + psi_e(i,:) + OCP(i,:);

               psi_s(i,N_neg+1:N_neg+N_sep) = zeros(1,N_sep);

                

           end
       end

      

       %% Calculate  i_e, i_s
       %%% Calculate i_e and i_s
       A = tril(ones(N_tot_active));
       B = [zeros(1,N_tot);A];
       i_e(i,:) = (zeros(N_tot_active+1,1) + B*((F*As.*del.*j_tot(i,:))'))' ;
       i_s(i,:) = (I_app(i)*ones(N_tot_active+1,1) - B*((F*As.*del.*j_tot(i,:))'))' ; 

       
       %% Calculate Time Derivatives %%%

       

       %% Calculate C_e
       
       %%% Use standard FVM with ode15s %%%

       % Make Extended Del Matrix
       del_ave_neg_sep = (del(N_neg) +  del(N_neg+1))/2 ;
       del_ave_sep_pos = (del(N_neg + N_sep) +  del(N_neg+N_sep+1))/2 ;
    
       del_extended = [del(1:N_neg) del_ave_neg_sep del(N_neg+2:N_neg+N_sep) del_ave_sep_pos del(N_neg+N_sep+1:N_tot_active)];
       
       % Creat Harmonic Mean at anode/seperator interface and seperator/cathode
       % interface
    
       beta_neg_sep = del(N_neg)/(del(N_neg) + del(N_neg+1));
       beta_sep_pos = del(N_neg+N_sep)/(del(N_neg+N_sep) + del(N_neg+N_sep+1));
       
       D_HM_neg_sep = (D_e(N_neg)*D_e(N_neg+1))/( beta_neg_sep*D_e(N_neg+1) + (1 - beta_neg_sep)*D_e(N_neg) );
       
       D_HM_sep_pos = (D_e(N_neg+N_sep)*D_e(N_neg+N_sep+1))/( beta_sep_pos*D_e(N_neg+N_sep+1) + (1 - beta_sep_pos)*D_e(N_neg+N_sep) );
       
       D_e_eff_ext = [D_e(1:N_neg) D_HM_neg_sep D_e(N_neg+2:N_neg+N_sep) D_HM_sep_pos D_e(N_neg+N_sep+1:N_tot_active)];
       
       %Caclulate using FVM
       op = odeset('reltol',tol);
       %[t,C_e_temp] = ode15s(@(t,C_e_temp) ( (1./eps_elec').*( ((1./del).*diff((D_e_eff_ext.*[ 0 ; diff(C_e_temp,1) ; 0 ]')./del_extended,1))' + ((1-t_li)*(As).*j_flux(i,:))' + (((As.*Cdl)./F).*(diff_t_psi_s_m(i,:)-diff_t_psi_e_m(i,:)))' )), [0 time_step], C_e(i,:)', op);
       
       [t,C_e_temp] = ode15s(@(t,C_e_temp) ( (1./eps_elec').*( ((1./del).*diff((D_e_eff_ext.*[ 0 ; diff(C_e_temp,1) ; 0 ]')./del_extended,1))' + ((1-t_li)*(As).*(j_tot(i,:)))' )), [0 time_step], C_e(i,:)', op);  %% Changed

       C_e(i+1,:) = C_e_temp(end,:); % update
    
       
       difft_C_e_m(i,:) = (1/time_step)*(C_e(i+1,:)-C_e(i,:));

       %% Calculate C_s, C_se and C_s_ave
       if order == 4
           %Calculate Solid Concentration Change In Anode
           N_flux_neg(i,:,:) = diff(C_s_full(i,:,1:N_neg),1,2)/dR_neg;
           
    
           N_flux_neg_temp = -squeeze(N_flux_neg(i,:,:))*diag(D_s(1:N_neg));
           M_neg =   diag(Sa_neg(1:end-1),0)*N_flux_neg_temp ;
           
           difft_C_s_neg = diag(1./dV_neg,0)*([zeros(1,N_neg) ; M_neg] - [M_neg ; zeros(1,N_neg)]); 
           
           %Calculate Solid Concentration Change In Cathode
           N_flux_pos(i,:,:) = diff(C_s_full(i,:,N_neg+N_sep+1:N_tot_active),1,2)/dR_pos;
    
           N_flux_pos_temp = -squeeze(N_flux_pos(i,:,:))*diag(D_s(N_neg+N_sep+1:N_tot_active));
           
           M_pos =   diag(Sa_pos(1:end-1),0)*N_flux_pos_temp ;
    
           difft_C_s_pos = diag(1./dV_pos,0)*([zeros(1,N_pos) ; M_pos] - [M_pos ; zeros(1,N_pos)]); 
            
           %Calculate Solid Concentration Change in Entire Cell
           difft_C_s = [difft_C_s_neg zeros(N_shell,N_sep) difft_C_s_pos] ;
    
           
           
           C_s_shells_next = difft_C_s*time_step;
           
           % Add j_flux effect
           C_s_shells_next(end,:) = C_s_shells_next(end,:) - time_step*j_flux(i,:).*[(Sa_neg(end)/dV_neg(end))*ones(1,N_neg) zeros(1,N_sep) (Sa_pos(end)/dV_pos(end))*ones(1,N_pos)];
            
           C_s_full(i+1,:,:) = C_s_full(i,:,:) +  permute(C_s_shells_next,[3 1 2]);
           
           C_se(i+1,:) = squeeze(C_s_full(i+1,N_shell,:))'; %Set Surface Concentration
          
           C_s_ave(i+1,:) = squeeze(mean(C_s_full(i+1,:,:)))' ;

       end

       if order == 5
            %C_s_shells = squeeze(C_s_full(i,:,:));

            %op = odeset('reltol',tol);
                
            %[t,C_s_shells_temp] = ode15s(@(t,C_s_shells_temp) ((D_s(1:N_neg)'*(1./dV_neg))).*[(Sa_neg(1).*diff(C_s_shells_temp(:,1:2),1,2)./dR_neg)' ;...
                %(diff(diag(Sa_neg(1:end-1))'*(diff(C_s_shells_temp,1,2)./dR_neg)')) ; (-Sa_neg(end).*(j_flux(i,1:N_neg)./D_s(1:N_neg))' - Sa_neg(end-1).*diff(C_s_shells_temp(:,end-1:end),1,2)./dR_neg)' ]',[0 time_step], C_s_shells(:,1:N_neg)',op);
             
            %[t,C_s_shells_temp] = ode15s(@(t,C_s_shells_temp) ((D_s(1:N_neg)'*(1./dV_neg))).*[(Sa_neg(1).*diff(C_s_shells_temp(:,1:2),1,2)./dR_neg)' ;...
                %(diff(diag(Sa_neg(1:end-1))'*(diff(C_s_shells_temp,1,2)./dR_neg)')) ; (-Sa_neg(end).*(j_flux(i,1:N_neg)./D_s(1:N_neg))' - Sa_neg(end-1).*diff(C_s_shells_temp(:,end-1:end),1,2)./dR_neg)' ]',[0 time_step], C_s_full(i,:,:)',op);
                
            C_s_shells = squeeze(C_s_full(i,:,:));

            C_s_shells_temp  = C_s_shells(:,1:N_neg)';

            diff_t_C_s_neg = ((D_s(1:N_neg)'*(1./dV_neg))).*[(Sa_neg(1).*diff(C_s_shells_temp(:,1:2),1,2)./dR_neg)' ;...
                        (diff(diag(Sa_neg(1:end-1))'*(diff(C_s_shells_temp,1,2)./dR_neg)')) ; (-Sa_neg(end).*(j_flux(i,1:N_neg)./D_s(1:N_neg))' - Sa_neg(end-1).*diff(C_s_shells_temp(:,end-1:end),1,2)./dR_neg)' ]';

            

            C_s_shells_temp  = C_s_shells(:,N_neg+N_sep+1:N_tot)';

            diff_t_C_s_pos = ((D_s(N_neg+N_sep+1:N_tot)'*(1./dV_pos))).*[(Sa_pos(1).*diff(C_s_shells_temp(:,1:2),1,2)./dR_pos)' ;...
                        (diff(diag(Sa_pos(1:end-1))'*(diff(C_s_shells_temp,1,2)./dR_pos)')) ; (-Sa_pos(end).*(j_flux(i,N_neg+N_sep+1:N_tot)./D_s(N_neg+N_sep+1:N_tot))' - Sa_pos(end-1).*diff(C_s_shells_temp(:,end-1:end),1,2)./dR_pos)' ]';

            

            next_shell(1,:,:) = time_step.*[diff_t_C_s_neg' zeros(N_shell,N_sep) diff_t_C_s_pos'];

            C_s_full(i+1,:,:) = C_s_full(i,:,:) + next_shell;

            C_se(i+1,:) = squeeze(C_s_full(i+1,N_shell,:))'; %Set Surface Concentration
          
            C_s_ave(i+1,:) = squeeze(mean(C_s_full(i+1,:,:)))' ;



            %C_s_shells_neg_next = C_s_shells_temp(end,:);
       end


       if order == 6

        for k = 1 : N_tot_active
            C_s_shell = (C_s_full(i,:,k));

            if k <N_neg+1

                op = odeset('reltol',tol);
                [t,C_s_shell_temp] = ode15s(@(t,C_s_shell_temp) (D_s(k)./dV_neg').*[Sa_neg(1)*diff(C_s_shell_temp(1:2))./dR_neg ;diff(Sa_neg(1:end-1)'.*diff(C_s_shell_temp)./dR_neg) ; -Sa_neg(end)*j_flux(i,k)./D_s(k) - Sa_neg(end-1)*diff(C_s_shell_temp(end-1:end))./dR_neg ],[0 time_step], C_s_shell',op);
                
                C_s_shell_next = C_s_shell_temp(end,:);

                C_se(i+1,k) = C_s_shell_next(end); %Set Surface Concentration
                C_s_ave(i,k) = mean(C_s_shell_next);
                C_s_full(i+1,:,k) = C_s_shell_next;
            end


            if k>N_neg+N_sep
                [t,C_s_shell_temp] = ode15s(@(t,C_s_shell_temp) (D_s(k)./dV_pos').*[Sa_pos(1)*diff(C_s_shell_temp(1:2))./dR_pos ;diff(Sa_pos(1:end-1)'.*diff(C_s_shell_temp)./dR_pos) ; -Sa_pos(end)*j_flux(i,k)./D_s(k) - Sa_pos(end-1)*diff(C_s_shell_temp(end-1:end))./dR_pos ],[0 time_step], C_s_shell',op);
                
                C_s_shell_next = C_s_shell_temp(end,:);

                C_se(i+1,k) = C_s_shell_next(end); %Set Surface Concentration
                C_s_ave(i,k) = mean(C_s_shell_next);
                C_s_full(i+1,:,k) = C_s_shell_next;

            end

        end
       end

       if order == 7

               %[t,C_s_ave_temp] = ode15s(@(t,C_s_ave_temp) ((-3./Rs).*j_flux(i,:))',[0 time_step], C_s_ave(i,:)',op);

               diff_C_s_ave_change = ((-3./Rs).*j_flux(i,:));
               C_s_ave(i+1,:) = C_s_ave(i,:) + time_step*diff_C_s_ave_change;

               %[t,C_se_temp] = ode15s(@(t,C_se_temp) ( C_s_ave(i+1,:) - ( (Rs./(5*D_s)).*j_flux(i,:)) )',[0 time_step], C_se(i,:)',op);


               diff_C_se_change = ( C_s_ave(i+1,:) - ( (Rs./(5*D_s)).*j_flux(i,:)) );
                C_se(i+1,:) = C_se(i,:) + time_step*diff_C_se_change;

       end

       %% Verlet Time Integration For C_e and Temperature
       if i >2

            %%% For C_e
                
            C_e_accel = (1/time_step)*(C_e(i+1,:) - 2*C_e(i,:) + C_e(i-1,:)) ;

            C_e(i+1,:) = 2*C_e(i,:) - C_e(i-1,:) + (((time_step^2)/2))*C_e_accel;

            %%% For C_s
    
            if order >= 4
                C_s_full_accel = C_s_full(i+1,:,:) - 2*C_s_full(i,:,:) + C_s_full(i-1,:,:);
                
                C_s_full(i+1,:,:) = 2*C_s_full(i,:,:) - C_s_full(i-1,:,:) + ((time_step)^2).*C_s_full_accel;

                C_se_accel = (1/time_step)*(C_se(i+1,:) - 2*C_se(i,:) + C_se(i-1,:)) ;
    
                C_se(i+1,:) = 2*C_se(i,:) - C_se(i-1,:) + ((time_step)^2)*C_se_accel;

                C_s_ave_accel = (1/time_step)*(C_s_ave(i+1,:) - 2*C_s_ave(i,:) + C_s_ave(i-1,:)) ;

                C_s_ave(i+1,:) = 2*C_s_ave(i,:) - C_s_ave(i-1,:) + ((time_step)^2)*C_s_ave_accel;
            
            else

                C_se_accel = (1/time_step)*(C_se(i+1,:) - 2*C_se(i,:) + C_se(i-1,:)) ;
    
                C_se(i+1,:) = 2*C_se(i,:) - C_se(i-1,:) + ((time_step)^2)*C_se_accel;

                C_s_ave_accel = (1/time_step)*(C_s_ave(i+1,:) - 2*C_s_ave(i,:) + C_s_ave(i-1,:)) ;

                C_s_ave(i+1,:) = 2*C_s_ave(i,:) - C_s_ave(i-1,:) + ((time_step)^2)*C_s_ave_accel;



            end
       end

       %% Convergence Testing

       if all(C_e(i+1,:) >0) == 0
            
           Flag_convergence =3;

           break;

        end

        if all(C_se(i+1,1:N_neg) >0) == 0
            
           Flag_convergence =4;

           break;

        end

        if all(C_se(i+1,N_neg+N_sep+1:N_tot) >0) == 0
            
           Flag_convergence =5;

           break;

        end

        %% Calculate V_cell, capcaity calculations and time vector update
        V_cell(i) = psi_s(i,N_tot_active) - psi_s(i,1);
        %capacity_neg = 100*mean(stoic(i,N_neg+N_sep+1:N_tot_active));
        time = time + time_step;
        time_vector(i+1) = time;

        
        %% Simulation Stop Condition

        if sim_type == 1
            if (or(any(100*x< Cut_off_capacity),V_cell(i)<V_cut_off_low))

                break;
            end


        end

        if sim_type == 2
            if (or(any(100*y< Cut_off_capacity),V_cell(i)>V_cut_off_high))

                break;
            end

        end

        if sim_type == 2.5
            if (or(any(100*y< Cut_off_capacity),V_cell(i)>V_cut_off_high))

                I_app_CV = 1;
            end

        end

        if sim_type == 3
            if (or(any(100*y< Cut_off_capacity),and(V_cell(i)>V_cut_off_high,polarity_change ==1)))

                break;
            end

        end


    end

    %% Post Processing

    stop_time = i;

    V_cell = V_cell(1:stop_time);
    time_vector = time_vector(1,1:stop_time);
    C_e = C_e(1:stop_time,:);
    C_se = C_se(1:stop_time,:);
    j_flux = j_flux(1:stop_time,:);
    C_s_full = C_s_full(1:stop_time,:,:);
    I_app = I_app(1:stop_time);
    i_e = i_e(1:stop_time,:);
    i_s = i_s(1:stop_time,:);


end