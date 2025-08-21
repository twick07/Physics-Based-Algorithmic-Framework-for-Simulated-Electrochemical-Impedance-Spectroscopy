function [V_cell,Flag_convergence,time_vector,C_e,C_s_full,C_se,c_SEI,c_li_pl,I_app,psi_s,psi_e,j_tot] = P2D_MDN_function_ODE_int(Model_Parameters,N_tot,N_shell,order,param,time_step,sim_type,V_cut_off_low,V_cut_off_high,Cut_off_capacity,Ageing_models,C_e_init,C_s_full_init,C_se_init,c_SEI_init,c_li_init,I_app_init,C_rate,psi_s_init,psi_e_init,j_tot_init,eps_sol_init)
    
    %% Psi_Delta_function
    function [dpsi_deltadt] = diff_psi_delta_dt_calculator(t,psi_delta,OCP_temp,alpha,F,R_const,T,ex_current,kappa,c_e,I_app,sigma,del_mid,del_x,t_li,As,Cdl_temp,boundary_conds,SEI_Settings)
        
        %% Get Psi_delta

        psi_delta = psi_delta';

        %% Get SEI Parameters
        SEI_con = SEI_Settings(1,1);
        
        if SEI_con > 0 
            OCP_SEI_func = SEI_Settings(1,2);
    
            alpha_SEI_func = SEI_Settings(1,3);
    
            R_film_func = SEI_Settings(2,:);
    
            j_prev_func = SEI_Settings(3,:);
    
            i_o_SEI_func = SEI_Settings(4,:);

            %% Get Fluxes
    
            eta_SEI_func = psi_delta - OCP_SEI_func - F*R_film_func.*j_prev_func;
                        
            j_SEI_func = -( (i_o_SEI_func./F).*exp( -(alpha_SEI_func*F./(R_const.*T)).*eta_SEI_func  ) );

    
            j_flux_func = (2.*ex_current).*sinh((alpha.*F./(R_const.*T)).*(psi_delta - OCP_temp - F*R_film_func.*j_prev_func));
    
            j_sum = j_flux_func + j_SEI_func;
        else
            
            j_sum = (2.*ex_current).*sinh((alpha.*F./(R_const.*T)).*(psi_delta - OCP_temp));
        end
        
        C_e_log_lbc = boundary_conds(2,1);

        C_e_log_rbc = boundary_conds(2,2);

        diffx_C_e_log = [C_e_log_lbc diff(log(c_e),1,2)./del_mid(2:end-1) C_e_log_rbc];


        kappa_d = (2.*kappa.*R_const.*T*(1-t_li)./F);

        cond = (kappa.*sigma)./(kappa + sigma);
        psi_delta_lbc = boundary_conds(1,1);

        psi_delta_rbc = boundary_conds(1,2);

        
        diff_psi_delta_x = [psi_delta_lbc diff(psi_delta)./del_mid(2:end-1) psi_delta_rbc];

        i_e_temp = cond.*(diff_psi_delta_x + (kappa_d./kappa).*diffx_C_e_log + I_app./sigma);

        dpsi_deltadt = ((1./(As.*Cdl_temp)).*(diff(i_e_temp,1,2)./del_x - As.*F.*j_sum))';

        
    end

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

    %% Anode SEI Specific Parameters
    L_SEI_param = Model_Parameters(30);
    kappa_SEI = Model_Parameters(31);
    M_SEI = Model_Parameters(32);
    rho_SEI = Model_Parameters(33);
    i_0_SEI = Model_Parameters(34);
    OCP_SEI = Model_Parameters(35);
    alpha_SEI = Model_Parameters(36);
    c_sol_0 = Model_Parameters(43);
    k_SEI = Model_Parameters(44);
    D_SEI = Model_Parameters(45);
    


    %% Anode Specific Plating Parameters
    k_li = Model_Parameters(37);
    alpha_pl_neg = Model_Parameters(38);
    alpha_pl_pos = Model_Parameters(39);
    kappa_li = Model_Parameters(40);
    M_Li = Model_Parameters(41);
    rho_Li = Model_Parameters(42);

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

    if sim_type == 0.5
        I_app = C_rate*I_1c.*I_app_init;
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
    if sum(eps_sol_init) == 0
        eps_elec = [eps_neg_elec*ones(1,N_neg) eps_sep_elec*ones(1,N_sep) eps_pos_elec*ones(1,N_pos)];
        D_e_eff = [D_neg_e_eff*ones(1,N_neg) D_sep_e_eff*ones(1,N_sep) D_pos_e_eff*ones(1,N_pos)];
        D_e = D_e_eff;

    else
        eps_elec = eps_sol_init;        
        D_e_eff(1:N_neg) = (eps_elec(1:N_neg).^brugg(1:N_neg)).*D_elec;
        D_e = D_e_eff;
    end
    eps_sol = [eps_neg*ones(1,N_neg) zeros(1,N_sep) eps_pos*ones(1,N_pos)];
    As = [As_neg*ones(1,N_neg) zeros(1,N_sep) As_pos*ones(1,N_pos)];
    Rs = [Rs_neg*ones(1,N_neg) zeros(1,N_sep) Rs_pos*ones(1,N_pos)];
    T = T_amb*ones(time_max,N_tot);
    sigma = [sigma_neg_eff*ones(1,N_neg) zeros(1,N_sep) sigma_pos_eff*ones(1,N_pos)];
    kappa_SEI_v = [kappa_SEI*ones(1,N_neg) zeros(1,N_sep) zeros(1,N_pos)];
    M_SEI_v = [M_SEI*ones(1,N_neg) zeros(1,N_sep) zeros(1,N_pos)];
    rho_SEI_v = [rho_SEI*ones(1,N_neg) zeros(1,N_sep) zeros(1,N_pos)];
    i_o_SEI = [i_0_SEI*ones(1,N_neg) zeros(1,N_sep) zeros(1,N_pos)];
    kappa_pl_v = [kappa_li*ones(1,N_neg) zeros(1,N_sep) zeros(1,N_pos)];
    M_li_v = [M_Li*ones(1,N_neg) zeros(1,N_sep) zeros(1,N_pos)];
    rho_li_v = [rho_Li*ones(1,N_neg) zeros(1,N_sep) zeros(1,N_pos)];
    i_o_li = [zeros(1,N_neg) zeros(1,N_sep) zeros(1,N_pos)];
    

    %% Dependent Variables for P2D-MDN
    V_cell = zeros(1,time_max+1);
    j_flux = zeros(time_max+1,N_tot_active);
    j_tot = zeros(time_max+1,N_tot_active);
    if sum(j_tot_init) == 0
        j_tot_init = j_tot(1,:);

    else

        j_tot(1,:) = j_tot_init;
    end
    
    j_dl  = zeros(time_max+1,N_tot_active);
    j_SEI = zeros(time_max+1,N_tot_active);
    j_li = zeros(time_max+1,N_tot_active);
    psi_e = zeros(time_max+1,N_tot_active);
    psi_s = zeros(time_max+1,N_tot_active);
    psi_delta = zeros(time_max+1,N_tot_active);
    eta = zeros(time_max+1,N_tot_active);
    i_e = zeros(time_max+1,N_tot_active+1);
    i_s = zeros(time_max+1,N_tot_active+1);
    ex_current = zeros(time_max+1,N_tot_active);
    OCP = zeros(time_max+1,N_tot_active);
    stoic = zeros(time_max+1,N_tot_active);
    
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
    
    %% Dependent Variables for SEI Model
    c_SEI = zeros(time_max+1,N_tot_active);
    L_SEI = zeros(time_max+1,N_tot_active);
    R_SEI = zeros(time_max+1,N_tot_active);

    if sum(c_SEI_init) == 0
        
        c_SEI(1,:) = [(L_SEI_param*rho_SEI/M_SEI).*As(1:N_neg) zeros(1,N_sep) zeros(1,N_pos)];
    else
        c_SEI(1,:) = c_SEI_init;
    end
    L_SEI(1,1:N_neg) = c_SEI(1,1:N_neg).*M_SEI_v(1:N_neg)./(rho_SEI_v(1:N_neg).*As(1:N_neg));
    
    R_SEI(1,1:N_neg) = L_SEI(1,1:N_neg)./kappa_SEI_v(1:N_neg);

    %% Dependent Variables for Li plating Model
    
    L_pl = zeros(time_max+1,N_tot_active);
    R_li = zeros(time_max+1,N_tot_active);
    if sum(c_li_init) == 0
        c_li_pl = zeros(time_max+1,N_tot_active);
    else
        c_li_pl = zeros(time_max+1,N_tot_active);
        c_li_pl(1,1:N_neg) = c_li_init(1:N_neg);
    end
    L_pl = zeros(time_max+1,N_tot_active);
    L_pl(1,1:N_neg) = c_li_pl(1,1:N_neg).*M_li_v(1:N_neg)./(rho_li_v(1:N_neg).*As(1:N_neg));
    
    R_li = zeros(time_max,N_tot_active);
    R_li(1,1:N_neg) = L_pl(1,1:N_neg)./kappa_pl_v(1:N_neg);
    
    R_film = zeros(time_max+1,N_tot_active);
    R_film(1,1:N_neg) = R_SEI(1,1:N_neg) + R_li(1,1:N_neg);

    L_film = zeros(time_max+1,N_tot_active);
    L_film(1,1:N_neg) = L_SEI(1,1:N_neg) + L_pl(1,1:N_neg);

    

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
    case_flux = 5;
    Convergence_max = 5000;
    %SEI = 1;
    %li_pl = 1;
    SEI = Ageing_models(1);
    li_pl = Ageing_models(2);
    OCP_li = 0; % OCP_li is li vs li+ which is fixed to 0
    %OCP_li = 0.14; % OCP_li is li vs li+ which is fixed to 0
    coupling = Ageing_models(3);

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
        
        %Calculate Plating Exchange current
        if li_pl == 1
            i_o_li(i,1:N_neg) = F*k_li;
        end
        if li_pl == 2
            i_o_li(i,1:N_neg) = F*k_li*(C_e_reduced(1:N_neg).^alpha_pl_neg);
        end

        %Set initial values
        if and(sum(psi_s_init) == 0 ,i ==1)
            psi_s(i,:) = OCP(i,:);
            psi_delta(i,:) = OCP(i,:);
            V_cell(i) = psi_s(i,end) - psi_s(i,1);

        end

        if and(sum(psi_s_init(end,:)) ~= 0 ,i ==1)
            psi_s(i,:) = psi_s_init(end,:);
            psi_e(i,:) = psi_e_init(end,:);
            psi_delta(i,:) = psi_s(i,:)-psi_e(i,:);
            V_cell(i) = psi_s(i,end) - psi_s(i,1);

        end

        %% psi_delta_calculator for negative electrode
        % Make extended parameters

        A = 0.5*eye(N_neg+1,N_neg);
            
        B = A + circshift(A,[1 0]);
        
        del_mid = B*del(1:N_neg)';
        
        del_mid = del_mid';

        del_mid(1) = del_mid(2);

        del_mid(end) = (del(N_neg)+del(N_neg+1))/2 ;

        kappa_n_extended = B*kappa(i,1:N_neg)';

        kappa_n_extended = kappa_n_extended';

        kappa_n_extended(1) = kappa_n_extended(2);

        kappa_n_extended(end) = (kappa(i,N_neg)+kappa(i,N_neg+1))/2 ;

        conductivity_n = (kappa_n_extended.*sigma(1))./(kappa_n_extended + sigma(1));

        psi_delta_lbc_def = -I_app(i)/sigma(1);

        kappa_n_s = (kappa(i,N_neg)+kappa(i,N_neg+1))/2;

        del_n_s = (del_neg + del_sep)/2 ;
        
        ln_c_e_n_s = diff(log(C_e_reduced(N_neg:N_neg+1)))./del_n_s ;

        %psi_delta_rbc_def = (1/conductivity_n(end)).*I_app(i) - I_app(i)./sigma(1) - (2.*kappa_n_s.*R_const.*T_amb*(1-t_li)./F).*ln_c_e_n_s;
        psi_delta_rbc_def = (1/conductivity_n(end)).*I_app(i) - I_app(i)./sigma(1) - (2.*R_const.*T_amb*(1-t_li)./F).*ln_c_e_n_s;
        %psi_delta_rbc_def = (1/conductivity_n(end)).*I_app(i) - I_app(i)./sigma(1) - (2.*R_const.*T_amb*(1-t_li)./F).*ln_c_e_n_s;

        %psi_delta_rbc_def = I_app(i)./kappa_n_extended(end) - (2.*kappa_n_s.*R_const.*T_amb*(1-t_li)./F).*ln_c_e_n_s;

        %grad_psi_delta_right_n = (1/pybamm.boundary_value(conductivity_n,"right"))*I_app - I_app/sigma_eff_n - pybamm.boundary_value(kappa_d_n,"right")*(1/pybamm.boundary_value(c_e_n,"right"))*grad_c_e_n_right



        ln_c_e_n_0 = 0;


        %i_e_lbc_def = 0;

        %i_e_rbc_def = I_app(i);

        boundary_conds =[psi_delta_lbc_def psi_delta_rbc_def; ln_c_e_n_0 ln_c_e_n_s];
        
        SEI_Settings = zeros(4,N_neg);
        
        SEI_Settings(1,1) = SEI*(I_app(i)<0);

        SEI_Settings(1,2) = OCP_SEI;
    
        SEI_Settings(1,3) = alpha_SEI ;

        SEI_Settings(2,:) = R_film(i,1:N_neg);

        j_ave_neg = I_app(i)/(As_neg*F*L_neg)  ;
        %SEI_Settings(3,:) = j_tot(i,1:N_neg);
        SEI_Settings(3,:) = j_ave_neg.*ones(1,N_neg);
        SEI_Settings(4,:) = i_o_SEI(1:N_neg) ;
        


        %[t,C_e_temp] = ode15s(@(t,C_e_temp) ( (1./eps_elec').*( ((1./del).*diff((D_e_eff_ext.*[ 0 ; diff(C_e_temp,1) ; 0 ]')./del_extended,1))' + ((1-t_li)*(As).*(j_tot(i,:)))' )), [0 time_step], C_e(i,:)', op);  %% Changed
        op = odeset('reltol',tol);

        %diff_psi_delta_dt_calculator(t,psi_delta,OCP,alpha,F,R_const,T,ex_current,cond,kappa,c_e,I_app,sigma,del_mid,t_li,As,Cdl_temp,boundary_conds)
        [t,psi_delta_temp_n] = ode15s(@(t,psi_delta_temp_n)  diff_psi_delta_dt_calculator(t,psi_delta_temp_n,OCP(i,1:N_neg),alpha,F,R_const,T_amb,ex_current(i,1:N_neg),kappa_n_extended,C_e_reduced(1:N_neg),I_app(i),sigma(1),del_mid(1:N_neg+1),del(1:N_neg),t_li,As(1:N_neg),Cdl(1:N_neg),boundary_conds,SEI_Settings), [0 time_step], psi_delta(i,1:N_neg)', op);

        psi_delta(i+1,1:N_neg) = psi_delta_temp_n(end,:); % update


        %% psi_delta_calculator for positive electrode
        % Make extended parameters

        A = 0.5*eye(N_pos+1,N_pos);
            
        B = A + circshift(A,[1 0]);
        
        del_mid = B*del(N_neg+N_sep+1:N_tot)';
        
        del_mid = del_mid';

        del_mid(1) = (del(N_neg+N_sep)+del(N_neg+N_sep+1))/2;

        del_mid(end) = del_mid(end-1) ;

        kappa_p_extended = B*kappa(i,N_neg+N_sep+1:N_tot)';

        kappa_p_extended = kappa_p_extended';

        kappa_p_extended(1) = (kappa(i,N_neg+N_sep)+kappa(i,N_neg+N_sep+1))/2 ;

        kappa_p_extended(end) = kappa_p_extended(end-1) ;

        conductivity_p = (kappa_p_extended.*sigma(end))./(kappa_p_extended + sigma(end));

        kappa_s_p = (kappa(i,N_neg+N_sep)+kappa(i,N_neg+N_sep+1))/2;

        del_s_p = (del_pos + del_sep)/2 ;
        
        ln_c_e_s_p = diff(log(C_e_reduced(N_neg+N_sep:N_neg+N_sep+1)))./del_s_p ;

        %psi_delta_lbc_def = (1/conductivity_p(1)).*I_app(i) - I_app(i)./sigma(end) - (2.*kappa_s_p.*R_const.*T_amb*(1-t_li)./F).*ln_c_e_s_p;
        psi_delta_lbc_def = (1/conductivity_p(1)).*I_app(i) - I_app(i)./sigma(end) - (2.*R_const.*T_amb*(1-t_li)./F).*ln_c_e_s_p;

        %psi_delta_lbc_def = I_app(i)./kappa_p_extended(1) - (2.*kappa_s_p.*R_const.*T_amb*(1-t_li)./F).*ln_c_e_s_p;

        psi_delta_rbc_def = -I_app(i)/sigma(end);
        %grad_psi_delta_right_n = (1/pybamm.boundary_value(conductivity_n,"right"))*I_app - I_app/sigma_eff_n - pybamm.boundary_value(kappa_d_n,"right")*(1/pybamm.boundary_value(c_e_n,"right"))*grad_c_e_n_right


        
        %i_e_lbc_def = 0;

        %i_e_rbc_def = I_app(i);

        ln_c_e_p_end = 0;

        boundary_conds =[psi_delta_lbc_def psi_delta_rbc_def; ln_c_e_s_p ln_c_e_p_end];

        
        SEI_Settings = zeros(4,N_pos);

        %[t,C_e_temp] = ode15s(@(t,C_e_temp) ( (1./eps_elec').*( ((1./del).*diff((D_e_eff_ext.*[ 0 ; diff(C_e_temp,1) ; 0 ]')./del_extended,1))' + ((1-t_li)*(As).*(j_tot(i,:)))' )), [0 time_step], C_e(i,:)', op);  %% Changed
        op = odeset('reltol',tol);

        %diff_psi_delta_dt_calculator(t,psi_delta,OCP,alpha,F,R_const,T,ex_current,cond,kappa,c_e,I_app,sigma,del_mid,t_li,As,Cdl_temp,boundary_conds)
        [t,psi_delta_temp_p] = ode15s(@(t,psi_delta_temp_p)  diff_psi_delta_dt_calculator(t,psi_delta_temp_p,OCP(i,N_neg+N_sep+1:N_tot),alpha,F,R_const,T_amb,ex_current(i,N_neg+N_sep+1:N_tot),kappa_p_extended,C_e_reduced(N_neg+N_sep+1:N_tot),I_app(i),sigma(end),del_mid(1:N_pos+1),del(N_neg+N_sep+1:N_tot),t_li,As(N_neg+N_sep+1:N_tot),Cdl(N_neg+N_sep+1:N_tot),boundary_conds,SEI_Settings), [0 time_step], psi_delta(i,N_neg+N_sep+1:N_tot)', op);

        psi_delta(i+1,N_neg+N_sep+1:N_tot) = psi_delta_temp_p(end,:); % update


        %% Calculate Fluxes
        if SEI == 0
            eta(i+1,:) = psi_delta(i+1,:) - OCP(i,:);
    
            j_flux(i+1,:) = 2.*ex_current(i,:).*sinh( ( alpha.*F./(R_const.*T(i,:)) ).*( eta(i+1,:) ) );
    
            j_dl(i+1,:) = (Cdl./F).*((psi_delta(i+1,:)-psi_delta(i,:))./time_step);
    
            j_tot(i+1,:) = j_flux(i+1,:) + j_dl(i+1,:);

        else

            %eta(i+1,:) = psi_delta(i+1,:) - OCP(i,:) - F.*R_film(i,:).*j_tot(i,:);
            eta(i+1,:) = psi_delta(i+1,:) - OCP(i,:) - F.*R_film(i,:).*j_ave_neg;
            
            %eta_SEI = psi_delta(i+1,:) - OCP_SEI - F.*R_film(i,:).*j_tot(i,:);
            eta_SEI = psi_delta(i+1,:) - OCP_SEI - F.*R_film(i,:).*j_ave_neg;

            j_flux(i+1,:) = 2.*ex_current(i,:).*sinh( ( alpha.*F./(R_const.*T(i,:)) ).*( eta(i+1,:) ) );
    
            j_dl(i+1,:) = (Cdl./F).*((psi_delta(i+1,:)-psi_delta(i,:))./time_step);

            j_SEI(i+1,:) = -(I_app(i)<0).*( (i_o_SEI./F).*exp( -(alpha_SEI*F./(R_const.*T(i,:))).*eta_SEI  ) );

            j_tot(i+1,:) = j_flux(i+1,:) + j_dl(i+1,:) + j_SEI(i+1,:);

        end


        %% Calculate Potentials

        % Calculate psi_e
        A = 0.5*eye(N_tot_active+1,N_tot_active);
            
        B = A + circshift(A,[1 0]);
        
        del_mid = B*del';
        
        del_mid = del_mid';
        
        diff_C_e_log = diff(log(C_e_reduced),1,2) ;

        diff_C_e_log = [0 diff_C_e_log 0]./del_mid;


        psi_e(i+1,1) = 0;
        diff_psi_e_m(i+1,1) = 0;
        for j = 1:N_tot_active-1
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


            diff_psi_e_m(i+1,j+1) = (1/kappa_k_plus_half)*((kappa_k_minus_half)*diff_psi_e_m(i+1,j) + (kd_k_plus_half*diff_C_e_log(j+1) - kd_k_minus_half*diff_C_e_log(j)) - (F*As(j)*del(j)*j_tot(i+1,j)));
    
            

            psi_e(i+1,j+1) = psi_e(i+1,j) + del(j)*diff_psi_e_m(i+1,j+1);
        
        
                

        end

        psi_s(i+1,:) = psi_delta(i+1,:) +psi_e(i+1,:);

        psi_s(i+1,N_neg+1:N_neg+N_sep) = zeros(1,N_sep);


       %% Calculate Concentrations
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
       
       [t,C_e_temp] = ode15s(@(t,C_e_temp) ( (1./eps_elec').*( ((1./del).*diff((D_e_eff_ext.*[ 0 ; diff(C_e_temp,1) ; 0 ]')./del_extended,1))' + ((1-t_li)*(As).*(j_tot(i+1,:)))' )), [0 time_step], C_e(i,:)', op);  %% Changed

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
                        (diff(diag(Sa_neg(1:end-1))'*(diff(C_s_shells_temp,1,2)./dR_neg)')) ; (-Sa_neg(end).*(j_flux(i+1,1:N_neg)./D_s(1:N_neg))' - Sa_neg(end-1).*diff(C_s_shells_temp(:,end-1:end),1,2)./dR_neg)' ]';

            

            C_s_shells_temp  = C_s_shells(:,N_neg+N_sep+1:N_tot)';

            diff_t_C_s_pos = ((D_s(N_neg+N_sep+1:N_tot)'*(1./dV_pos))).*[(Sa_pos(1).*diff(C_s_shells_temp(:,1:2),1,2)./dR_pos)' ;...
                        (diff(diag(Sa_pos(1:end-1))'*(diff(C_s_shells_temp,1,2)./dR_pos)')) ; (-Sa_pos(end).*(j_flux(i+1,N_neg+N_sep+1:N_tot)./D_s(N_neg+N_sep+1:N_tot))' - Sa_pos(end-1).*diff(C_s_shells_temp(:,end-1:end),1,2)./dR_pos)' ]';

            

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

       %% Calculate SEI growth
        %%% Kinetic Limited SEI Growth
        if SEI ==0
            %%% Add SEI Growth
            c_SEI(i+1,1:N_neg) = c_SEI(i,1:N_neg) - time_step*As(1:N_neg).*j_SEI(i+1,1:N_neg);
            
            L_SEI(i+1,1:N_neg) = (1./As(1:N_neg)).*(c_SEI(i+1,1:N_neg).*M_SEI_v(1:N_neg)./rho_SEI_v(1:N_neg));
            
            R_SEI(i+1,1:N_neg)= L_SEI(i+1,1:N_neg)./kappa_SEI_v(1:N_neg);
        end

        if li_pl == 0

            c_li_pl(i+1,:) = c_li_pl(i,:) - time_step*As.*j_li(i+1,:);

            L_pl(i+1,1:N_neg) = (1./As(1:N_neg)).*(c_li_pl(i+1,1:N_neg).*M_li_v(1:N_neg)./rho_li_v(1:N_neg));

            R_li(i+1,1:N_neg)= L_pl(i+1,1:N_neg)./kappa_pl_v(1:N_neg);

        end

        if coupling == 0

            R_film(i+1,1:N_neg) = R_SEI(i+1,1:N_neg) + R_li(i+1,1:N_neg);

            L_film(i+1,1:N_neg) = L_SEI(i+1,1:N_neg) + L_pl(i+1,1:N_neg);

        end


        if SEI ==1
            %%% Add SEI Growth
            c_SEI(i+1,1:N_neg) = c_SEI(i,1:N_neg) - time_step*As(1:N_neg).*j_SEI(i+1,1:N_neg);
            
            L_SEI(i+1,1:N_neg) = (1./As(1:N_neg)).*(c_SEI(i+1,1:N_neg).*M_SEI_v(1:N_neg)./rho_SEI_v(1:N_neg));
            
            R_SEI(i+1,1:N_neg)= L_SEI(i+1,1:N_neg)./kappa_SEI_v(1:N_neg);
        end

        if li_pl == 1

            c_li_pl(i+1,:) = c_li_pl(i,:) - time_step*As.*j_li(i+1,:);

            L_pl(i+1,1:N_neg) = (1./As(1:N_neg)).*(c_li_pl(i+1,1:N_neg).*M_li_v(1:N_neg)./rho_li_v(1:N_neg));

            R_li(i+1,1:N_neg)= L_pl(i+1,1:N_neg)./kappa_pl_v(1:N_neg);

        end

        if coupling == 1

            R_film(i+1,1:N_neg) = R_SEI(i+1,1:N_neg) + R_li(i+1,1:N_neg);

            L_film(i+1,1:N_neg) = L_SEI(i+1,1:N_neg) + L_pl(i+1,1:N_neg);

        end

        if coupling == 2
            
             %% Update Degredation Parameters

            R_film(i+1,1:N_neg) = R_SEI(i+1,1:N_neg) + R_li(i+1,1:N_neg);

            L_film(i+1,1:N_neg) = L_SEI(i+1,1:N_neg) + L_pl(i+1,1:N_neg);

            %eps_sol(1:N_neg) = eps_sol(1:N_neg) - As(1:N_neg).*(L_film(i+1,1:N_neg) - L_film(i,1:N_neg));
            eps_elec(1:N_neg) = eps_elec(1:N_neg)- As(1:N_neg).*(L_film(i+1,1:N_neg) - L_film(i,1:N_neg));

            %eps_sol_store(i+1,1:N_neg) = eps_elec(1:N_neg);

           
            %eps_sol_store(i+1,1:N_neg) = eps_sol(1:N_neg);

            %%  Update all internal parameters
            D_e_eff(1:N_neg) = (eps_elec(1:N_neg).^brugg(1:N_neg)).*D_elec;
            D_e = D_e_eff;


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
        V_cell(i+1) = psi_s(i+1,N_tot_active) - psi_s(i+1,1);
        %capacity_neg = 100*mean(stoic(i,N_neg+N_sep+1:N_tot_active));
        time = time + time_step;
        time_vector(i+1) = time;

        
        %% Simulation Stop Condition

        if sim_type == 1
            if (or(any(100*x< Cut_off_capacity),V_cell(i+1)<V_cut_off_low))

                break;
            end


        end

        if sim_type == 2
            if (or(any(100*y< Cut_off_capacity),V_cell(i+1)>V_cut_off_high))

                break;
            end

        end

        if sim_type == 2.5
            if (or(any(100*y< Cut_off_capacity),V_cell(i+1)>V_cut_off_high))

                I_app_CV = 1;
            end

        end

        if sim_type == 3
            if (or(any(100*y< Cut_off_capacity),and(V_cell(i+1)>V_cut_off_high,polarity_change ==1)))

                break;
            end

        end

    end

   
    %% Post Processing

    stop_time = i;
    %V_cell_init = psi_s_init(end) - psi_s_init(1);
    %V_cell = [V_cell_init V_cell(1:stop_time)];
    V_cell = V_cell(1,1:stop_time+1);
    time_vector = time_vector(1,1:stop_time+1);
    C_e = C_e(1:stop_time+1,:);
    C_se = C_se(1:stop_time+1,:);
    C_s_full = C_s_full(1:stop_time+1,:,:);
    %j_tot = [j_tot_init; j_tot(1:stop_time,:)];
    j_tot = j_tot(i:stop_time+1,:);
    psi_s = psi_s(1:stop_time+1,:);
    psi_e = psi_e(1:stop_time+1,:);
    %psi_s = [psi_s_init;psi_s(1:stop_time,:)];
    %psi_e = [psi_e_init;psi_e(1:stop_time,:)];
    %if SEI > 0
    c_SEI = c_SEI(1:stop_time+1,:);
    %end

    %if li_pl > 0
    c_li_pl = c_li_pl(1:stop_time+1,:);
    %end
    if sum(j_tot_init) ==0
        I_app_t_0 = 0;
        I_app = [I_app_t_0 I_app(1:stop_time)];
    else
        I_app = [I_app(1:stop_time) I_app(stop_time)];
    end

end