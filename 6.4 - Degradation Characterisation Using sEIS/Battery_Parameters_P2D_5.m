if param == 1

    % From Marquis et al (LCO/C6) Cell
    %% Cell properties and constants
    current_collector = 0;  % Model has no current collectors
    A = 1;  % [m^2]
    
    % Constants
    R_const = 8.314;  % Ideal gas constant [J/(K*mol)]
    F = 96487;  % Faraday Constant [C/mol]
    %T_c = 25.15;  % Temperature in Celsius
    T_c = 25.15;  % Temperature in Celsius
    T = 273 + T_c;  % Temperature in Kelvin
    T_amb = T;
    
    % Bruggman constants
    brugg = 1.5;
    brugg_neg = brugg;
    brugg_sep = brugg;
    brugg_pos = brugg;
    
    % Cell dimensions
    L_ncc = 0;
    L_neg = 0.0001;  % [m]
    L_sep = 2.5e-05;  % [m]
    L_pos = 0.0001;  % [m]
    L_pcc = 0;
    L_tot_active = L_neg + L_sep + L_pos;  % [m]
    L_tot = L_ncc + L_tot_active + L_pcc;  % [m]
    
    %% Separator/Electrolyte Parameters
    C_e_initial = 1000.0;  % Initial electrolyte lithium concentration [mol/m^3]
    D_elec = 5.34e-10;  % Diffusion Coefficient for Lithium in electrolyte [m^2/s] 
    D_sep_e_eff = D_elec;
    t_li = 0.363;  % Lithium Ion Transference Number [unitless]
    eps_sep_elec = 0.4;  % Electrolyte volume fraction
    kappa_const = 0.6013; % Function Lookup
    dlnfdlnce_sep = 0;
    
    %% Anode Parameters
    %soc_initial_neg = 0.53;
    soc_initial_neg = 0.8;
    soc_initial_neg_100 = 0.8;
    soc_initial_neg_80 = 0.7948;
    soc_initial_neg_60 = 0.6722;
    soc_initial_neg_40 = 0.6181;
    soc_initial_neg = soc_initial_neg_100;
    eps_neg_elec = 0.3;  % Electrolyte volume fraction (porosity)
    eps_neg = 0.6; % Solid volume fraction
    C_max_neg = 24983.2619938437;  % Maximum solid lithium concentration [mol/m^3]
    C_neg_initial = soc_initial_neg * C_max_neg;
    Rs_neg = 10e-6;  % Electrode particle size [m]
    D_neg_s = 3.9e-14;  % Solid diffusion coefficient [m^2/s] 
    D_neg_eff_s = D_neg_s;
    D_neg_e_eff = D_elec;
    sigma_neg = 100;  % Electrical conductivity [S/m]
    m_ref_neg = 2e-5;
    %K_0_neg = 1.94e-11;  % Negative rate constant #Function
    K_0_neg = m_ref_neg/F ;
    As_neg = 3 * (eps_neg / Rs_neg);  % Surface area to volume ratio [m^-1]
    %Cdl_neg = 100;
    Cdl_neg = 0.2;
    
    
    
    %% Cathode Parameters
    %soc_initial_pos = 0.17;
    soc_initial_pos = 0.6;
    soc_initial_neg_100 = 0.6;
    soc_initial_neg_80 = 0.5877;
    soc_initial_neg_60 = 0.6427;
    soc_initial_neg_40 = 0.6841;
    soc_initial_pos = soc_initial_neg_100;
    eps_pos_elec = 0.3;  % Electrode Porosity
    eps_pos = 0.5;  % Solid volume fraction
    C_max_pos = 51217.9257309275;  % Maximum solid lithium concentration [mol/m^3]
    C_pos_initial = soc_initial_pos * C_max_pos;
    Rs_pos = 10e-6;  % Electrode particle size [m]
    D_pos_s = 1e-13;  % Solid diffusion coefficient [m^2/s] 
    D_pos_eff_s = D_pos_s;
    D_pos_e_eff = D_elec;
    sigma_pos = 10;  % Electrical conductivity [S/m]
    m_ref_pos = 6e-7;
    %K_0_pos = 2.16e-11;  % Cathode rate constant #Function
    K_0_pos = m_ref_pos/F;
    As_pos = 3 * (eps_pos / Rs_pos);  % Surface area to volume ratio [m^-1]
    %Cdl_pos = 0;
    Cdl_pos = 0.2;

    %% Anode SEI Specific Parameters
    L_SEI = 1e-9;%%"Initial SEI Thickness [m]"
    kappa_SEI = 5e-6 ;%"SEI ionic conductivity [S/m]"
    M_SEI = 0.162 ;%"SEI Molar Mass [kg/mol]"
    rho_SEI = 1.69e3 ;% "SEI layer density [kg/m^3]"
    i_0_SEI_tune = 1.02;
    i_0_SEI = i_0_SEI_tune*2e-03;
    OCP_SEI = 0 ;
    alpha_SEI = 0.5 ;%"Charge transfer coefficient, SEI"
    kappa_li = 1/(92.8e-9); % Lithium Conductivity [S/m], taken from wikipedia
    M_Li = 6.94e-3; % Molar Mass of Lithium from RSC website [Kg/mol]
    rho_Li = 0.534e3; % Density of Lithiun from RSC website [Kg/m^3]    
    k_0_SEI_tune = 1.02;    
    k_SEI = k_0_SEI_tune*1.0e-11; %"Rate constant adjusted, SEI [m/s]"
    c_sol_0 = 4541; % Concentration of Solvent in bulk electrolyte [mol/m^3]
    D_SEI =  1e-21 ;% "EC diffusion coefficient [m^2/s]"
    


    %% Anode Plating Specific Parameters
    k_li_tuning = 1.28;
    k_li = k_li_tuning*6.0e-10; % Lithium plating rate constant
    alpha_pl_neg = 0.3; % Anodic Plating Transfer Coefficient
    alpha_pl_pos = 0.7; % Cathodic Plating Transfer Coefficient
    OCP_li =0.165; % OCP Plating Potential
    


end

if param == 2

    % From Marquis et al (LCO/C6) Cell
    %% Cell properties and constants
    current_collector = 0;  % Model has no current collectors
    A = 1;  % [m^2]
    
    % Constants
    R_const = 8.314;  % Ideal gas constant [J/(K*mol)]
    F = 96487;  % Faraday Constant [C/mol]
    T_c = 25.15;  % Temperature in Celsius
    T = 273 + T_c;  % Temperature in Kelvin
    T_amb = T;
    
    % Bruggman constants
    brugg = 1.5;
    brugg_neg = brugg;
    brugg_sep = brugg;
    brugg_pos = brugg;
    
    % Cell dimensions
    
    L_neg = 100e-6;  % [m]
    L_sep = 52e-6;  % [m]
    L_pos = 183e-6;  % [m]
    
    L_tot_active = L_neg + L_sep + L_pos;  % [m]
    L_tot =  L_tot_active;  % [m]
    
    %% Separator/Electrolyte Parameters
    C_e_initial = 2000.0;  % Initial electrolyte lithium concentration [mol/m^3]
    D_elec = 7.5e-11;  % Diffusion Coefficient for Lithium in electrolyte [m^2/s] 
    D_sep_e_eff = D_elec;
    t_li = 0.363;  % Lithium Ion Transference Number [unitless]
    eps_sep_elec = 1.0;  % Electrolyte volume fraction
    kappa_const = 0.6013; % Function Lookup
    dlnfdlnce_sep = 0;
    
    %% Anode Parameters
    %soc_initial_neg = 0.53;
    soc_initial_neg = 0.58;
    soc_initial_neg_100 = 0.58;
    soc_initial_neg_80 = 0.7948;
    soc_initial_neg_60 = 0.6722;
    soc_initial_neg_40 = 0.6181;
    soc_initial_neg = soc_initial_neg_100;
    eps_neg_elec = 0.357;  % Electrolyte volume fraction (porosity)
    eps_neg = 0.471; % Solid volume fraction
    C_max_neg = 26390;  % Maximum solid lithium concentration [mol/m^3]
    C_neg_initial = soc_initial_neg * C_max_neg;
    Rs_neg = 12.5e-6;  % Electrode particle size [m]
    D_neg_s = 3.9e-14;  % Solid diffusion coefficient [m^2/s] 
    D_neg_eff_s = D_neg_s;
    D_neg_e_eff = D_elec;
    sigma_neg = 100;  % Electrical conductivity [S/m]
    K_0_neg = 2e-10 ;
    As_neg = 3 * (eps_neg / Rs_neg);  % Surface area to volume ratio [m^-1]
    Cdl_neg = 0.2;
    
    
    %% Cathode Parameters
    soc_initial_pos = 0.19;
    soc_initial_neg_100 = 0.19;
    soc_initial_neg_80 = 0.5877;
    soc_initial_neg_60 = 0.6427;
    soc_initial_neg_40 = 0.6841;
    soc_initial_pos = soc_initial_neg_100;
    eps_pos_elec = 0.4440;  % Electrode Porosity
    eps_pos = 0.2970;  % Solid volume fraction
    C_max_pos = 22860;  % Maximum solid lithium concentration [mol/m^3]
    C_pos_initial = soc_initial_pos * C_max_pos;
    Rs_pos = 8.5e-6;  % Electrode particle size [m]
    D_pos_s = 1e-13;  % Solid diffusion coefficient [m^2/s] 
    D_pos_eff_s = D_pos_s;
    D_pos_e_eff = D_elec;
    sigma_pos = 3.8;  % Electrical conductivity [S/m]
    K_0_pos = 2.0e-10;
    As_pos = 3 * (eps_pos / Rs_pos);  % Surface area to volume ratio [m^-1]
    Cdl_pos = 0.2;

    %% Anode SEI Specific Parameters
    L_SEI =  1e-9;%%"Initial SEI Thickness [m]"
    kappa_SEI = 5e-6 ;%"SEI ionic conductivity [S/m]"
    M_SEI = 0.162 ;%"SEI Molar Mass [kg/mol]"
    rho_SEI = 1.69e3 ;% "SEI layer density [kg/m^3]"
    i_0_SEI_tune = 1.02;
    i_0_SEI = i_0_SEI_tune*8.8e-4 ;%"SEI Exchange Current [A/m^2]"
    OCP_SEI = 0 ; %SEI OCP [V]
    alpha_SEI = 0.5 ;%"Charge transfer coefficient, SEI"

    kappa_li = 1/(92.8e-9); % Lithium Conductivity [S/m], taken from wikipedia
    M_Li = 6.94e-3; % Molar Mass of Lithium from RSC website [Kg/mol]
    rho_Li = 0.534e3; % Density of Lithiun from RSC website [Kg/m^3]

    k_0_SEI_tune = 1.07;
    k_SEI = 3e-10;
    c_sol_0 = 4541; % Concentration of Solvent in bulk electrolyte [mol/m^3]
    D_SEI =  3.5e-20 ;% "EC diffusion coefficient [m^2/s]"
    epss_SEI = 0.03;
    c_sol_0 = epss_SEI*c_sol_0;

    %% Anode Plating Specific Parameters
    k_li = 3.0e-10; % Lithium plating rate constant
    alpha_pl_neg = 0.3; % Anodic Plating Transfer Coefficient
    alpha_pl_pos = 0.7; % Cathodic Plating Transfer Coefficient
    OCP_li =0; % OCP Plating Potential
    


end