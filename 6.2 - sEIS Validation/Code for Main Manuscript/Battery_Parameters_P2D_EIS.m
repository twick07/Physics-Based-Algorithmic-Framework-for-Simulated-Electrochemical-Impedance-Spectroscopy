
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

    if SoC == 100
        soc_initial_neg = soc_initial_neg_100;
    end
    if SoC == 80
        soc_initial_neg = soc_initial_neg_80;
    end
    if SoC == 60
        soc_initial_neg = soc_initial_neg_60;
    end
    if SoC == 40
        soc_initial_neg = soc_initial_neg_40;
    end

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
    R_film_neg = 0; % Film Resistance [Ohm m^2]
    
    
    %% Cathode Parameters
    %soc_initial_pos = 0.17;
    soc_initial_pos = 0.6;
    soc_initial_pos_100 = 0.6;
    soc_initial_pos_80 = 0.5877;
    soc_initial_pos_60 = 0.6427;
    soc_initial_pos_40 = 0.6841;
    if SoC == 100
        soc_initial_pos = soc_initial_pos_100;
    end
    if SoC == 80
        soc_initial_pos = soc_initial_pos_80;
    end
    if SoC == 60
        soc_initial_pos = soc_initial_pos_60;
    end
    if SoC == 40
        soc_initial_pos = soc_initial_pos_40;
    end
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
    R_film_pos =0; % Film Resistance [Ohm m^2]

    %% Anode SEI Specific Parameters
    %L_SEI =  1e-9;%%"Initial SEI Thickness [m]"
    L_SEI = 0;
    kappa_SEI = 5e-6 ;%"SEI ionic conductivity [S/m]"
    M_SEI = 0.162 ;%"SEI Molar Mass [kg/mol]"
    rho_SEI = 1.69e3 ;% "SEI layer density [kg/m^3]"
    i_0_SEI_tune = 1.02;
    i_0_SEI = i_0_SEI_tune*8e-4 ;%"SEI Exchange Current [A/m^2]"
    OCP_SEI = 0 ; %SEI OCP [V]
    alpha_SEI = 0.5 ;%"Charge transfer coefficient, SEI"

    kappa_li = 1/(92.8e-9); % Lithium Conductivity [S/m], taken from wikipedia
    M_Li = 6.94e-3; % Molar Mass of Lithium from RSC website [Kg/mol]
    rho_Li = 0.534e3; % Density of Lithiun from RSC website [Kg/m^3]

    %k_SEI = 1.37e-5 ;%"Rate constant, SEI [m/s]"
    k_0_SEI_tune = 1.07;
    k_SEI = k_0_SEI_tune*1.0e-11; %"Rate constant adjusted, SEI [m/s]"
    c_sol_0 = 4541; % Concentration of Solvent in bulk electrolyte [mol/m^3]
    %D_SEI =  1e-13 ;% "EC diffusion coefficient [m^2/s]"
    D_SEI =  1e-20 ;% "EC diffusion coefficient [m^2/s]"


    %% Anode Plating Specific Parameters
    k_li_tuning = 1.22;
    k_li = k_li_tuning*3.0e-10; % Lithium plating rate constant
    alpha_pl_neg = 0.3; % Anodic Plating Transfer Coefficient
    alpha_pl_pos = 0.7; % Cathodic Plating Transfer Coefficient
    OCP_li =0; % OCP Plating Potential
    %OCP_li =0.165; % OCP Plating Potential
    


end



if param == 2
    
    % From Mohtat et al (NMC111/C6) Cell
    % Cell properties and constants
    
    A = 1.0;  % [m^2]
    
    % Constants
    R_const = 8.314;     % Ideal gas constant [J/(K*mol)]
    F = 96487;           % Faraday Constant [C/mol]
    T_c = 25.15;         % Temperature in Celsius
    T = 273 + T_c;       % Temperature in Kelvin
    T_amb = T;
    
    % Bruggman constants
    brugg = 1.5;
    brugg_neg = brugg;
    brugg_sep = brugg;
    brugg_pos = brugg;
    
    % Cell dimensions
     %{
    L_neg = 6.2e-05;             % [m]
    L_sep = 1.2e-05;             % [m]
    L_pos = 6.7e-05;             % [m]
   
   
    L_neg = 110e-6;%75e-6; % [m]
    L_sep = 24e-6;%17e-6; % [m] 
    L_pos = 129e-6;%40e-6; % [m]
    %}
    L_neg = 120e-6;
    L_sep = 30e-06;
    L_pos = 160e-6;
    %L_neg = 175e-6;  % [m]
    %L_sep = 2.5e-05;  % [m]
    %L_pos = 150e-6;  % [m]
     
    L_tot_active = L_neg + L_sep + L_pos;  % [m]
    L_tot = L_tot_active;
    dlnfdlnce = 1;
    
    %% Separator/Electrolyte Parameters
    C_e_initial = 1000.0;       % [mol/m^3]
    D_elec = 5.35e-10;          % [m^2/s]
    D_sep_e_eff = D_elec;
    t_li = 0.38;                % Transference Number
    %t_li =0.2594;                % Transference Number
    eps_sep_elec = 0.4;         % Electrolyte volume fraction
    %eps_sep_elec = 0.5;         % Electrolyte volume fraction
    kappa_const = 0.6013;       % Function Lookup
    dlnfdlnce_sep = 0;
    
    %% Anode Parameters
    soc_initial_neg_100 = 0.9;
    %soc_initial_neg_100 = 0.9691;
    soc_initial_neg_80 = 0.8268;
    soc_initial_neg_60 = 0.7198;
    soc_initial_neg_40 = 0.5879;

    if SoC == 100
        soc_initial_neg = soc_initial_neg_100;
    end
    if SoC == 80
        soc_initial_neg = soc_initial_neg_80;
    end
    if SoC == 60
        soc_initial_neg = soc_initial_neg_60;
    end
    if SoC == 40
        soc_initial_neg = soc_initial_neg_40;
    end
    
    eps_neg_elec = 0.3;         % Electrolyte porosity
    eps_neg = 0.61;             % Solid volume fraction
    C_max_neg = 28746.0;        % [mol/m^3]
    C_neg_initial = soc_initial_neg * C_max_neg;
    Rs_neg = 2.5e-06;           % [m]
    %Rs_neg = 10e-06;           % [m]
    D_neg_s = 5.0e-15;          % [m^2/s]
    D_neg_eff_s = D_neg_s;
    D_neg_e_eff = D_elec;
    sigma_neg = 100;            % [S/m]
    m_ref_neg = 1.061e-6;
    K_0_neg = m_ref_neg / F;
    As_neg = 3 * (eps_neg / Rs_neg);  % [m^-1]
    Cdl_neg = 0.2;
    R_film_neg = 0;
    
    %% Cathode Parameters
    soc_initial_pos_100 = 0.27; 
    %soc_initial_pos_100 = 0.1988; 

    soc_initial_pos_80 = 0.3311;
    soc_initial_pos_60 = 0.4205;
    soc_initial_pos_40 = 0.5307;

    if SoC == 100
        soc_initial_pos = soc_initial_pos_100;
    end
    if SoC == 80
        soc_initial_pos = soc_initial_pos_80;
    end
    if SoC == 60
        soc_initial_pos = soc_initial_pos_60;
    end
    if SoC == 40
        soc_initial_pos = soc_initial_pos_40;
    end
    
    eps_pos_elec = 0.3;         % Electrode Porosity
    eps_pos = 0.445;            % Solid volume fraction
    C_max_pos = 35380.0;        % [mol/m^3]
    C_pos_initial = soc_initial_pos * C_max_pos;
    %Rs_pos = 3.5e-06;           % [m]
    Rs_pos = 2.0e-06;
    %D_pos_s = 8.0e-15;          % [m^2/s]
    D_pos_s = 1.2e-14; 
    D_pos_eff_s = D_pos_s;
    D_pos_e_eff = D_elec;
    sigma_pos = 100;            % [S/m]
    m_ref_pos = 4.824e-6;
    K_0_pos = m_ref_pos / F;
    As_pos = 3 * (eps_pos / Rs_pos);  % [m^-1]
    Cdl_pos = 0.2;
    R_film_pos = 0;

    %% Anode SEI Specific Parameters
    %L_SEI =  1e-9;%%"Initial SEI Thickness [m]"
    L_SEI = 0;
    kappa_SEI = 5e-6 ;%"SEI ionic conductivity [S/m]"
    M_SEI = 0.162 ;%"SEI Molar Mass [kg/mol]"
    rho_SEI = 1.69e3 ;% "SEI layer density [kg/m^3]"
    i_0_SEI_tune = 1.02;
    i_0_SEI = i_0_SEI_tune*8e-4 ;%"SEI Exchange Current [A/m^2]"
    OCP_SEI = 0 ; %SEI OCP [V]
    alpha_SEI = 0.5 ;%"Charge transfer coefficient, SEI"

    kappa_li = 1/(92.8e-9); % Lithium Conductivity [S/m], taken from wikipedia
    M_Li = 6.94e-3; % Molar Mass of Lithium from RSC website [Kg/mol]
    rho_Li = 0.534e3; % Density of Lithiun from RSC website [Kg/m^3]

    %k_SEI = 1.37e-5 ;%"Rate constant, SEI [m/s]"
    k_0_SEI_tune = 1.07;
    k_SEI = k_0_SEI_tune*1.0e-11; %"Rate constant adjusted, SEI [m/s]"
    c_sol_0 = 4541; % Concentration of Solvent in bulk electrolyte [mol/m^3]
    %D_SEI =  1e-13 ;% "EC diffusion coefficient [m^2/s]"
    D_SEI =  1e-20 ;% "EC diffusion coefficient [m^2/s]"


    %% Anode Plating Specific Parameters
    k_li_tuning = 1.22;
    k_li = k_li_tuning*3.0e-10; % Lithium plating rate constant
    alpha_pl_neg = 0.3; % Anodic Plating Transfer Coefficient
    alpha_pl_pos = 0.7; % Cathodic Plating Transfer Coefficient
    OCP_li =0; % OCP Plating Potential
    %OCP_li =0.165; % OCP Plating Potential
    

   

end