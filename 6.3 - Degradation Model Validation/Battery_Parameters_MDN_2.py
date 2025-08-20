import numpy as np

def Battery_Parameters_P2D(param,n_parameters):
    Model_Parameters = np.zeros((n_parameters))
    

    if param == 1:
        # From Marquis et al
        # Cell properties and constants
        current_collector = 0  # Model has no current collectors
        A = 1  # [m^2]
        
        # Constants
        R_const = 8.314  # Ideal gas constant [J/(K*mol)]
        F = 96487  # Faraday Constant [C/mol]
        T_c = 25.15  # Temperature in Celsius
        T = 273 + T_c  # Temperature in Kelvin
        T_amb = T

        # Bruggman constants
        
        brugg = 1.5
        brugg_neg, brugg_sep, brugg_pos = brugg, brugg, brugg

        # Cell dimensions
        L_ncc = 0
        L_neg = 0.0001  # [m]
        L_sep = 2.5e-05  # [m]
        L_pos = 0.0001  # [m]
        L_pcc = 0
        L_tot_active = L_neg + L_sep + L_pos  # [m]
        L_tot = L_ncc + L_tot_active + L_pcc  # [m]

        # Separator/Electrolyte Parameters
        C_e_initial = 1000.0  # Initial electrolyte lithium concentration [mol/m^3]
        D_elec = 5.34e-10  # Diffusion Coefficient for Lithium in electrolyte [m^2/s] 
        D_sep_e_eff = D_elec
        t_li = 0.363  # Lithium Ion Transference Number [unitless]
        eps_sep_elec = 0.4  # Electrolyte volume fraction
        kappa_const = 0.6013 # Function Lookup
        dlnfdlnce = 1

        # Anode Parameters
        #soc_initial_neg = 0.53
        soc_initial_neg =0.8
        soc_initial_neg_100 = 0.8
        soc_initial_neg_80 = 0.7948
        soc_initial_neg_60 = 0.6722
        soc_initial_neg_40 = 0.6181
        eps_neg_elec = 0.3  # Electrolyte volume fraction (porosity)
        eps_neg = 0.6 # Solid volume fraction
        C_max_neg = 24983.2619938437  # Maximum solid lithium concentration [mol/m^3]
        C_neg_initial = soc_initial_neg * C_max_neg
        Rs_neg = 10e-6  # Electrode particle size [m]
        D_neg_s = 3.9e-14  # Solid diffusion coefficient [m^2/s] 
        D_neg_eff_s = D_neg_s
        D_neg_e_eff = D_elec
        sigma_neg = 100  # Electrical conductivity [S/m]
        #K_0_neg = 1.94e-11  # Negative rate constant #Function
        m_ref_neg = 2e-5
        K_0_neg = m_ref_neg/F
        As_neg = 3 * (eps_neg / Rs_neg)  # Surface area to volume ratio [m^-1]
        #Cdl_neg = 100
        Cdl_neg = 0.2

        # Anode SEI Growth Parameters
        L_SEI =  1e-9
        #L_SEI =  5e-11
        kappa_SEI = 5e-6 
        M_SEI = 0.162 
        rho_SEI = 1.69e3 
        #i_0_SEI = 8e-4 
        i_0_SEI = 2e-3 
        #OCP_SEI = 0.4 
        #OCP_SEI = 0.2 
        OCP_SEI = 0
        alpha_SEI = 0.5 
        #alpha_SEI = 0.4
        c_SEI_init = (L_SEI*rho_SEI/M_SEI)*As_neg
        
        R_film_neg = L_SEI/kappa_SEI

        k_SEI = 1.0e-11 # SEI Reaction Rate constant [m/s]
        #k_SEI = 1.0e-10
        c_sol_0 = 4541  # Concentration of Solvent in bulk electrolyte [mol/m^3]
        #D_SEI =  5e-22  # "EC diffusion coefficient [m^2/s]"
        D_SEI =  1e-21  # "EC diffusion coefficient [m^2/s]"
        #OCP_li = 0  # SEI OCP [V]
        OCP_li =0.165 # OCP Plating Potential adjusted
        # Anode Plating Specific Parameters
        k_li = 3.0e-10 # Lithium plating rate constant
        alpha_pl_neg = 0.3 # Anodic Plating Transfer Coefficient
        alpha_pl_pos = 0.7 # Cathodic Plating Transfer Coefficient
        c_li_0 = 0  # Initial Concentration of plated lithium [mol/m^3]
        kappa_li = 1 / (92.8e-9)  # Lithium Conductivity [S/m]
        M_Li = 6.94e-3  # Molar Mass of Lithium [Kg/mol]
        rho_Li = 0.534e3  # Density of Lithium [Kg/m^3]

        # Cathode Parameters
        #soc_initial_pos = 0.17
        soc_initial_pos =0.6
        soc_initial_neg_100 = 0.6
        soc_initial_neg_80 = 0.5877
        soc_initial_neg_60 = 0.6427
        soc_initial_neg_40 = 0.6841
        eps_pos_elec = 0.3  # Electrode Porosity
        eps_pos = 0.5  # Solid volume fraction
        C_max_pos = 51217.9257309275  # Maximum solid lithium concentration [mol/m^3]
        C_pos_initial = soc_initial_pos * C_max_pos
        Rs_pos = 10e-6  # Electrode particle size [m]
        D_pos_s = 1e-13  # Solid diffusion coefficient [m^2/s] 
        D_pos_eff_s = D_pos_s
        D_pos_e_eff = D_elec
        sigma_pos = 10  # Electrical conductivity [S/m]
        #K_0_pos = 2.16e-11  # Cathode rate constant #Function
        m_ref_pos = 6e-7
        K_0_pos = m_ref_pos/F
        As_pos = 3 * (eps_pos / Rs_pos)  # Surface area to volume ratio [m^-1]
        #Cdl_pos = 0
        Cdl_pos = 0.2
        R_film_pos = 0

    

    

    if param == 2:

         # From Mohtat et al (NMC111/C6) Cell
        # Cell properties and constants

    
        A = 1.0  # [m^2]

        # Constants
        R_const = 8.314  # Ideal gas constant [J/(K*mol)]
        F = 96487  # Faraday Constant [C/mol]
        T_c = 25.15  # Temperature in Celsius
        T = 273 + T_c  # Temperature in Kelvin
        T_amb = T

        # Bruggman constants
        brugg = 1.5
        brugg_neg = brugg
        brugg_sep = brugg
        brugg_pos = brugg

        # Cell dimensions
        #L_neg = 6.2e-05  # [m]
        #L_sep = 1.2e-05  # [m]
        #L_pos = 6.7e-05  # [m]

        L_neg = 120e-6 # [m]
        L_sep = 30e-06  # [m]
        L_pos = 160e-6 # [m]
        #L_neg = 175e-6 # [m]
        #L_sep = 2.5e-05  # [m]
        #L_pos = 150e-6 # [m]
        L_tot_active = L_neg + L_sep + L_pos  # [m]
        L_tot = L_tot_active
        dlnfdlnce = 1

        # Separator/Electrolyte Parameters
        C_e_initial = 1000.0  # Initial electrolyte lithium concentration [mol/m^3]
        D_elec = 5.35 * 10 ** (-10)  # Diffusion Coefficient [m^2/s]
        D_sep_e_eff = D_elec
        t_li = 0.38  # Transference Number
        #t_li = 0.2594  # Transference Number
        eps_sep_elec = 0.4  # Electrolyte volume fraction
        kappa_const = 0.6013  # Function Lookup
        dlnfdlnce_sep = 0

        # Anode Parameters
        soc_initial_neg_100 = 0.9
        #soc_initial_neg_100 = 0.9691
        soc_initial_neg_80 = 0.8375
        soc_initial_neg_60 = 0.7349
        soc_initial_neg_40 = 0.6187
        soc_initial_neg = soc_initial_neg_100

        eps_neg_elec = 0.3  # Electrolyte porosity
        eps_neg = 0.61 # Solid volume fraction
        C_max_neg = 28746.0  # Maximum solid lithium concentration [mol/m^3]
        C_neg_initial = soc_initial_neg * C_max_neg
        Rs_neg = 2.5e-06  # Electrode particle size [m]
        D_neg_s = 5.0 * 10 ** (-15) # Solid diffusion coefficient [m^2/s]
        D_neg_eff_s = D_neg_s
        D_neg_e_eff = D_elec
        sigma_neg = 100  # Electrical conductivity [S/m]
        m_ref_neg = 1.061 * 10 ** (-6)
        K_0_neg = m_ref_neg/F
        As_neg = 3 * (eps_neg / Rs_neg)  # Surface area to volume ratio [m^-1]
        Cdl_neg = 0.2
        R_film_neg = 0

        # Cathode Parameters
        soc_initial_pos_100 = 0.27
        #soc_initial_pos_100 = 0.1988
        soc_initial_pos_80 = 0.3294
        soc_initial_pos_60 = 0.4268
        soc_initial_pos_40 = 0.5371
        soc_initial_pos = soc_initial_pos_100

        eps_pos_elec = 0.3  # Electrode Porosity
        eps_pos = 0.445  # Solid volume fraction
        C_max_pos = 35380.0  # Maximum solid lithium concentration [mol/m^3]
        C_pos_initial = soc_initial_pos * C_max_pos
        #Rs_pos = 3.5e-06  # Electrode particle size [m]
        Rs_pos = 2.0e-06  # Electrode particle size [m]
        #D_pos_s = 8 * 10 ** (-15)  # Solid diffusion coefficient [m^2/s]
        D_pos_s = 1.2 * 10 ** (-14)  # Solid diffusion coefficient [m^2/s]
        D_pos_eff_s = D_pos_s
        D_pos_e_eff = D_elec
        sigma_pos = 100  # Electrical conductivity [S/m]
        m_ref_pos = 4.824 * 10 ** (-6)
        K_0_pos = m_ref_pos/F
        As_pos = 3 * (eps_pos / Rs_pos)
        Cdl_pos = 0.2
        R_film_pos = 0

        # Anode SEI Growth Parameters
        L_SEI =  1e-9
        #L_SEI =  5e-11
        kappa_SEI = 5e-6 
        M_SEI = 0.162 
        rho_SEI = 1.69e3 
        #i_0_SEI = 8e-4 
        i_0_SEI = 2e-3 
        #OCP_SEI = 0.4 
        #OCP_SEI = 0.2 
        OCP_SEI = 0
        alpha_SEI = 0.5 
        #alpha_SEI = 0.4
        c_SEI_init = (L_SEI*rho_SEI/M_SEI)*As_neg
        
        R_film_neg = L_SEI/kappa_SEI

        k_SEI = 1.0e-11 # SEI Reaction Rate constant [m/s]
        #k_SEI = 1.0e-10
        c_sol_0 = 4541  # Concentration of Solvent in bulk electrolyte [mol/m^3]
        #D_SEI =  5e-22  # "EC diffusion coefficient [m^2/s]"
        D_SEI =  1e-21  # "EC diffusion coefficient [m^2/s]"
        #OCP_li = 0  # SEI OCP [V]
        OCP_li =0.165 # OCP Plating Potential adjusted
        # Anode Plating Specific Parameters
        k_li = 3.0e-10 # Lithium plating rate constant
        alpha_pl_neg = 0.3 # Anodic Plating Transfer Coefficient
        alpha_pl_pos = 0.7 # Cathodic Plating Transfer Coefficient
        c_li_0 = 0  # Initial Concentration of plated lithium [mol/m^3]
        kappa_li = 1 / (92.8e-9)  # Lithium Conductivity [S/m]
        M_Li = 6.94e-3  # Molar Mass of Lithium [Kg/mol]
        rho_Li = 0.534e3  # Density of Lithium [Kg/m^3]

    # Output Parameters
    # Separator Parameters
    
    Model_Parameters[0] = A  # [m] Cell Area
    Model_Parameters[1] = brugg
    Model_Parameters[2] = L_sep  # Length of separator [m]
    Model_Parameters[3] = C_e_initial  # Initial Electrolyte Lithium Concentration [mol/[m^3]]
    Model_Parameters[4] = D_elec  # Diffusion Coefficient for Lithium in electrolyte [[m^2]/s]
    Model_Parameters[5] = t_li  # Lithium Ion Transference Number [unitless]
    Model_Parameters[6] = eps_sep_elec  # Electrolyte Volume Fraction
    Model_Parameters[28] = dlnfdlnce  # Activity

    # Anode Parameters
    Model_Parameters[7] = L_neg
    Model_Parameters[8] = soc_initial_neg
    Model_Parameters[9] = eps_neg_elec  # Electrolyte Volume Fraction (porosity)
    Model_Parameters[10] = eps_neg  # Solid Volume Fraction
    Model_Parameters[11] = C_max_neg  # Maximum Solid Lithium Concentration [mol/[m^3]]
    Model_Parameters[12] = Rs_neg  # Electrode Particle Size [m]
    Model_Parameters[13] = D_neg_s  # Solid Diffusion Coefficient [m^2/s]
    Model_Parameters[14] = sigma_neg  # Electrical Conductivity [S/m]
    Model_Parameters[15] = K_0_neg  # Negative Rate Constant
    Model_Parameters[26] = Cdl_neg # Negative Electrode Double Layer Capacitance
    Model_Parameters[29] = R_film_neg # SEI Resistance Negative

    # Anode SEI Parameters
    
    Model_Parameters[31] = c_SEI_init
    Model_Parameters[32] = kappa_SEI 
    Model_Parameters[33] = M_SEI 
    Model_Parameters[34] = rho_SEI 
    Model_Parameters[35] = i_0_SEI 
    Model_Parameters[36] = OCP_SEI 
    Model_Parameters[37] = alpha_SEI 
    Model_Parameters[38] = L_SEI 
    Model_Parameters[39] = k_SEI
    Model_Parameters[40] = c_sol_0
    Model_Parameters[41] = D_SEI

    # Anode Plating Parameters
    Model_Parameters[42] = k_li
    Model_Parameters[43] = kappa_li
    Model_Parameters[44] = alpha_pl_neg
    Model_Parameters[45] = alpha_pl_pos
    Model_Parameters[46] = c_li_0
    Model_Parameters[47] = OCP_li
    Model_Parameters[48] = M_Li
    Model_Parameters[49] = rho_Li



    # Cathode Parameters
    Model_Parameters[16] = L_pos
    Model_Parameters[17] = soc_initial_pos
    Model_Parameters[18] = eps_pos_elec  # Electrolyte Volume Fraction
    Model_Parameters[19] = eps_pos  # Solid Volume Fraction
    Model_Parameters[20] = C_max_pos  # Maximum Solid Lithium Concentration [mol/[m^3]]
    Model_Parameters[21] = Rs_pos  # Electrode Particle Size [m]
    Model_Parameters[22] = D_pos_s  # Solid Diffusion Coefficient [m^2/s]
    Model_Parameters[23] = sigma_pos  # Electrical Conductivity [S/m]
    Model_Parameters[24] = K_0_pos  # Negative Rate Constant
    Model_Parameters[25] = kappa_const
    Model_Parameters[27] = Cdl_pos # Positive Electrode Double Layer Capacitance
    Model_Parameters[30] = R_film_pos # SEI Resistance Positive


    return Model_Parameters
        