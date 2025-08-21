function [kappa_const] = kappa_calculator(C_e_reduced,param,kappa_in,eps_elec,brugg)

    %% Calculate Kappa Function
    

    if param == 1

        kappa_const = ( ...
                0.0911 ...
                + 1.9101 .* (C_e_reduced ./ 1000) ...
                - 1.052 .* (C_e_reduced ./ 1000) .^ 2 ...
                + 0.1554 .* (C_e_reduced ./ 1000) .^ 3);

       kappa_const = kappa_const.*(eps_elec.^brugg);
    end

    if param == 2

        kappa_const = 1.0793e-2 + (6.7461e-4).*(C_e_reduced) - (5.2245e-7).*(C_e_reduced.^2) + (1.3605e-10).*(C_e_reduced.^3) - (1.1724e-14).*(C_e_reduced.^4);

        %kappa_const = 1.0793e-2 + ((6.7461e-4) * C_e_reduced) - ((5.2245e-7) * C_e_reduced.^2) + ((1.3605e-10) * C_e_reduced.^3) - ((1.1724e-14) * C_e_reduced.^4);

        kappa_const = kappa_const.*(eps_elec.^brugg);
    end


end