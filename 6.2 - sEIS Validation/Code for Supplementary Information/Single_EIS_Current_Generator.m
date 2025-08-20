function [time_vector, I_app_EIS,idx_last_cycle] = Single_EIS_Current_Generator(Delta_I,Freq,time_step,Cycles_per_freq,last_x_cycles,Points_Per_cycle)


    time_max = (Cycles_per_freq*(1/Freq));

    if (1/Freq)> Points_Per_cycle

        Points_Per_cycle = round(time_max/time_step);


    else

        Points_Per_cycle = Points_Per_cycle*Cycles_per_freq;

    end
    
    time_vector = linspace(0,time_max,Points_Per_cycle);

    I_app_EIS = Delta_I*sin(2*pi*Freq*time_vector);

    %%% Find indexes for last_x cycles

    time_point_end_cycle = (Cycles_per_freq - last_x_cycles)*(1/Freq);

    idxs = find(time_vector>time_point_end_cycle);

    idx_last_cycle = idxs(1);




end