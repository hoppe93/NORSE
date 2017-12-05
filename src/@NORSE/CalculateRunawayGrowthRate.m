function rate = CalculateRunawayGrowthRate(o)
    % Calculates the runaway growth rate n^{-1} dn_r/dt from the
    % runaway fraction at the saved time steps.
    %
    % Usage:
    %   growthRate = CalculateRunawayGrowthRate()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(o.runawayFraction)
        error('Cannot calculate the runaway growth rate before the runaway fraction is known!');
    end            

    rate = [0,diff(o.runawayFraction)]./[1,diff(o.times)];
    %%% OBS! This depends on what time steps are saved! For
    %%% accurate results we need to save a lot of steps!
end
