function TEff = CalculateEffectiveTemperature(o,energy,density)
    % Determines an effective temperature by finding the Maxwellian
    % with the same energy content as the distribution. Naturally,
    % the effective temperature is only a good approximation if the
    % distribution is close to Maxwellian. The input energy and
    % density should be in units internal to NORSE.
    %
    % Usage:
    %   TEff = CalculateEffectiveTemperature(energy,density)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~energy || ~density
        TEff = 0;
        return
    end
    constant = energy/density + 1;            
    optFunc = @(Theta) besselk(3,1./Theta,1)./besselk(2,1./Theta,1)...
                       -Theta - constant;
    try
        TEff = fzero(optFunc,[1e-8,1000]); %Theta is very likely to stay within these boundaries
    catch
        %If the optimization fails to converge 
        TEff = 0;
    end
    %Convert from Theta to eV
    TEff = TEff * o.constants.m*(o.constants.c)^2/o.constants.e;
end
