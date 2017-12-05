function status = InterruptIntegration(p,flag,pMax)
    % Aborts the calculation of the separatrix if the calculation
    % is close to reaching the end of the p grid. If status=1, the
    % ode solver aborts the calculation. See the documentation of
    % odeset for a description of the flag property.
    %
    % Usage:
    %   status = InterruptIntegration(p,flag,pMax)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     status = 0;
     if isempty(flag) && p(end)>0.8*pMax
         status = 1;
     end
end
