function Initialize(o)
    % Calculates various quantities used in the calculation
    % of the avalanche source term.
    %
    % Usage:
    %   Initialize()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    oN = o.norse;
    o.Ec0 = oN.CalculateEc(oN.T(0), oN.n(0));
    
    o.BuildSourceMatrix();
    
    % Set the avalanche source
    o.S = o.intOperator;

end