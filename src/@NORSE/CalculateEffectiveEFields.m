function [EEc,EED] = CalculateEffectiveEFields(o,varargin)
    % Uses the effective temperature based on the energy content of
    % the bulk of the distribution to calculate effective values
    % for E/E_c(T_eff) and E/E_D(T_eff).
    %
    % Usage: 
    %   [EOverEc,EOverED] = CalculateEffectiveEFields()
    %   [EOverEc,EOverED] = CalculateEffectiveEFields(type)
    %
    % type -- Optional string argument, either 'applied' or 'inductive'.
    %         Determines whether to use the applied or consistently
    %         calculated inductive electric field when calculating E/E_c
    %         and E/E_D. Only has an effect when the inductive electric
    %         field is used. 
    %         Default: 
    %           Standard NORSE run  -- 'applied'
    %           Inductive NORSE run -- 'inductive'
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(o.effectiveBulkTemperature)
        error('You need to run NORSE first!');
    end

    e = o.constants.e;            
    m = o.constants.m;
    c2 = (o.constants.c)^2;
    
    %Check the type of electric field to use
    if o.timeAdvanceMode == 3 && ...
                        (nargin >= 2 && ~strcmpi(varargin{1},'applied'))
        %Use the inductive electric field
        EHat = o.EHatInductive;
    else
        %Use the applied electric field
        EHat = o.EHat(o.times); %Convert from time-dep. param. to vector
    end

    %Get the E-field in V/m
    EcOld = o.CalculateEc();
    E = EHat.*EcOld(o.times);

    %Calculate the new normalized E-fields
    TEff = o.effectiveBulkTemperature;
    if all(o.bulkDensity==0) %Use the bulk density unless it is unavailable
        n = o.density; 
    else
        n = o.bulkDensity;
    end
    Ec = o.CalculateEc(TEff,n);                        
    EEc = E./Ec;
    EED = EEc.*TEff*e/(m*c2);
end
