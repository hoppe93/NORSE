function ProcessPhysicalParameters(o)
    % Verifies the provided physical parameters and calculates
    % aggreagate quantities needed in NORSE.
    %
    % Usage:
    %   ProcessPhysicalParameters()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t0 = o.referenceTime;

    %Initiate TimeDependentParameter objects for the physical
    %properties or use predefined ones
    if any([isempty(o.T),isempty(o.n),isempty(o.Z),isempty(o.EHat),isempty(o.B)])
        error('Not all physical parameters are defined.');
    end
    if any(~[isscalar(o.T),isscalar(o.n),isscalar(o.Z),isscalar(o.EHat),isscalar(o.B)])
        error('The physical parameters must be scalar or TimeDependentParameter objects.');
    end            
    %Create TimeDependent Parameter objects from any scalar input
    if ~isa(o.T,'TimeDependentParameter') 
        o.T = TimeDependentParameter(t0,o.T);                        
    end 
    if ~isa(o.n,'TimeDependentParameter') 
        o.n = TimeDependentParameter(t0,o.n);
    end 
    if ~isa(o.Z,'TimeDependentParameter') 
        o.Z = TimeDependentParameter(t0,o.Z);
    end 
    if ~isa(o.EHat,'TimeDependentParameter')
        o.EHat = TimeDependentParameter(t0,o.EHat);
    end  
    if ~isa(o.B,'TimeDependentParameter')
        o.B = TimeDependentParameter(t0,o.B);
    end 
    if ~isa(o.LHat,'TimeDependentParameter') 
        o.LHat = TimeDependentParameter(t0,o.LHat);                        
    end

    %Print some info
    if ~all([o.T.isScalar,o.n.isScalar,o.Z.isScalar,o.EHat.isScalar,...
                                            o.B.isScalar,o.LHat.isScalar])
        o.Print('   Some physical parameters are time-dependent.\n');
    else
        o.Print('   Using constant physical parameters.\n');
    end

    %Calculate some aggregate quantities            
    m = o.constants.m;       
    o.Theta = o.T*( o.constants.e/(m*(o.constants.c)^2) ); %T in eV
    o.kappa = besselk(2,1/o.Theta,1); %exp(1/Theta)*K_2(1/Theta)
    o.fM0 = o.n(t0) / (4*pi*o.Theta(t0)*o.kappa(t0)); %Normalization of f
    lnLambda = o.CalculateLnLambda();
    o.sigma = 2/3*o.B.^2*o.constants.eps0/(m*o.n(t0)*lnLambda(t0)); 
                                %Ratio of collision time to 
                                %synchrotron radiation time scale

    o.nBar = o.n/o.n(t0);
    o.lnLambdaBar = lnLambda/lnLambda(t0);
    o.ThetaBar = o.Theta/o.Theta(t0);
    o.kappaBar = o.kappa/o.kappa(t0);

    o.maxwellianPreFactor = o.nBar/(o.ThetaBar*o.kappaBar);            
end
