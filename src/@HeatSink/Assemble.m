function Assemble(o,varargin)
    % Builds the operator describing the heat sink, as well as
    % quantities needed to determine the heat sink magnitude based 
    % on the distribution and forces or parameter changes.
    %
    % Usage: 
    %   Assemble()
    %   Assemble(t)
    %
    % The optional argument is the time at which to evaluate the
    % operator in case of time-dependent parameters. If
    % unspecified, the reference time is used.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;

    if nargin >= 2
        t = varargin{1};
    else
        t = oN.referenceTime; 
    end

    %%% Build the heat sink %%%
    sz = o.grid.matSize;
    p = o.grid.pBig;                
    switch oN.heatSinkFunction
        case 0 %Relativistic Maxwellian
            width = 1; %Small: wide sink, large: narrow sink
            SH = oN.maxwellianPreFactor(t) ...
                    * exp(-width*(o.grid.gammaBig-1)/oN.Theta(t));
            dSHdp = -width/oN.Theta(t)*p./o.grid.gammaBig.*SH;
            SHOverPAt0 = 0; %Fudge the point at p=0 to avoid div by 0                        
        case 1 %p-weighted relativistic Maxwellian - leads to a small triangular peak at p=0
            width = 1;
            SH = p.*oN.maxwellianPreFactor(t) ... 
                    .* exp(-width*(o.grid.gammaBig-1)/oN.Theta(t));
            dSHdp = exp(-width*(o.grid.gammaBig-1)/oN.Theta(t)) ...
                    - width/oN.Theta(t)*p./o.grid.gammaBig.*SH;
            SHOverPAt0 = 1;
        case 2 %p-weighted relativistic Maxwellian - leads to a small triangular peak at p=0
            width = 0.2;
            SH = p.^2.*oN.maxwellianPreFactor(t) ... 
                    .* exp(-width*(o.grid.gammaBig-1)/oN.Theta(t));
            dSHdp = 2.*p.*exp(-width*(o.grid.gammaBig-1)/oN.Theta(t)) ...
                    - width/oN.Theta(t)*p./o.grid.gammaBig.*SH;
            SHOverPAt0 = 0;
        case 3
            v0 = 0.004;
            SH = p.^2 ./ (v0.^2 + (p ./ o.grid.gammaBig).^2);
            dSHdp = 2 * p .* (v0.^2 + (p./o.grid.gammaBig).^4) ...
                    ./ (v0^2 + (p./o.grid.gammaBig).^2).^2;
            SHOverPAt0 = 0;
        otherwise
            error('Invalid choice of heat sink p dependence')
    end

    %Build the sink operator
    fTerm = spdiags(2./p.*SH + dSHdp,0,sz,sz);
    fTerm(end,end) = 2*SHOverPAt0 + dSHdp(end); 
    dfdpTerm = spdiags(SH,0,sz,sz)*o.grid.ddpMat;
    o.heatSinkOperator = fTerm + dfdpTerm;                

    %Precalculate the integrals we need to determine the energy
    %change due to the various terms, and the density in the
    %specified region. These will later multiply f (or some more
    %complicated quantity).
    if oN.heatSinkCutOff
        %Include only the region close to the bulk in the
        %integration
        sourceCalcCutOff = oN.heatSinkCutOff*sqrt(2*oN.Theta(t));% Units of thermal momenta
        intMask = zeros(o.grid.nP,o.grid.nXi); 
        intMask(o.grid.p<=sourceCalcCutOff,:) = 1;
        intMask = o.grid.MapGridToBigVector(intMask); 
        o.energyChangeInt = (intMask.*o.grid.intdpdxi.*(o.grid.gammaBig-1))';
        o.densityInt = (intMask.*o.grid.intdpdxi)';
    else
        %Include the entire grid                
        o.energyChangeInt = (o.grid.intdpdxi.*(o.grid.gammaBig-1))';
        o.densityInt = o.grid.intdpdxi';
    end            
end
