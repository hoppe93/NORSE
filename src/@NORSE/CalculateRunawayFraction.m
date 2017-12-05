function varargout = CalculateRunawayFraction(o,iteration)
    % Defines the runaway region and uses it to calculate the
    % runaway fraction of the distribution at a given time. Also
    % optionally returns the fraction of the total energy carried
    % by the runaways.
    %
    % Usage:
    %   runawayFraction = CalculateRunawayFraction(iteration)
    %   [runawayFraction,runawayEnergyFraction] ...
    %                        = CalculateRunawayFraction(iteration)
    %
    % iteration is the index in the save arrays (iSave).
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    time = o.times(iteration);
    EEc = o.EHat(time);    

    %Recalculate the lower boundary of the runaway region if EHat
    %has changed or the non-linearly valid definition is used
    %(since the separatrix then depends on the distribution).
    calculateMask = (EEc ~= o.runawayRegion.EHat) || ...
                    (o.runawayRegionMode == 3 && EEc>1);
    if calculateMask
        mask = zeros(o.grid.nP,o.grid.nXi); 
        pcs = Inf*ones(o.grid.nXi,1);

        switch o.runawayRegionMode
            case 0 %Force-balance, isotropic                    
                if EEc>1
                    pc = 1/sqrt(o.EHat(time)-1);
                    idBound = find(o.grid.p>pc,1);                            
                    mask(idBound:end,:) = 1; 
                    pcs = pc*ones(size(o.runawayRegion.pcs));
                end                        
            case 1 %Force-balance, xi-dependent
                if EEc>1
                    xiEEc = o.grid.xi*EEc;
                    id1 = find(xiEEc>1,1);
                    pcs = [Inf*ones(id1-1,1); 1./sqrt(xiEEc(id1:end)-1)];
                    for i = o.grid.xi0Id:numel(pcs)
                        id = find(o.grid.p>pcs(i),1); 
                        mask(id:end,i) = 1;
                    end                                                     
                end
            case 2 %Trajectory-based, similar to Smith et al. [PoP, 12 (2005) 122505]
                if EEc>1
                    xiEEc = 0.5*(o.grid.xi+1)*EEc;
                    id1 = find(xiEEc>1,1);
                    pcs = [Inf*ones(id1-1,1); 1./sqrt(xiEEc(id1:end)-1)];
                    pcs(pcs>o.grid.p(end)) = Inf;
                    for i = id1:numel(pcs)
                        id = find(o.grid.p>pcs(i),1); 
                        mask(id:end,i) = 1;
                    end                                                     
                end              
            case 3 %Trajectory, based on actual distribution
                %Find p_c on the parallel axis (xi=1)
                pc = FindParallelPCrit(o,iteration);
                switch pc
                    case Inf
                        %No information on Pi is available, 
                        %EEc<1 or pc>pMax -- nothing to do
                    case 0
                        %The entire compuational grid is in the runaway region
                        mask(:) = 1;
                        pcs(:) = o.grid.pMax;
                    otherwise                
                        %Integrate 'backwards' in xi to trace out
                        %the entire separatrix. This returns the
                        %separatrix points in reverse order.
                        [~,pS] = o.IntegrateSeparatrix(pc,iteration);

                        pcs(1:numel(pS)) = pS;                        
                        for i = 1:numel(pS)
                            id = find(o.grid.p>pcs(i),1); 
                            mask(id:end,i) = 1;
                        end
                        %Flip the vectors, so that the points agree
                        %with the order of the grid points
                        mask = fliplr(mask); 
                        pcs = flipud(pcs);
                end
            otherwise
                error('Invalid runaway region mode');
        end

        %Save the mask 
        o.runawayRegion.mask = o.grid.MapGridToBigVector(mask);
        o.runawayRegion.pcs = pcs;
        o.runawayRegion.EHat = EEc;
        o.runawayRegion.time = time;
    end

    %Calculate the runaway fraction using the calculated mask
    f = o.f(:,iteration);
    maskedF = (o.runawayRegion.mask.*f);
    runawayFraction = (o.grid.intdpdxi'*maskedF)/(o.grid.intdpdxi'*f);
    switch nargout
        case 0                    
            varargout = {};
        case 1
            varargout = {runawayFraction};
        case 2
            energyInt = (o.grid.intdpdxi.*(o.grid.gammaBig-1))';
            runawayEnergy = energyInt*maskedF;
            totEnergy = energyInt*f;                    
            varargout = {runawayFraction,runawayEnergy/totEnergy};
        otherwise
            error('Invalid number of output arguments');
    end            
end
