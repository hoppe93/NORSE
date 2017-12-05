function CalculateEnergyChangeMagnitude(o,f,t,tOld)
    % Calculates the appropriate magnitude of the heat sink to
    % counteract the Ohmic heating, synchrotron losses, and effects
    % of collisions, as well as carry out temperature or
    % density-change-related energy changes.
    %
    % Usage:
    %   CalculateEnergyChangeMagnitude(f,t,tOld)
    %
    % f and t are the distribution and time at the current time
    % step. tOld is the time at the previous time step.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sMagn = o.magnitudeStructTemplate;               
    oN = o.norse;

    SHInt = o.energyChangeInt*(o.heatSinkOperator*f); %Energy moment of the sink            
    sMagn.normalization = -o.rateNormalization*SHInt; 
        %Conversion factor between the source magnitudes and W/m^3

    %%% Temperature change
    if o.includeTempChanges && (oN.Theta(t) ~= oN.Theta(tOld)) 
        %The temperature has changed, perform the corresponding
        %change in energy content

        %We need to rebuild the heat sink to reflect the
        %new temperature. This doesn't seem to do much in
        %terms of accuracy, but it's not much of a cost
        %either...
        o.Assemble(t); 

        %Calculate the energy change corresponding to the
        %difference in energy content between two Maxwellian at
        %the the relevant temperatures (T(t2) and T(t1)) 
        oMinG = 1-o.grid.gammaBig;           
        preNew = o.preFactor(t);
        preOld = o.preFactor(tOld);
        expNew = exp(oMinG/oN.Theta(t));
        expOld = exp(oMinG/oN.Theta(tOld));
        dW = o.energyChangeInt*(preNew*expNew - preOld*expOld);


        %Determine and save the heat sink magnitude
        sMagn.tempChange = dW/((t-tOld)*SHInt);
        sMagn.sink = sMagn.sink+sMagn.tempChange;                
    end

    %%% Energy change associated with density change
    if o.calculateDensityRelatedChanges
        dn = oN.n(t)-oN.n(tOld);
        if dn     
            %The density has changed, calculate the corresponding
            %heat change.                    
            sMagn.densityChange = o.particleSource.energyContentWBar(t) * ...
                                  dn/oN.fM0/((t-tOld)*SHInt); 
            %Only change the energy content if the user desires it
            if o.includeDensityRelatedChanges
                sMagn.sink = sMagn.sink + sMagn.densityChange; 
            end
        end
    end

    %%% Compensate for E field, collisions and synchrotron losses                            
    EInt = o.energyChangeInt*(oN.EHat(t)*oN.EFieldOperator*f);
    CInt = o.energyChangeInt*(oN.collisionOperator.C*f);
    RRInt = o.energyChangeInt*(oN.sigma(t)*oN.synchrotronOperator*f);                                    
    sMagn.E = EInt/SHInt;
    sMagn.C = CInt/SHInt;
    sMagn.synch = RRInt/SHInt;
    %Only do the compensation if the user desires it
    if o.includeHeatSink
        sMagn.sink = sMagn.sink - sMagn.E - sMagn.C - sMagn.synch;            
    end            

    %Add a correction if the energy content of the current
    %distribution does not agree with the specified value
    if o.enforceStrictHeatConservation 
        %Calculate the energy discrepancy. Make sure we keep the
        %temperature (not energy) constant if the bulk loses
        %particles.                               
        nBarEff = oN.fM0*(o.densityInt*f)/oN.n(tOld); 
        oMinG = 1-o.grid.gammaBig;
        preOld = o.preFactor(tOld);
        expOld = exp(oMinG/oN.Theta(tOld));                    
        tempCorr = o.energyChangeInt*...
                        (nBarEff*preOld*expOld-f/oN.nBar(tOld));

        %Save the correction magnitude
        sMagn.correction = tempCorr/((t-tOld)*SHInt);
        sMagn.sink = sMagn.sink+sMagn.correction;
    end

    %Model a maximum energy outflow rate by limiting the maximum
    %heat sink magnitude.
    if o.restrictHeatSinkRate
        %Do not include the effect of density changes in the
        %restriction, since in that case that energy is carried by
        %the particles themselves
        magnToRestrict = sMagn.sink - sMagn.densityChange;

        excess = (magnToRestrict + o.rateCutOff/SHInt);
        if excess > 0
            sMagn.sink = sMagn.sink-excess;
            sMagn.restriction = -excess;
        end
    end

    o.energyChangeMagnitude = sMagn;
end 
