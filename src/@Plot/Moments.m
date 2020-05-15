function Moments(o)
    % Plots various moments of the distribution, showing for
    % instance the conservation properties and runaway generation.
    %
    % Usage:
    %   Moments()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;

    figure(o.GetFigId('Moments'));
    nRows = 3;
    nCols = 2;            
    clf;

    %%% Left column
    subplot(nRows,nCols,1)            
    o.FunctionOfTime(oN.density/oN.density(1),'Normalized density',1);
    title('\bf{General moments}','Interpreter','latex')

    subplot(nRows,nCols,3)
    o.FunctionOfTime(oN.energy/oN.energy(1),'Normalized energy',1);

    subplot(nRows,nCols,5)
    o.FunctionOfTime(oN.currentDensity,'Current density (A/m$^2$)',0);

    %%% Right column            
    subplot(nRows,nCols,2)            
    o.FunctionOfTime(log10(max(0,oN.runawayFraction)),'$\log_{10}(n_r/n)$',1);   
    title('\bf{Runaway-related quantities}','Interpreter','latex')

    subplot(nRows,nCols,4)
    %Positive growth rate:
    rate = oN.runawayGrowthRate;
    isNeg = rate<=0;             
    data = rate;
    data(isNeg) = 0;            
    o.FunctionOfTime(log10(data),' ',1);
    hold on;
    %Negative growth rate:
    data = -rate;
    data(~isNeg) = 0;
    o.FunctionOfTime(log10(data),...
                    '$\log_{10}(\mathrm {d}n_r/\mathrm {d}\tau)$',1,1);
    if any(data)
        hL = legend('Positive','Negative','Location','Best');
        hL.Interpreter = 'latex';
    end

    subplot(nRows,nCols,6)
    o.FunctionOfTime(oN.effectiveBulkTemperature,'$T_{\mathrm {eff,bulk}}$ (eV)',0);   
    title('\bf{Other}','Interpreter','latex')
end
