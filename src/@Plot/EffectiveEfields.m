function EffectiveEfields(o)
    % Plots the effective E/E_c and E/E_D values during the NORSE
    % run, based on the effective temperature of the distribution.
    % Also shows the effective temperature for comparison.
    %
    % Clicking on a point in the plots produces a new 1D plot with
    % the parallel distribution at that time.
    %
    % Usage:
    %   EffectiveEfields()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;                        
    if oN.timeAdvanceMode == 3
        isInductive = true;
    else
        isInductive = false;
    end
    
    [EEc,EED] = oN.CalculateEffectiveEFields('applied');
    if isInductive
        [EEcI,EEDI] = oN.CalculateEffectiveEFields('inductive');
    end

    if any(EEc)  
        figure(o.GetFigId('EffectiveEfields'));
        clf;

        %E/E_c
        subplot(3,1,1);
        plot(oN.times,EEc,'ButtonDownFcn',@o.PlotAtTimeX);
        if isInductive
           hold on;
           plot(oN.times,EEcI,'--','ButtonDownFcn',@o.PlotAtTimeX);
        end
        xlabel('$\tau$','Interpreter','latex');
        ylabel('$E/E_{\mathrm{c}}$','Interpreter','latex');
        if isInductive 
            legend('Applied','Inductive','Location','Best');
        end

        %E/E_D
        subplot(3,1,2);
        plot(oN.times,EED,'ButtonDownFcn',@o.PlotAtTimeX);
        hold on;
        if isInductive           
           plot(oN.times,EEDI,'--','ButtonDownFcn',@o.PlotAtTimeX);
        end
        xL = get(gca,'XLim');
        plot(xL,0.215*[1,1],'--k');
        if any(EED>0.1)
            ylim([0,max([1.1*max(EED),0.25])]) 
        else
            ylim([0.9*min(EED),1.1*max(EED)])
        end
        xlabel('$\tau$','Interpreter','latex');
        ylabel('$E/E_{\mathrm{D}}$','Interpreter','latex');

        %Effective T
        subplot(3,1,3);
        plot(oN.times,oN.effectiveTemperature,'ButtonDownFcn',@o.PlotAtTimeX);
        xlabel('$\tau$','Interpreter','latex');
        ylabel('$T_{\mathrm{eff}}$ (eV)','Interpreter','latex');
    else
        fprintf('No electric field was applied.\n');
    end
end
