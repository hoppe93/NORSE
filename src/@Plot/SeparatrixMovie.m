function SeparatrixMovie(o,varargin)
    % Shows a movie overlaying the evolution of the 2D distribution
    % and the runaway-region separatrix. The separatrix is
    % calculated using the runaway-region mode of the NORSE run,
    % unless an additional argument is passed (-1 also gives the
    % mode used in the run). A vector of contour lines (in
    % log10(F)) for the plot of the distribution can also be
    % specified (Default: linspace(-20,1,15)).
    %
    % Usage: 
    %   o.SeparatrixMovie()
    %   o.SeparatrixMovie(runawayRegionMode)
    %   o.SeparatrixMovie(runawayRegionMode,contours)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    useAutoTrim = 0;
    contours    = linspace(-20,1,15);            

    %%%%%%%%%%%%

    oN = o.norse;                        

    %Handle input
    switch nargin
        case 1
            rRMode   = -1; 
        case 2
            rRMode   = varargin{1};                    
        case 3
            rRMode   = varargin{1};
            contours = varargin{2};
        otherwise
            error('Invalid number of input arguments.')
    end

    % Prepare switches for the separatrix calculation
    if rRMode<0
        rRMode = oN.runawayRegionMode; 
    end
    actualRERegionMode    = oN.runawayRegionMode;
    actualRERegion        = oN.runawayRegion;
    oN.runawayRegionMode  = rRMode;                
    oN.runawayRegion.EHat = -1; %Make sure the region is recalculated        

    %Prepare figure
    hF = figure(o.GetFigId('SeparatrixMovie'));
    clf;            
    colormap(MyColorMap(numel(contours)));

    %Loop through the time steps and plot the dist, separatrix and
    %force balance
    for i = 1:size(oN.f,2)                               
        clf(hF);    

        %%% The dist and separatrix %%%
        subplot(3,1,1:2);

        %Prepare and plot f
        f = o.grid.MapBigVectorToGrid(oN.f(:,i));
        f(f<0) = realmin;
        f = log10(f);
        contourf(o.grid.pPara2D,o.grid.pPerp2D,f,contours);
        hold on;

        %Calculate and plot separatrix
        oN.CalculateRunawayFraction(i);
        if any(oN.runawayRegion.pcs < o.grid.pMax)
            pPara = o.grid.xi.*oN.runawayRegion.pcs;
            pPerp = oN.runawayRegion.pcs.*sqrt(1-o.grid.xi.*o.grid.xi);
            hSep = plot(pPara,pPerp,'--r','Linewidth',2,...
                                    'DisplayName','Separatrix');
        end

        %Label the plot
        hL = legend(hSep, 'Location','NorthEast');
        hL.Interpreter = 'latex';
        colorbar
        axis equal
        xlabel('$p_{||}$','Interpreter','latex');
        ylabel('$p_{\perp}$','Interpreter','latex');
        title(sprintf('$\\tau=%.2f$, runaway-region mode: %d',...
                            oN.times(i),oN.runawayRegionMode),...
                                            'Interpreter','latex');

        %Set plot limits
        if useAutoTrim
            %Automatically restrict the plot to show only the
            %interesting parts. Use both f and the last separatrix
            %(nonlinear traj.) to determine suitable boundaries.

            cutOff = contours(1);

            %Clean up the separatrix point arrays
            pPara(isinf(pPara)) = [];
            pPara(isnan(pPara)) = [];
            pPerp(isinf(pPerp)) = [];
            pPerp(isnan(pPerp)) = [];

            %Restrict in positive p_para direction                    
            pCutPos = 1.05*max(o.grid.pPara2D(f>cutOff)); %based on f
            pCutPos = max([pCutPos,max(pPara)]); %based on the separatrix

            %Restrict in negative p_para direction
            pCutNeg = 1.05*min(o.grid.pPara2D(f>cutOff)); 
            pCutNeg = min([pCutNeg,min(pPara)]);

            %Restrict in p_perp direction
            pPerpCut = 1.05*max(o.grid.pPerp2D(f>cutOff)); 
            pPerpCut = max([pPerpCut,max(pPerp)]);

            %Apply the restrictions
            xlim([pCutNeg,pCutPos])
            ylim([0,pPerpCut])
        else
            %Show the entire computational grid
            xlim([-o.grid.pMax,o.grid.pMax])
            ylim([0,o.grid.pMax])
        end



        %%% The force balance at xi=1 %%%
        subplot(3,1,3);

        %Calculate and plot force balance
        [pc,dpdtOnGrid] = oN.FindParallelPCrit(i);
        plot(o.grid.p,dpdtOnGrid,'-b','LineWidth',2);
        hold on;
        plot([0,o.grid.pMax],[0,0],':k','LineWidth',1);
        hPc = plot([pc,pc],get(gca,'YLim'),'--r','LineWidth',1);

        %Label the plot
        xlabel('$p_{||}$','Interpreter','latex');
        ylabel('Sum of forces ($\xi=1$)','Interpreter','latex');
        if useAutoTrim
            xlim([0,pCutPos]);
        else
            xlim([0,o.grid.pMax]);
        end
        hL = legend(hPc,'$p_c$','Location','SouthEast');
        hL.Interpreter = 'latex';

        drawnow
    end

    %Restore the final calculation state
    if oN.runawayRegionMode ~= actualRERegionMode;
        oN.runawayRegionMode = actualRERegionMode;
        oN.runawayRegion = actualRERegion;                
    end
end
