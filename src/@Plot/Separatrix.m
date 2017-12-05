function Separatrix(o,varargin)
    % Overlays the 2D distribution and the various definitions of
    % the runaway-region separatrix in momentum space. By default,
    % the final time step is used, but a desired step can also be
    % provided as an additional argument (-1 gives the final step,
    % -2 the next-to-final step, etc). A lower cut-off (in
    % log10(F)) for the plot of the distribution can also be
    % specified (Default: -20).
    %
    % Usage:
    %   Separatrix()
    %   Separatrix(step)
    %   Separatrix(step,cutOff)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cutOff = -20;
    nContours = 15;
    useAutoTrim = 1;

    %%%%%%%%%%%%

    oN = o.norse;
    warning off MATLAB:ode45:IntegrationTolNotMet

    %Handle input
    if nargin == 1
        iteration = -1;
    elseif nargin == 2
        iteration = varargin{1};
    elseif nargin == 3
        iteration = varargin{1};
        cutOff = varargin{2};
    else
        error('Wrong number of input arguments')
    end            
    if iteration <= 0
        iteration = numel(oN.times)+1+iteration; 
    end

    %Save the actual runawayRegion used, so that we can restore it
    %later
    actualRERegionMode = oN.runawayRegionMode;
    actualRERegion = oN.runawayRegion;            

    %%%%%%%%%%


    %Prepare figure
    figure(o.GetFigId('Separatrix')+1);
    clf;
    colormap(o.ColorMap(nContours));
    colors = o.ColorMap(5);         %for the separatrices
    styles = {':','-.','-','--'};   %     -- || --       


    %%% Plot with all the different runaway regions %%%            
    subplot(3,1,1:2);            

    %Prepare and plot distribution
    f = o.grid.MapBigVectorToGrid(oN.f(:,iteration));
    f(f<0) = realmin;
    f = log10(f);
    contours = linspace(cutOff,max(f(:)),nContours);            
    contourf(o.grid.pPara2D,o.grid.pPerp2D,f,contours);
    hold on;

    %Plot the various different runaway separatrices
    hs = [];
    for iMode = 0:3
        %Calculate the separatrix
        oN.runawayRegionMode = iMode;
        oN.runawayRegion.EHat = 0; %Makes sure that the pc is recalculated
        oN.CalculateRunawayFraction(iteration);

        %Plot it
        pPara = o.grid.xi.*oN.runawayRegion.pcs;
        pPerp = oN.runawayRegion.pcs.*sqrt(1-o.grid.xi.*o.grid.xi);
        pPerp(pPerp>o.grid.pMax) = NaN;
        hs(iMode+1) = plot(pPara,pPerp,'Color',colors(iMode+1,:),...
                           'LineStyle',styles{iMode+1},'Linewidth',2);                
    end

    %Label plot
    hL = legend(hs,{'Isotropic f.b.','$\xi$-dep. f.b.','Linear traj.','Nonlinear traj.'},...
                                    'Location','Best');
    hL.Interpreter = 'latex';
    xlabel('$p_{||}$','Interpreter','latex');
    ylabel('$p_{\perp}$','Interpreter','latex');            
    colorbar
    axis equal

    %Set plot limits
    if useAutoTrim
        %Automatically restrict the plot to show only the
        %interesting parts. Use both f and the last separatrix
        %(nonlinear traj.) to determine suitable boundaries.

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
        xlim([pCutNeg ,pCutPos])
        ylim([0,pPerpCut])
    else
        %Show the entire computational grid
        xlim([-o.grid.pMax,o.grid.pMax])
        ylim([0,o.grid.pMax])
    end

    %Restore the final runaway region obtained in the NORSE run
    oN.runawayRegionMode = actualRERegionMode;
    oN.runawayRegion = actualRERegion;




    %%% Force balance at xi=1 %%%
    subplot(3,1,3);

    %Calculate force balance
    [pc,dpdtOnGrid] = oN.FindParallelPCrit(iteration);

    %Plot it
    plot(o.grid.p,dpdtOnGrid,'-b','LineWidth',2);
    hold on;
    plot([0,o.grid.pMax],[0,0],':k','LineWidth',1);
    plot([pc,pc],get(gca,'YLim'),'--r','LineWidth',1);

    %Label the plot
    xlabel('$p_{||}$','Interpreter','latex')
    ylabel('Sum of forces ($\xi=1$)','Interpreter','latex')
    if useAutoTrim
        xlim([0,pCutPos]);
    else
        xlim([0,o.grid.pMax]);
    end



    warning on MATLAB:ode45:IntegrationTolNotMet            
end
