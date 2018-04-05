function Dist2D(o,varargin)
    % Generates a contour plot of the distribution in both (p,xi)
    % and (p_para,p_perp) space. 
    %
    % Clicking at a point in the plots produces a new 1D plot with
    % the distribution at that xi (and all p).
    %
    % Usage:
    %   Dist2D()
    %   Dist2D(idTime)
    %   Dist2D(idTime,contours)
    %   Dist2D(idTime,contours,autoTrim)
    %
    % idTime:   time step to plot (-1 = last step, -2 = next-to- 
    %           last step, etc). Default: -1
    % contours: vector of contours to use. 
    %           Default: linspace(-20,1,25)
    % autoTrim: specifies whether to automatically crop the plots
    %           to remove excess white space. Default: true
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contours = linspace(-20,1,25);
    useAutoTrim = 1;

    %%%%%%%%%%%

    %Handle input
    oN = o.norse;
    iteration = -1;
    switch nargin
        case 1
            %Nothing to do
        case 2
            iteration   = varargin{1};
        case 3
            iteration   = varargin{1};
            contours    = varargin{2};
        case 4
            iteration   = varargin{1};
            contours    = varargin{2};
            useAutoTrim = varargin{3};
        otherwise
            error('Wrong number of input arguments');
    end
    if iteration<0
        iteration = size(oN.f,2)+1+iteration;
    end

    %Prepare the distribution
    cutOff = contours(1);            
    f = o.grid.MapBigVectorToGrid(oN.f(:,iteration));
    f(f<0) = realmin;
    f = log10(f);
    f(f<(cutOff+0.1)) = cutOff;            

    if useAutoTrim %Determine what white space to remove from the figure 
        pCut = 1.05*max(o.grid.p2D(f>cutOff));

        %Do further restriction for the (p_para,p_perp) plot                                
        %--Negative xi
        pCutOnNegSide = 1.05*min(o.grid.pPara2D(f>cutOff));
        %--p_perp
        pPerpCut = 1.05*max(o.grid.pPerp2D(f>cutOff));
    else                
        pCut = o.grid.pMax;                
        pCutOnNegSide = -pCut;
        pPerpCut = pCut;
    end



    %-- (p_para,p_perp)-space -------------------------------------            
    hF = figure(o.GetFigId('Dist2D'));
    clf;
    colormap(o.ColorMap(numel(contours)));

    %Contours of the distribution
    btnFcn = @(obj,event) o.PlotAtXiCyl(iteration,obj,event); 
    contourf(o.grid.pPara2D,o.grid.pPerp2D,f,...
                                    contours,'ButtonDownFcn',btnFcn);
    hold on;

    %Runaway region
    if ~all(oN.runawayRegion.pcs == o.grid.pMax)
        %Do not show it if it is not on the grid

        ppa = o.grid.xi.*oN.runawayRegion.pcs;
        ppe = oN.runawayRegion.pcs.*sqrt(1-o.grid.xi.*o.grid.xi);
        hPc = plot(ppa,ppe,'--r','Linewidth',2,'DisplayName','$p_c(\xi)$');
        hL = legend(hPc,'Location','NorthEast');
        hL.Interpreter = 'latex';
    end

    axis equal
    xlim([pCutOnNegSide,pCut]);
    ylim([0,pPerpCut]);

    %Label the plot            
    xlabel('$p_{||}$','Interpreter','latex')
    ylabel('$p_\perp$','Interpreter','latex')                        

    hC = colorbar;
    hC.Label.Interpreter = 'Latex';
    hC.Label.String = '$\log_{10}(F)$';

    title(sprintf('$\\tau = %.3g$ (relativistic coll. times)',...
                        oN.times(iteration)),'Interpreter','latex');
    %--------------------------------------------------------------

    %-- (p,xi)-space ----------------------------------------------
    hF = figure(o.GetFigId('Dist2D')+1);
    clf;
    colormap(o.ColorMap(numel(contours)));            

    btnFcn = @(obj,event) o.PlotAtXiSph(iteration,obj,event);
    contourf(o.grid.xi2D,o.grid.p2D,f,...
                                contours,'ButtonDownFcn',btnFcn);
    hold on;            
    plot(o.grid.xi,oN.runawayRegion.pcs,'--r','Linewidth',2);

    ylim([0,pCut]);

    xlabel('$\xi$','Interpreter','latex');
    ylabel('$p$','Interpreter','latex');

    hC = colorbar;
    hC.Label.Interpreter = 'Latex';
    hC.Label.String = '$\log_{10}(F)$';

    title(sprintf('$\\tau = %.3g$ (relativistic coll. times)',...
                        oN.times(iteration)),'Interpreter','latex');
    %--------------------------------------------------------------            
end
