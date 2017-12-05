function LegendreModes(o,varargin)
    % Plots a number of Legendre modes of the distribution.
    %
    % Usage:
    %   PlotLegendreModes() - The final time step
    %   PlotLegendreModes(step) - A specific time step
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nModes = 8;

    %%%%%%%%%%%

    oN = o.norse;
    p = o.grid.p;

    %Handle inputs
    if nargin == 1
        iteration = numel(oN.times);
    elseif nargin==2
        iteration = varargin{1};
    else
        error('Wrong number of input arguments')
    end

    %Get the Legendre modes
    fls = oN.MapBigVectorToLegModes(oN.f(:,iteration));            

    %Map the modes onto a form suitable for plotting (one mode in
    %each column)
    flsMat = zeros(o.grid.nP,oN.nL);            
    flsMat(2:end,:) = reshape(fls(1:end-1),o.grid.nP-1,oN.nL);            
    flsMat(1,1) = fls(end);            
    fls = flsMat;            

    % Determine which modes to plot (make sure modes 0,1 & 2 are
    % always shown)
    nModes = min(oN.nL,nModes);
    ls = unique([0,1,round(linspace(2,oN.nL-1,nModes-2))]);
    nModes = numel(ls);
    nRows = 2;
    nCols = ceil(nModes/2);

    figure(o.GetFigId('LegendreModes'))
    clf;            
    %Plot the selected modes
    for iL = 1:numel(ls)
        subplot(nRows,nCols,iL)
        semilogy(p,fls(:,ls(iL)+1),'-b',...
                 p,-fls(:,ls(iL)+1),':b');                
        xlabel('$p$','Interpreter','latex')
        ylabel('$F_l$','Interpreter','latex')
        title(sprintf('$l = %d$',ls(iL)),'Interpreter','latex')
        if iL == 1
            hL = legend('$F_l$','$-F_l$'); 
            hL.Interpreter = 'latex';
        end
    end
end
