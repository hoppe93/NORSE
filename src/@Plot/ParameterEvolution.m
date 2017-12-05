function ParameterEvolution(o)
    % Compares the prescribed temperature and density evolution to
    % the actual ones obtained in NORSE.
    %
    % Clicking on a point in the plots produces a new 1D plot with
    % the parallel distribution at that time.
    %
    % Usage:
    %   ParameterEvolution()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;            
    figure(o.GetFigId('ParameterEvolution'));
    clf;            

    %Density
    subplot(2,1,1)
    plot(oN.times,oN.n(oN.times),'-',oN.times,oN.density,'--',...
                                    'ButtonDownFcn',@o.PlotAtTimeX);
    ylabel('Density','Interpreter','latex')
    hL = legend('Prescribed $n$','Density moment of $F$','Location','best');
    hL.Interpreter = 'latex';

    %Temperature
    subplot(2,1,2)
    plot(oN.times,oN.T(oN.times),'-',oN.times,oN.effectiveBulkTemperature,'--',...
                                        'ButtonDownFcn',@o.PlotAtTimeX)
    hL = legend('Prescribed $T$','Effective $T$ of $F_{\mathrm{bulk}}$',...
                                                'Location','best');
    hL.Interpreter = 'latex';
    ylabel('Bulk temperature','Interpreter','latex')
    xlabel('$\tau$','Interpreter','latex');
end
