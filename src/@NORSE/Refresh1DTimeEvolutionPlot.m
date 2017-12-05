function Refresh1DTimeEvolutionPlot(o,f)
    % Update the plot showing the parallel distribution during
    % runtime. 
    %
    % Usage:
    %   Refresh1DTimeEvolutionPlot(f)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Prepare the distribution
    F = o.grid.MapBigVectorToGrid(f);
    F = [flipud(F(:,1)); F(:,end)];
    Fpos = 1e-100*ones(size(F));
    Fneg = Fpos;
    Fpos(F>0) = F(F>0);
    Fneg(F<=0) = F(F<=0);

    %Update the line in the plot
    set(o.timeEvoLine{1},'YData',Fpos)
    set(o.timeEvoLine{2},'YData',-Fneg)            

    axis([-o.grid.pMax/2 o.grid.pMax 1e-25 1])
    drawnow
end
