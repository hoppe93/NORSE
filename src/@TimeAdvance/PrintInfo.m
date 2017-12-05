function PrintInfo(o,isARestart,infoStr)
    % Prints the preliminary info about the time advance scheme.
    %
    % Usage:
    %   PrintInfo(isARestart,infoStr)
    %
    % The first argument is a boolean, the second is a string
    % describing the time-advance scheme.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    oN = o.norse;
    
    oN.Print('----------------------------------------------------------------------\n'); 
    if isARestart 
        oN.Print(' Continuing the calculation.\n'); 
    end
    oN.Print(infoStr);
    if oN.show1DTimeEvolution
        oN.Initialize1DTimeEvolutionPlot();
    end   
end
