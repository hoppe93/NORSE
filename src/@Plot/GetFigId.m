function id = GetFigId(o,str)
    % Returns an appropriate figure number, determined by the
    % calling function and the figure offset of the NORSE object.
    % The argument is the name of the calling plotting function.
    %
    % Usage:
    %   figId = GetFigId(str)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    id = o.GetFigOffset(); 
    switch str
        case 'Dist1D'
            plotId = 100;                               
        case 'Dist2D'
            plotId = 110;
        case 'DistMovie1D'
            plotId = 105;
        case 'DistMovie2D'
            plotId = 115;
        case 'DistVsTime'
            plotId = 120;
        case 'Dist1DAtXi'
            plotId = 130;
        case 'LegendreModes'
            plotId = 140;

        case 'Moments'
            plotId = 200;
        case 'ParameterEvolution'
            plotId = 210;
        case 'HeatSink'
            plotId = 220;                
        case 'GMRESInfo'
            plotId = 230;
        case 'TimeStep'
            plotId = 240;

        case 'EffectiveEfields'
            plotId = 300;                
        case 'Separatrix'
            plotId = 310;
        case 'SeparatrixMovie'
            plotId = 315;        
        otherwise
            %Unspecified plot
            plotId = 1;
    end     
    id = id + plotId;
end
