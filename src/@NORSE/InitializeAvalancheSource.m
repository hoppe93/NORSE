function InitializeAvalancheSource(o)
    % Initializes the AvalancheSource to use (if any)
    %
    % Usage:
    %   InitializeAvalancheSource()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch o.avalancheSourceMode
        case 0 %No avalanche source term
            o.avalancheSource = ASChiuHarvey.empty;
            str = 'not used';
        case 1 %Chiu-Harvey
            warning('The avalanche source has not been tested and should not be used.');
            o.avalancheSource = ASChiuHarvey(o);
            str = 'Chiu-Harvey';
        case 2 %Chiu-Harvey-Embreus
            warning('The avalanche source has not been tested and should not be used.');
            o.avalancheSource = ASChiuHarvey(o);
            o.avalancheSource.includeEmbreusTerms = 1;
            str = 'Chiu-Harvey with Embreus terms';
            error('The Chiu-Harvey-Embreus avalanche source has not been implemented yet.');
        otherwise
            error('Invalid avalanche source term specified.')
    end
    
    o.Print('   The avalanche source is %s.\n',str);
end
