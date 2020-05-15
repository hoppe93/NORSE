classdef ASChiuHarvey < AvalancheSource
    % ASCHIUHARVEYEMBREUS -- Implementation of the Chiu-Harvey avalanche
    %                        source, published in
    %                          [Chiu et al, Nucl. Fusion 38 (1998) 1711]
    %                          [Harvey et al, Phys. Plasmas 7 (2000) 4590]
    %                        An improved, momentum- and particle conserving
    %                        source was also given in
    %                          [Embréus et al, J. Plasma Phys. 84 (2018) 905840102]
    %                        and is implemented as an option here.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Usage:
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialization:
    %       AvalancheSource() -- empty object
    %
    %   In NORSE:
    %       AvalancheSource(oNORSE)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Interface methods and output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        Assemble(o,t,fOld)
        % Overloads the superclass method. Assembles the precalculated
        % avalanche source parts into one source matrix.
        BuildSourceMatrix(o,pc);
        % Builds the avalanche source matrix that is used to 
        Initialize(o);
        % Initializes the avalanche source object.
    end

    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Settings %%%
    properties
        % properties
        % 1 = Include Embreus' particle- and momentum-conserving terms in
        %     the Chiu-Harvey operator
        % 0 = Not 1
        includeEmbreusTerms = 0
    end
    
    %%% Internal variables %%%
    properties
        Ec0
        % Critical electric at t = 0, as defined in
        % Connor & Hastie [NF, 15 (1975) 415].
        lnLambda0
        % Coulomb logarithm at t = 0
        intOperator
        % The full avalanche operator, evaluated in all points
        % of the NORSE grid.
        
        % Source term (computed by this class)
        S
    end

    methods
        function o = ASChiuHarvey(varargin)
            % Constructor.
            %
            % Usage:
            %   ASChiuHarveyEmbreus() -- empty object
            %   ASChiuHarveyEmbreus(oNORSE) -- call to use in NORSE
            
            o@AvalancheSource(varargin{:});
        end
    end

end
