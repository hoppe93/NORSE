classdef AvalancheSource < matlab.mixin.Copyable
    % AVALANCHECH -- Implementation of the Chiu-Harvey avalanche
    %                source, published in
    %                  [Chiu et al, Nucl. Fusion 38 (1998) 1711]
    %                  [Harvey et al, Phys. Plasmas 7 (2000) 4590]
    %                an improved, momentum- and particle conserving
    %                source was also given in
    %                  [Embréus et al, J. Plasma Phys. 84 (2018) 905840102]
    %                and is implemented as an option here.
    %
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialize:
    %   AvalancheSource() -- empty object
    %
    %   In NORSE:
    %       AvalancheCH(oNORSE)
    %       InitializeParts()
    %		SynchSettings()
    %
    %   The avalanche source at time t is built via a call to
    %       Assemble(t).
    %   It is then accessible through the class property avaS.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        avaS
    end
    
    methods (Abstract)
        Assemble(o,t,fOld)
		% Overloads the superclass method. Assembles the precalculated
		% avalanche source parts into one source matrix.
        Initialize(o);
        % Initializes the avalanche source object.
    end

	methods
		function SynchSettings(o)
			% Make sure settings are consistent with those in NORSE.
			% Should be called before the calculation starts (or at restart).

			%Nothing to do here in general
		end
	end

	%%% Internal properties and methods
	properties
		%properties
		norse
		grid
	end

	methods
		function o = AvalancheSource(varargin)
			% Constructor
			%
			% Usage:
			%   AvalancheSource() -- empty object
			%   AvalancheSource(oNORSE) -- call to use in NORSE

			switch nargin
				case 0
					% Nothing to do
				case 1
					if isa(varargin{1},'NORSE')
						o.norse = varargin{1};
						o.grid = o.norse.grid;
					else
						error('The first argument must be a NORSE object.');
					end
				otherwise
					error('Invalid number of input arguments');
			end
		end

		GetSize(o)
		% Calculates the size in memory of the AVALANCHESOURCE object
		% by looping through the fields and summing up the variable
		% sizes, and prints it to the console.
	end
    
end
