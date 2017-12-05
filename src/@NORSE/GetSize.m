function GetSize(o,varargin) 
    % Calculates the size in memory of the NORSE object by looping
    % through the fields and summing up the variable sizes, and
    % prints it to the console. Also calculates the size of the
    % other objects associated with a NORSE calculation. If an
    % argument is passed, it determines whether to show the size of
    % each property of the NORSE class.
    %
    % Usage: 
    %   GetSize()
    %   GetSize(doShowIndividualVariableSizes)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    props = properties(o); 
    totSize = 0; 
    if nargin>=2
        printAll = varargin{1};
    else
        printAll = false;
    end
    for i = 1:length(props) 
        currentProperty = getfield(o, char(props(i))); 
        s = whos('currentProperty'); 
        totSize = totSize + s.bytes;

        if printAll
            %Print the size of each field:
            fprintf(1, '%s: %.1f MB\n', char(props(i)),s.bytes/(1024*1024));
        end
    end
    fprintf(1, 'NORSE: %.1f MB\n', totSize/(1024*1024));

    % Also calculate the size of the referenced objects
    o.grid.GetSize();
    o.timeAdvance.GetSize();
    if ~isempty(o.potentials)                
        o.potentials.GetSize();
        o.collisionOperator.GetSize();
        o.heatSink.GetSize();
        o.particleSource.GetSize();	
    else
        %CleanUp has been called, so there are no objects to query
    end            
end
