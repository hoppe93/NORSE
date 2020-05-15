function GetSize(o) 
    % Calculates the size in memory of the AVALANCHESOURCE object
    % by looping through the fields and summing up the variable
    % sizes, and prints it to the console.
    %
    % Usage: 
    %   GetSize()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    props = properties(o);
    totSize = 0;
    for i = 1:length(props)
        currentProperty = getfield(o, char(props(i)));
        s = whos('currentProperty');
        totSize = totSize + s.bytes;
    end
    fprintf(1, 'AvalancheSource: %.1f MB\n', totSize/(1024*1024));
end
