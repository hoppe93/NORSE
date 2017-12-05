function bools = HasChanged(o,ts)
    % Takes an array of times ts and determines if the parameter has
    % changed inbetween its entries. Returns a logical array the same size
    % as ts. The first entry of the array will always be true.
    %
    % Usage: 
    %   bools = HasChanged(ts)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if o.isScalar
        bools = false(size(ts));
        bools(1) = true;
    else
        vals  = o.Evaluate(ts);
        bools = [1,logical(diff(vals))];
    end
end
