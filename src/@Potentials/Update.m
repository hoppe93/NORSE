function Update(o,fls) 
    % Updates the potentials based on the distribution fls in the
    % finite-difference--Legendre-mode basis.
    %
    % Usage:
    %   Update(fls)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    o.Upsilon0 = o.norse.MapLegModesToBigVector(o.Upsilon0Matrix*fls);
    o.Upsilon1 = o.norse.MapLegModesToBigVector(o.Upsilon1Matrix*fls);
    o.Upsilon2 = o.norse.MapLegModesToBigVector(o.Upsilon2Matrix*fls);             

    o.UpsilonPlus = 4*o.Upsilon2 + o.Upsilon1;
    o.UpsilonMinus = 4*o.Upsilon2 - o.Upsilon1;
    o.Pi = o.norse.MapLegModesToBigVector(o.PiMatrix*fls);
end
