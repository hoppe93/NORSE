function Update(o,fls) 
    % Updates the potentials based on the distribution fls in the
    % finite-difference--Legendre-mode basis.
    %
    % Usage:
    %   Update(fls)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	FLS = fls;

	if o.norse.ignoreZerothLegMode
		maxwellian = o.norse.maxwellianPreFactor(0) * ...
			exp( (1-o.grid.gammaBig)/o.norse.Theta(0) );
		maxwellianL = o.norse.MapBigVectorToLegModes(maxwellian);

		FLS(1:o.norse.nP-1) = maxwellianL(1:o.norse.nP-1);
		FLS(end) = maxwellianL(end);
	end

	o.Upsilon0 = o.norse.MapLegModesToBigVector(o.Upsilon0Matrix*FLS);
	o.Upsilon1 = o.norse.MapLegModesToBigVector(o.Upsilon1Matrix*FLS);
	o.Upsilon2 = o.norse.MapLegModesToBigVector(o.Upsilon2Matrix*FLS);             

	o.UpsilonPlus = 4*o.Upsilon2 + o.Upsilon1;
	o.UpsilonMinus = 4*o.Upsilon2 - o.Upsilon1;
	o.Pi = o.norse.MapLegModesToBigVector(o.PiMatrix*FLS);
end
