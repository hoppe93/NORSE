function SetTimeDependentParameters(o,Tt,Tv,nt,nv,Zt,Zv,Et,Ev,Bt,Bv)
    %Convenient method for setting all time-dependent parameters in
    %one call, if time and data vectors are available. The units
    %are: T: eV, n: m^(-3), E: V/m, B: T.
    %NOTE that the property referenceTime should be set before this
    %function is called, since its value is used here.
    %
    % Usage:
    %   o.SetTimeDependentParameters(TTimes,TVals,nTimes,nVals,...
    %                       ZTimes,ZVals,ETimes,EVals,BTimes,BVals)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    isStepWiseConst = 0;
    isConstOutsideRange = 1;

    o.T = TimeDependentParameter(Tt,Tv,isStepWiseConst,isConstOutsideRange);
    o.n = TimeDependentParameter(nt,nv,isStepWiseConst,isConstOutsideRange);
    o.Z = TimeDependentParameter(Zt,Zv,isStepWiseConst,isConstOutsideRange);
    o.B = TimeDependentParameter(Bt,Bv,isStepWiseConst,isConstOutsideRange);            
    E = TimeDependentParameter(Et,Ev,isStepWiseConst,isConstOutsideRange);

    Ec = o.CalculateEc();
    o.EHat = E/Ec(o.referenceTime);            
end
