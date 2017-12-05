function CalculateWeightedPotentialMatrices(o,t) 
    % Weighs the matrices used to calculate the potentials from f in the
    % Legendre basis by an exponential factor. This way, the number of
    % nonzero elements in the field-particle term can be reduced
    % drastically and full matrices can be avoided. t is the time.
    %
    % Usage:
    %   CalculateWeightedPotentialMatrices(t)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    oN  = o.norse;
    oG  = o.grid;
    oP  = o.potentials;
    pCO = o.potentialCutOff;
    
    %First, make a weighing factor that we can apply in
    %the L-basis (i.e. each L-mode is multiplied by the exponential)
    lSize   = size(oP.PiMatrix);
    expDiag = exp(1/oN.Theta(t)*(1-oG.gamma(2:end))) * ones(1,oN.nL);    
    expFactor = spdiags([expDiag(:);1],0,lSize(1),lSize(2)); 
    
    %Now, make weighted mapping matrices and remove points that won't
    %affect the result to reduce nnz
    o.weightedU0Mat    = expFactor*oP.Upsilon0Matrix;
    o.weightedU0Mat( abs(o.weightedU0Mat)<pCO )       = 0;
    
    o.weightedU2Mat    = expFactor*oP.Upsilon2Matrix;
    o.weightedU2Mat( abs(o.weightedU2Mat)<pCO )       = 0;
    
    o.weightedUMinMat  = expFactor*o.UpsilonMinMatrix;
    o.weightedUMinMat( abs(o.weightedUMinMat)<pCO )   = 0;
    
    o.weightedUPlusMat = expFactor*o.UpsilonPlusMatrix;
    o.weightedUPlusMat( abs(o.weightedUPlusMat)<pCO ) = 0;
    
    o.weightedPiMat    = expFactor*oP.PiMatrix;
    o.weightedPiMat( abs(o.weightedPiMat)<pCO )       = 0;    
end
