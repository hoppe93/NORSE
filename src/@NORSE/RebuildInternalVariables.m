function RebuildInternalVariables(o)
    % Re-initializes operators and objects that have been cleared
    % using CleanUp(). Makes it possible to restart the
    % calculation.
    %
    % Usage:
    %   RebuildInternalVariables()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    o.AssembleEFieldOperator();                                                
    o.AssembleSynchrotronOperator();

    o.potentials = Potentials(o);            
    o.potentials.GeneratePotentialMatrices();

    o.particleSource = ParticleSource(o);
    o.heatSink = HeatSink(o);

    if isempty(o.plot)
        o.plot = Plot(o);
    end
end
