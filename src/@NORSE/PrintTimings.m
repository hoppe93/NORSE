function PrintTimings(o)
    % Prints information on the time used to perform various tasks
    % during the NORSE run. May sometimes show incorrect result
    % after a restart.
    %
    % Usage:
    %   PrintTimings()
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~o.silent
        t = o.timing;
        %Print timings --------------------------------------------------------  
        tPerTimeStep = t.timeAdvance/o.nTimeSteps;
        tBuild = t.collOp/o.nTimeSteps;
        fprintf('----------------------------------------------------------------------\n');    
        fprintf(' Total calculation time: %f s\n',t.total);
        fprintf('  - Time to initialize: %f s\n',t.initialization);                
        fprintf('    - Time to build operators for the potentials: %f s\n',t.potentialMats);                
        fprintf('  - Time to advance system in time: %f s\n',t.timeAdvance);                
        fprintf('    - Average time per timestep: %f s\n',tPerTimeStep);                
        fprintf('      - Average time to build the collision operator: %f s (%.1f%%)\n',tBuild,100*tBuild/tPerTimeStep);
        switch o.timeAdvanceMode 
            case 0 
                tInvert = t.matrixInversion/o.nTimeSteps;
                fprintf('      - Average time to invert the matrix: %f s (%.1f%%)\n',...
                                        tInvert,100*tInvert/tPerTimeStep);
            case {1,2}
                tPreCon = t.matrixFactorization/o.nTimeSteps;
                tSolve = t.GMRES/o.nTimeSteps;
                fprintf('      - Average time (per time step) to generate the pre-conditioner: %f s (%.1f%%)\n',...
                                        tPreCon,100*tPreCon/tPerTimeStep);
                fprintf('      - Average time to solve the system iteratively: %f s (%.1f%%)\n',...
                                          tSolve,100*tSolve/tPerTimeStep);
            case 3
                tInvert = t.matrixInversion/o.nTimeSteps;
                tNewton = tInvert/o.nNewtonSteps;
                fprintf('      - Average time per Newton iteration: %f s\n',...
                                        tNewton);
                fprintf('      - Average time to solve the system: %f s (%.1f%%)\n',...
                                        tInvert,100*tInvert/tPerTimeStep);
            case 4
                tInvert = t.matrixInversion/o.nTimeSteps;
                tAnalyse = t.mumpsAnalysis/o.nTimeSteps;
                fprintf('      - Average time to solve the system: %f s (%.1f%%)\n',...
                                        tInvert,100*tInvert/tPerTimeStep);
                fprintf('        - Average time to analyse the system: %f s (%.1f%%)\n',...
                                        tAnalyse,100*tAnalyse/tPerTimeStep);
            otherwise
                warning('Cannot print timings -- unknown time advance mode');
        end
        fprintf('  - Time to process and save data: %f s\n',t.postProcess);
        fprintf('    - Time to calculate the runaway region separatrix: %f s\n',t.separatrix);
        fprintf('*********************************************************************\n\n');
        %----------------------------------------------------------------------
    end 
end        
