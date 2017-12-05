classdef NORSE < matlab.mixin.Copyable
    % NORSE -- Main class for performing a NORSE calculation.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   Use the provided parameters and perform the calculation
    %   automatically (using default settings):
    %       NORSE(nP,nXi,nL,pMax,dt,tMax,T,n,Z,EHat) 
    %       NORSE(nP,nXi,nL,pMax,dt,tMax,T,n,Z,EHat,B)
    %
    %
    %   Initialize an empty object but do no calculation:
    %       NORSE() -- Settings can be specified individually, or using
    %                    SetParameters(nP,nXi,nL,pMax,dt,tMax,T,n,Z,EHat),
    %                    SetParameters(nP,nXi,nL,pMax,dt,tMax,T,n,Z,EHat,B) 
    %                       or
    %                    SetTimeDependentParameters(TTimes,TVals,
    %                                           nTimes,nVals,ZTimes,ZVals,
    %                                           ETimes,EVals,BTimes,BVals)
    %   Then perform the calculation using 
    %       PerformCalculation() 
    %   or do the steps therein explicitly for increased control.
    %
    %   Continue a previously completed calculation:
    %       ContinueCalculation()
    %   tMax must be changed before the call to reflect the new _total_
    %   calculation time (see also the description of the
    %   ContinueCalculation and Rewind methods).
    %
    %   The result of the calculation is available in the NORSE object.
    %   Some useful fields for examining the output are mentioned below.
    %
    %   The methods CleanUp() and AggressiveCleanUp() are provided to
    %   remove internal large variables before saving. When using the
    %   former, the calculation can be continued after loading, but only
    %   after RebuildInternalVariables() has been called. 
    %
    %   The parameters lnLambda, the critical electric field E_c, the
    %   relativistic collision time, and the normalized inductance LHat can
    %   be conveniently calculated using methods of this class. A method
    %   for converting times and momenta in 'thermal' units into
    %   'relativistic' units is also provided.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %   Settings:
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % %%%  Physical parameters %%%
    %
    %   The following parameters can be specified as either
    %   scalar values or predefined TimeDependentParameter objects.
    %
    %   T    -- Bulk plasma temperature (eV). Used to define an initial
    %           state and a normalization.
    %   n    -- Electron number density (particles per m^3).
    %   EHat -- E/E_c, electric field normalized to the critical value of
    %           Connor & Hastie [NF, 15 (1975) 415].
    %   Z    -- Effective charge of the plasma.
    %   B    -- Magnetic field (T).
    %   LHat -- Inductance in units internal to NORSE (see the class method
    %           NormalizeInductance). Default: 0
    %  
    %   referenceTime -- The parameter values at this time are used for the 
    %                    normalization of the distribution and shape of the
    %                    heat sink. Only matters when using time-dependent
    %                    physical parameters. Default: 0
    %
    %
    % %%%  Numerical parameters %%%
    %
    %   nP   -- Number of points in the p grid.
    %   nXi  -- Number of points in the xi grid (xi = cos(theta)).    
    %   nL   -- Number of Legendre modes to use in the calculation of
    %           the potentials of the non-linear collision operator.
    %   pMax -- Maximum value of the p grid (p = gamma*v/c).
    %   dt   -- Time step (in units of collision times for relativistic
    %           electrons). If timeAdvanceMode = 2 is used, this is the
    %           inital time step.
    %   tMax -- Total calculation time (in the same units as dt).
    %
    %   The class method GetRelativisticNormalization can be used to
    %   convert dt and tMax from thermal collision times to relativistic
    %   collision times.
    %
    %
    % %%%  Time advance %%%
    %
    %   timeAdvanceMode -- Scheme for time advance (default: 1)
    %                         0: Constant time step, direct solution of
    %                            linear system (using \).
    %                         1: Constant time step, indirect solution 
    %                            (GMRES with preconditioner from LU).
    %                         2: Adaptive time step, based on 1, step
    %                            length modified according to number of
    %                            iterations used by GMRES.
    %                         3: Simultaneous consistent solution of the
    %                            distribution and inductive electric-field 
    %                            evolution, using Newtons method in each 
    %                            time step.
    %   GMRESTolerance  -- Convergence tolerance used in calls to GMRES.
    %                      The smaller, the more accurate the solution.
    %                      Default: 1e-12
    %   nStepsBetweenFactorizations -- Steps taken between recalculation of
    %                                  the preconditioner in schemes 1 & 2.
    %                                  Default: 2
    %   dtIncreaseLimit             -- Maximum increase in step length 
    %                                  compared to the initial time step
    %                                  (scheme 2 only). A value of 0 means
    %                                  no limit is enforced. Default: 0
    %   nNewtonSteps    -- Number of Newton iterations to take in each time
    %                      step in scheme 3. A higher number should lead to
    %                      better converged results.
    %
    %
    % %%% Collision operator %%%
    %   collisionOperatorMode -- The electron-electron collision operator  
    %                            to use (default:0)
    %                               0: Nonlinear relativistic operator
    %                               1: Linear relativistic operator
    %                               2: Test-particle term of the linear
    %                                  relativistic operator
    %   potentialCutOff   -- Value at which to cut off the exponential-
    %                        weighted potentials in the field particle term
    %                        of the linear collision operator. This avoids
    %                        C becoming a full matrix, and significantly 
    %                        reduces the computational expense. 
    %                        Default: 1e-30
    %   specialFunctionEvaluation -- Method for calculating the special 
    %                                functions of the test-particle term of
    %                                the linear collision operator.
    %                                   0: Accurate  (using integral)
    %                                   1: Efficient (using Simpson's rule)
    %                                Default: 1
    %
    % %%%  Program flow and behavior %%%
    %
    %   initialDistribution   -- Initial condition for the calculation
    %                             0: A Maxwellian (with parameters at t=0).
    %                             1: A weighted Maxwellian, y^2*fM where 
    %                                y=gamma*v/v_th=p/sqrt(2*Theta).
    %                             2: A shifted Maxwellian. Specified T and
    %                                n are interpreted as being in the
    %                                Maxwellian's rest frame
    %                             3: Two Maxwellians, shifted in different
    %                                directions. Each Maxwellian has T and 
    %                                n in its rest frame.
    %                            Default: 0
    %   maxwellianShift       -- Shift (in units of y) in modes 2 & 3. 
    %                            Default: 3   
    %   runawayRegionMode     -- Definition of runaway region
    %                             0: Force-balance, isotropic
    %                                ( p_c=1/sqrt(EHat-1) for all xi ).
    %                             1: Force-balance, xi-dependent
    %                             2: Trajectory-based (see Smith et al.
    %                                [PoP, 12 (2005) 122505]).
    %                             3: Trajectory, based on NORSE
    %                                distribution, includes non-linear and
    %                                synchrotron effects. Slower than the
    %                                others.
    %                             Default: 2
    %   abortWhenOnlyRunaways -- Aborts the calculation if the distribution
    %                            enters the slide-away regime (all
    %                            particles are in runaway region). Requires
    %                            runawayRegionMode = 3. Default: 0 
    %   GeneralAbortFunction  -- Handle to a function which can be used to
    %                            check for a user-supplied abort criterion.
    %                            The function, which is called in each time
    %                            step, must take the five arguments
    %                            o,t,isSaveStep,iSave and iteration, and
    %                            return true (abort) or false (do not
    %                            abort). Here
    %                               o is the NORSE object
    %                               t is the time (in rel. coll. times)
    %                               isSaveStep determines if data has been
    %                                   saved in this time step
    %                               iSave is the id into the save arrays
    %                                   where data was last saved
    %                               iteration is the time step id
    %                            Default: empty
    %
    % %%%  Grid  %%%
    %
    %   pGridMode  -- Specifies the mapping for the p grid (default: 1)
    %                   0: Uniform
    %                   1: Approximately quadratically increasing spacing 
    %                   2: Approximately cubically increasing spacing 
    %                   3: Approximately quartic increasing spacing
    %                   4: A grid with a tanh step in spacing, giving a 
    %                       dense grid at low p and a sparse grid at high p
    %                       (this option also has some parameters
    %                       associated with it, see Grid.m)
    %   xiGridMode -- Specifies the mapping for the xi grid (default: 1)
    %                   0: Uniform
    %                   1: A sigmoid, denser close to xi=1 
    %                   2: Uniform in the polar angle (arccos(\xi))
    %
    %   See Grid.m for more details.
    %
    % %%%  Sources and sinks %%%
    %
    %   includeHeatSink  -- Determines whether to include a heat sink which
    %                       keeps the bulk energy content constant
    %                       (counteracts heat changes due to the electric
    %                       field, synchrotron losses and collisions).
    %                       Density or temperature changes (time-dependent
    %                       parameters) do not make sense without it.
    %                       Default: 0 (but automatically enabled for
    %                                   time-dependent T or n)
    %   heatSinkCutOff   -- Defines region in which the sink magnitude is
    %                       calculated. Given in units of y=gamma*v/v_th
    %                       (thermal momenta). A value of 0 includes the
    %                       entire momentum grid. Default: 5
    %   heatSinkFunction -- Form of the heat sink (in p; it is isotropic in 
    %                       xi). Default: 0
    %                         0: Maxwellian, fM
    %                         1: p*fM
    %   enforceStrictHeatConservation      -- Tries to maintain the 
    %                                         temperature by making sure
    %                                         that any discrepancy in the
    %                                         heat content is removed,
    %                                         including numerical heating.
    %                                         Default: 0
    %   maxHeatSinkRate                    -- Model a maximum energy outflow
    %                                         rate by limiting the maximum
    %                                         heat sink magnitude.
    %                                         Overrides the strict
    %                                         conservation above. The unit
    %                                         is W/m^3. A value of 0 means
    %                                         no restriction. Default: 0
    %   keepTempConstantWhenDensityChanges -- Determines whether to remove
    %                                         or add energy corresponding
    %                                         to the change in the number
    %                                         of particles. Essentially
    %                                         keeps the bulk temperature,
    %                                         rather than the energy
    %                                         content, constant, which
    %                                         usually is the desired
    %                                         behavior. Default: 1 
    %
    %   dampeningStrength -- Strength of the artificial damping at p=pMax.
    %                        Nonzero values can help reduce ringing
    %                        originating from the grid boundary. Default: 0
    %   dampeningWidth    -- Width in per cent of pMax of the artifial
    %                        dampening. Determines how large part of the
    %                        computational region close to pMax which is
    %                        affected by the damping. Default: 1.5 per cent
    %
    % %%%  Saving and printing %%%
    %
    %   nSaveSteps          -- Number of steps to retain after the
    %                          calculation. Affects memory requirements.
    %                          May be adjusted during runtime to be
    %                          consistent with the final output. A value of
    %                          0 retains all steps (if timeAdvanceMode=0 or
    %                          1). Default: 10
    %   timeStepSaveMode    -- Determines how the save steps are chosen
    %                          from all the computational steps
    %                           0: linearly 
    %                           1: logarithmically increasing separation 
    %                              (useful for cases with an initial 
    %                              transient and subsequent slower 
    %                              relaxation)
    %                          Default: 0
    %   savePotentials      -- Determines whether the relativistic
    %                          potentials are saved. Default = 0
    %   show1DTimeEvolution -- Determines whether to show a plot of the
    %                          distribution during calculation. Default = 0
    %   silent              -- Determines whether to print information to
    %                          the console. Default = 0 (verbose)
    %   figOffset           -- Makes it possible to distinguish figures 
    %                          from different runs of NORSE. Default = 1000
    %
    %
    %
    %
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %   Output & results:
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   plot      -- Object with a number of pre-defined plots for
    %                examining the distribution and its moments. See Plot.m
    %                for details.
    %   grid      -- Object containing the finite-difference grid and
    %                related quantities. See Grid.m for details. 
    %   constants -- Struct containing the physical constants e (elemental
    %                charge), m (electron mass), c (speed of light) and
    %                eps_0 (vacuum permittivity).
    %
    % %%% Times and time steps %%%
    %   times      -- Times at the output steps (in relativistic collision
    %                 times). Corresponds to entries in f and the various
    %                 moment vectors.
    %   nSaveSteps -- Number of output steps.
    %   allTimes   -- Times at all of the calculation steps.
    %   dtsUsed    -- Time step used at all calculation steps.
    %   nTimeSteps -- Total number of calculation steps.    
    %
    % %%% Distribution and potentials %%%
    %
    %   --These all have the dimensions (matrixSize x nSaveSteps) --
    %   f          -- Electron distribution
    %   Upsilon0   -- Relativistic potential
    %   Upsilon1   --       -- || --
    %   Upsilon2   --       -- || --
    %   Pi         --       -- || --
    %
    %   If savePotentials = false, only f is saved.
    %
    % %%% Moments of the distribution %%% 
    %
    %   -- These all have the dimensions (1 x nSaveSteps) --
    %   bulkDensity              -- Number density moment of the non-runaway 
    %                               part of f (m^{-3}) 
    %   currentDensity           -- (A/m^2)
    %   density                  -- Number density moment of f (m^{-3}) 
    %   effectiveTemperature     -- Temperature of a Maxwellian with the 
    %                               same energy content as f (eV)
    %   effectiveBulkTemperature -- Temperature of a Maxwellian with the 
    %                               same energy content as the non-runaway 
    %                               part of f (eV)
    %   energy                   -- Energy moment (J/m^3)
    %   energyPerParticle        -- (eV)
    %   energyOfNonRelMaxwellian -- 3/2 n*T (J/m^3)
    %   momentum                 -- Parallel momentum moment (kg m/(s^2 m^3))
    %   runawayEnergyFraction    -- Fraction of total energy content of the 
    %                               distribution carried by particles in 
    %                               the runaway region p>p_c
    %   runawayFraction          -- Fraction of particles in the runaway
    %                               region
    %   runawayGrowthRate        -- n^{-1} dn_r/dtau
    %   runawayRegion            -- struct with information about the
    %                               runaway region separatrix used in the
    %                               final time step.
    %
    % %%% Inductive electric field %%% 
    %
    %   EHatInductive -- The electric field calculated consistently given
    %                    the applied electric field (EHat), the inductance
    %                    (LHat), and the distribution evolution.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Settings, parameters and interface methods %%%%%%%%%%%%%%%%%%%%%%%%
    
    properties %%% Settings %%%
        
        %Time advance
        timeAdvanceMode = 1 
        GMRESTolerance  = 1e-12
        nStepsBetweenFactorizations = 2
        dtIncreaseLimit = 0
        nNewtonSteps = 2
        
        %Collision operator
        collisionOperatorMode     = 0
        potentialCutOff           = 1e-30
        specialFunctionEvaluation = 1
        
        %Program flow and behavior 
        initialDistribution   = 0
        maxwellianShift       = 3
        runawayRegionMode     = 2
        abortWhenOnlyRunaways = 0
        GeneralAbortFunction
        
        %Grid
        pGridMode      = 1
        xiGridMode     = 1
        pGridParameter = 0.2;
                
        %Sources and sinks
        includeHeatSink  = 0 
        heatSinkCutOff   = 5 
        heatSinkFunction = 0
        enforceStrictHeatConservation = 0
        maxHeatSinkRate = 0
        keepTempConstantWhenDensityChanges = 1 
        dampeningStrength = 0
        dampeningWidth = 1.5
        
        %Saving and output
        nSaveSteps = 0
        silent     = false
        figOffset  = 1000
        timeStepSaveMode    = 0
        savePotentials      = 0
        show1DTimeEvolution = 0        
    end
    
    properties %%% Physical and numerical parameters %%% 
        nP
        nXi
        nL
        pMax
        dt
        tMax        
        
        T                
        n               
        Z
        EHat        
        B = 0 
        
        LHat = 0
        
        referenceTime = 0        
    end    
    
    methods %%% Main program components / interface functions %%% 
        function o = NORSE(varargin)
            % Constructor.
            %
            % Usage:
            %   NORSE() -- Initializes an empty object, does no calculation
            %
            % Use the provided settings and perform the calculation
            % automatically:
            %   NORSE(nP,nXi,nL,pMax,dt,tMax,T,n,Z,EHat)
            %   NORSE(nP,nXi,nL,pMax,dt,tMax,T,n,Z,EHat,B)
            %            
                        
            switch nargin
                case 0
                    o.grid = Grid(); %Empty grid object
                case {10,11} 
                    o.SetParameters(varargin{:});
                    o.PerformCalculation(); 
                otherwise
                    error('Wrong number of input arguments');
            end                        
        end
        
        SetParameters(o,nP,nXi,nL,pMax,dt,tMax,T,n,Z,EHat,varargin)
            %This function is handy for setting several parameters in one
            %call. The physical parameters can be scalars or
            %TimeDependentParameter objects.
        SetTimeDependentParameters(o,Tt,Tv,nt,nv,Zt,Zv,...
                                              Et,Ev,Bt,Bv)
            %Convenient method for setting all time-dependent parameters in
            %one call, if time and data vectors are available. The units
            %are: T: eV, n: m^(-3), E: V/m, B: T.
            %NOTE that the property referenceTime should be set before this
            %function is called, since its value is used here.        
        PerformCalculation(o, varargin)
            %Carries out the various steps in a standard NORSE calculation.            
        ContinueCalculation(o)
            % Restart/continue the calculation. Use existing settings
            % and/or time-dependent parameters (the user may change
            % settings before this call).
        Rewind(o,iSave)
            % Make an intermediate save step the final step, i.e. "rewind"
            % the calculation to a previous state. 
        Initialize(o)
            % Initializes variables and sets up parameters and flags.            
        AdvanceInTime(o,varargin)
            % Interface method for the various underlying time-advance
            % schemes.                
        PostProcess(o,varargin)
            % Calculates moments of the distribution and other
            % post-processing work.            
        CleanUp(o)
            % Clear large internal variables that are not needed when
            % saving the NORSE object. The object can still be used to
            % continue the calculation (after calling
            % RebuildInternalVariables).       
        AggressiveCleanUp(o)
            % Keep only the bare minimum which still makes the plots etc.
            % work, to reduce memory consumption. The object cannot be used
            % to continue the calculation!
        RebuildInternalVariables(o)
            % Re-initializes operators and objects that have been cleared
            % using CleanUp(). Makes it possible to restart the
            % calculation.
    end
    
    methods %%% Other, misc. %%%
        varargout = GetRelativisticNormalization(o,tMax,dt,T,varargin)
            %Converts tMax & dt given in thermal collision times to the
            %relativistic collision times used in NORSE.         
        LHat = NormalizeInductance(o,L,a,R)
            % Normalizes the inductance L in Henry to the internal LHat
            % used in NORSE.
        lnLambda = CalculateLnLambda(o,varargin) 
            %Calculates the Coulomb logarithm. The formula used is taken
            %from the NRL formulary.                       
        Ec = CalculateEc(o,varargin)
            % Calculates the critical field of Connor & Hastie [NF, 15
            % (1975) 415] in V/m.            
        nu = CalculateRelativisticCollFreq(o,varargin)
            % Calculates the relativistic collision frequency nu.            
        [EEc,EED] = CalculateEffectiveEFields(o,varargin)
            % Uses the effective temperature based on the energy content of
            % the bulk of the distribution to calculate effective values
            % for E/E_c(T_eff) and E/E_D(T_eff).             
        GetSize(o,varargin) 
            % Calculates the size in memory of the NORSE object by looping
            % through the fields and summing up the variable sizes, and
            % prints it to the console. Also calculates the size of the
            % other objects associated with a NORSE calculation.
    end
    
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties %%% Internal properties %%%
        
        %Objects and structs
        grid
        potentials
        collisionOperator
        heatSink
        particleSource
        timeAdvance
        plot
        
        %%%%%%%%%%%%%%%%%
        
        %Aggregate quantities
        Theta
        kappa
        fM0
        sigma
        nBar
        lnLambdaBar
        ThetaBar
        kappaBar
        maxwellianPreFactor
        
        %%%%%%%%%%%%%%%%%
                
        %Legendre mode mapping
        legModesToBigVectorMap
        bigVectorToLegModeMap
        
        %%%%%%%%%%%%%%%%%

        %Elements of the kinetic equation        
        identityWithBoundaryConditions
        EFieldOperator                         
        synchrotronOperator        
                
        %%%%%%%%%%%%%%%%%
                
        %Times and time steps        
        times        
        nTimeSteps        
        
        %%%%%%%%%%%%%%%%%
                        
        %Distribution and potentials
        f
        Upsilon0
        Upsilon1
        Upsilon2
        Pi      
        
        %Moments of the distribution
        density
        energy
        energyPerParticle
        energyOfNonRelMaxwellian
        effectiveTemperature
        effectiveBulkTemperature
        bulkDensity
        momentum
        currentDensity
        runawayFraction
        runawayEnergyFraction
        runawayGrowthRate
        runawayRegion
        
        %Source and sink magnitudes
        energyChangeMagnitudes
        densityChangeMagnitudes
        
        %Inductive electric field        
        EHatInductive
        
        %%%%%%%%%%%%%%%%%
                
        %Misc.
        timing
        timeEvoLine        
        isRestart
    end
    
    properties (Constant) %%% Physical constants %%%
        constants = struct('c',   2.997925e8,...
                           'e',   1.602176e-19,...
                           'm',   9.109380e-31,...
                           'eps0',8.854188e-12);
    end
    
    methods %%% Components of the kinetic equation %%% 
        AssembleEFieldOperator(o)
            % Builds the operator describing the electric field
            % acceleration.
        AssembleSynchrotronOperator(o)
            % Builds the operator describing synchrotron-radiation-reaction
            % losses.
    end
    
    methods %%% Mapping to Legendre modes %%%
        BuildMappingMatrices(o)
            % Constructs matrices for mapping between the 2D
            % finite-difference grid and finite-difference--Legendre-mode
            % representations.
        function fls = MapBigVectorToLegModes(o,f)
            % Maps a quantity defined on the 2D finite-difference grid to
            % the finite-difference--Legendre-mode space
            %
            % Usage:
            %   fls = MapBigVectorToLegModes(f)    

            fls = o.bigVectorToLegModeMap*f;
        end        
        function f = MapLegModesToBigVector(o,fls)
            % Maps a quantity in the finite-difference--Legendre-mode space
            % to the 2D finite-difference grid.
            %
            % Usage:
            %   f = MapLegModesToBigVector(fls)
            
            f = o.legModesToBigVectorMap*fls;             
        end
    end
    methods (Static) %   ...   %%%
        outPls = LegendrePolynomials(l,x)
            % Calculates the legendre polynomials P_i(x) for i=0,1,...,l
            % using Bonnet's recursion formula
    end
    
    methods %%% Runaway region %%%
        varargout = CalculateRunawayFraction(o,iteration)
            % Defines the runaway region and uses it to calculate the
            % runaway fraction of the distribution at a given time. Also
            % optionally returns the fraction of the total energy carried
            % by the runaways.
        varargout = FindParallelPCrit(o,varargin)
            % Finds where dp/dt at xi=1 vanishes -- this is the end of the
            % separatrix. Optionally returns the sum of forces at all p on
            % the parallel grid (at xi=1) as the array dpdt.
        [xiS,pS] = IntegrateSeparatrix(o,pc,iteration)
            % Strating from p=pc, integrate 'backwards' in xi (from xi=1 to
            % xi=-1) to trace out the entire separatrix.
        rate = CalculateRunawayGrowthRate(o)
            % Calculates the runaway growth rate n^{-1} dn_r/dt from the
            % runaway fraction at the saved time steps.
    end    
    methods (Static) %   ...   %%%
        dpdx = SeparatrixIntegrationEquation(xi,p,TNorm,EHat,...
                                                        sigma,dPidp,dPidxi)
            % Describes the equation to integrate in the calculation of the
            % separatrix of the runaway region. The inputs are the
            % corresponding quantites in the method IntegrateSeparatrix.
        status = InterruptIntegration(p,flag,pMax)
            % Aborts the calculation of the separatrix if the calculation
            % is close to reaching the end of the p grid. 
    end
    
    methods %%% Plotting, output %%%
        Print(o,str,varargin)
            %Prints text to the console if silent mode is not
            %enabled. Uses the syntax of fprintf.
        PrintTimings(o)
            % Prints information on the time used to perform various tasks
            % during the NORSE run. 
        ResetTimings(o)
            % Resets the various timers
        Initialize1DTimeEvolutionPlot(o)
            % Initializes a plot for showing a parallel cut of the
            % distribution during runtime.
        Refresh1DTimeEvolutionPlot(o,f)
            % Update the plot showing the parallel distribution during
            % runtime. 
    end
    
    methods %%% Other, misc. %%%
        ProcessPhysicalParameters(o)
            % Verifies the provided physical parameters and calculates
            % aggreagate quantities needed in NORSE.
        InitializeTAandC(o)
            % Initializes the proper TimeAdvance and CollisionOperator
            % objects, depending on the settings.
        TEff = CalculateEffectiveTemperature(o,energy,density)
            % Determines an effective temperature by finding the Maxwellian
            % with the same energy content as the distribution. 
        function delete(o) %%% Destructor
            % Destructor. Clears large variables from memory, since Matlab 
            % doesn't do it well.
            %
            % Usage:
            %   delete()
            
            if ~isempty(o.timeAdvance)
                %If the object has not really been used, we do not need to
                %be exessive with the clean-up (which will give errors)
                o.AggressiveCleanUp();
            end

            delete(o.grid);            
            delete(o.plot);
            delete(o.heatSink);
            delete(o.particleSource);
            delete(o.timeAdvance);
                        
            o.energyChangeMagnitudes = [];            
                   
            o.times = [];              
            
            o.f = [];            
            o.Pi = [];   
            
            o.density = [];
            o.energy = [];
            o.energyPerParticle = [];
            o.energyOfNonRelMaxwellian = [];
            o.momentum = [];
            o.currentDensity = [];
            o.runawayFraction = [];
            o.runawayGrowthRate = [];
            o.runawayRegion = [];                              
        end
    end
    
    methods (Access = protected)
        oCopy = copyElement(o)
            % Overrides Matlab's standard functionality for making
             % (shallow) copies of handle objects to ensure that copies of
             % the referenced handle objects are also made.
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
