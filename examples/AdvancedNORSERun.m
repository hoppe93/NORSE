%%% This script shows an example of a more involved NORSE run, including
%%% how to change settings, use time-dependent parameters, and continue the
%%% calculation after it has finished. The runtime is approximately 30
%%% seconds. The obtained final distribution (corresponding to the
%%% final column of the field o.f) is available in outputAdvanced.mat.
%%%
%%% The NORSE files must be in the Matlab path.
%%% 
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Written by Adam Stahl, 2017
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%% Resolution and parameters
%Constant physical parameters
n     = 5e19; %m^{-3}
EHat  = 6;    %E/E_c
Z     = 1;    %Z_eff
B     = 3;    %T

%Numerical parameters
nP    = 175;
nXi   = 35;
yMax  = 14;   % Thermal momenta (gamma v/v_th)
nL    = 7;
dt    = 0.2;  % These are now specified in thermal collision times!
tMax  = 40;   %       -- || --

%Define a time-dependent temperature
tT    = tMax*[0,0.2,0.7,1];
valsT = 1500*[1,1,1.03,1.02]; %eV
T     = TimeDependentParameter(tT,valsT);

%% Plot the temperature evolution
T.Visualize();

%% Set up NORSE

%Initialize a NORSE object
o = NORSE();       

%Change some settings (see NORSE.m for a complete list)
o.nSaveSteps                    = 30;   %0 = save the distribution at all timesteps
o.includeHeatSink               = 1; 
o.enforceStrictHeatConservation = 1;
o.show1DTimeEvolution           = 1;

%Convert times and grid maximum from 'thermal' units to NORSE units
%(relativistic collision times and p). The third argument means that the
%temperature at time t=0 is used to define v_th 
[tMax,dt,pMax,tScaleFactor] = o.GetRelativisticNormalization(tMax,dt,T(0),yMax);

%Make sure the temperature changes on the correct time scale
T = T.Rescale(tScaleFactor); 

o.SetParameters(nP,nXi,nL,pMax,dt,tMax,T,n,Z,EHat,B);

%% Run NORSE
o.PerformCalculation();

%% Visualize the result
% See Plot.m for a complete list of plots

o.plot.Dist2D();
o.plot.ParameterEvolution();
o.plot.HeatSink();

%% Resume evolving the system with a weaker electric field
o.EHat       = 0.5*o.EHat;
o.tMax       = 2*o.tMax;        %This is the total time!
o.nSaveSteps = 2*o.nSaveSteps;  %This is the total number of save steps!
o.ContinueCalculation();

%% Visualize the end result
o.plot.DistVsTime();
o.plot.Moments();
o.plot.ParameterEvolution();
o.plot.HeatSink();

%% Compare to the expected final distribution

expected = load('outputAdvanced.mat');
if max(abs(expected.f - o.f(:,end))) < 1e-12
    fprintf('\n  The resulting distribution agrees with the expected value!\n');
else
    fprintf('\n  The resulting distribution does not agree with the expected value!\n');
    fprintf('  Have you changed some setting? \n');
end
