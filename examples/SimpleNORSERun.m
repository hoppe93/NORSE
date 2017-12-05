%%% This script shows a simple example of a NORSE run. The runtime is
%%% approximately 20 seconds. The obtained final distribution
%%% (corresponding to the last column of the field o.f) is available in
%%% outputSimple.mat.
%%%
%%% The NORSE files must be in the Matlab path.
%%%  
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Written by Adam Stahl, 2017
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%% Resolution and parameters
%Physical parameters
T     = 1500; %eV
n     = 5e19; %m^{-3}
EHat  = 40;   %E/E_c
Z     = 1;    %Z_eff
B     = 3;    %T

%Numerical parameters
nP    = 200;
nXi   = 55;
pMax  = 1;
nL    = 15;
dt    = 3e-5; 
tMax  = 5e-3; 

%% Set up and run NORSE with default settings
o = NORSE(nP,nXi,nL,pMax,dt,tMax,T,n,Z,EHat,B);

%% Visualize the result
% See Plot.m for a complete list of plots

o.plot.Dist2D();
o.plot.DistVsTime();
o.plot.Moments();

%% Compare to the expected final distribution

expected = load('outputSimple.mat');
if max(abs(expected.f - o.f(:,end))) < 1e-12
    fprintf('\n  The resulting distribution agrees with the expected value!\n');
else
    fprintf('\n  The resulting distribution does not agree with the expected value!\n');
    fprintf('  Have you changed some setting? \n');
end
