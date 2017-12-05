classdef Plot < matlab.mixin.Copyable
    % PLOT -- Class that implements various plots and visualizations for
    %         the NORSE solution or its properties.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Usage: 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialize:
    %   Plot()       -- empty object
    %   Plot(oNORSE) -- call to use in NORSE     
    %
    % After the NORSE calculation has finished, any number of plots can be
    % generated by calling the methods of this class. 
    % 
    % Example:
    %   oN.plot.Dist1D(),
    % where oN is the NORSE object, plots the distribution in the parallel
    % direction in the final time step.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    %%% Interface methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    methods %%% Plots of the distribution %%% 
        varargout = Dist1D(o,varargin)
            % Plot of the distribution in the parallel direction. 
        Dist2D(o,varargin)
            % Generates a contour plot of the distribution in both (p,xi)
            % and (p_para,p_perp) space.
        DistMovie1D(o,varargin)
            % Plays a "movie" of the distribution evolution in the parallel
            % direction. 
        DistMovie2D(o,varargin)
            % Plays a "movie" of the 2D distribution evolution in
            % (p_para,p_perp) space.
        DistVsTime(o,varargin)
            % Contour plot of the distribution in the parallel direction as
            % a function of time.
        Dist1DAtXi(o,idT,idXi,varargin)
            % Plots the distribution at a specified xi (id on the xi-grid)
            % at a given time step. 
        LegendreModes(o,varargin)
            % Plots a number of Legendre modes of the distribution.
    end
   
    methods %%% Information about the solution %%% 
        Moments(o)
            % Plots various moments of the distribution, showing for
            % instance the conservation properties and runaway generation.
        ParameterEvolution(o)
            % Compares the prescribed temperature and density evolution to
            % the actual ones obtained in NORSE.
        HeatSink(o,varargin)
            % Shows the magnitude of the heat sink, as well as the
            % contributions from its various components, as functions of
            % time. Also shows the momentum-space shape of the heat sink in
            % the final time step. 
        GMRESInfo(o)
            % Plots information about the performance of GMRES in the
            % iterative and adaptive time-advancement schemes. 
        TimeStep(o)
            % Plots the time step used by the adaptive time step scheme, as
            % a function of time (normalized to the initial time step).
    end
    
    methods %%% Other plots, misc. %%%        
        EffectiveEfields(o)
            % Plots the effective E/E_c and E/E_D values during the NORSE
            % run, based on the effective temperature of the distribution.
            % Also shows the effective temperature for comparison.
        Separatrix(o,varargin)
            % Overlays the 2D distribution and the various definitions of
            % the runaway-region separatrix in momentum space. 
        SeparatrixMovie(o,varargin)
            % Shows a movie overlaying the evolution of the 2D distribution
            % and the runaway-region separatrix. 
    end
   
    %%% Internal properties and methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties 
        norse
        grid 
    end
    
    methods %%% Constructor %%%
        function o = Plot(varargin)
            % Constructor.
            %
            % Usage:
            %   Plot() -- empty object
            %   Plot(oNORSE) -- call to use in NORSE 
            
           switch nargin
               case 0
                   %Nothing to do
               case 1
                   if isa(varargin{1},'NORSE')
                       o.norse = varargin{1};
                       o.grid = o.norse.grid;                       
                   else
                       error('The first argument must be a NORSE object.'); 
                   end
               otherwise
                   error('Invalid number of input arguments');
           end
           warning off MATLAB:Axes:NegativeDataInLogAxis            
       end              
    end
    
    methods %%% Helper functions %%% 
        PlotAtTimeX(o,~,eventdata)
            % Plots the 1D distribution at the time selected with the
            % mouse, when the time coordinate is on the x axis.
        PlotAtTimeY(o,~,eventdata)
            % Plots the 1D distribution at the time selected with the
            % mouse, when the time coordinate is on the y axis.
        PlotAtXiCyl(o,idT,~,eventdata)
            % Make a 1D plot at a given xi, specified by a click in a 2D 
            % plot in (p_para, p_perp) space.
        PlotAtXiSph(o,idT,~,eventdata)
            % Make a 1D plot at a given xi, specified by a click in a 2D 
            % plot in (p,xi) space.
        FunctionOfTime(o,vals,labelT,varargin)
            % Plots a vector as a function of time, with some convenient
            % labels and additions. The length of the vector must coincide
            % with the number of saved time step.
        id = GetFigId(o,str)
            % Returns an appropriate figure number, determined by the
            % calling function and the figure offset of the NORSE object.
        figOffset = GetFigOffset(o)
            % Returns the figure offset specified in the NORSE object.
    end
    
    methods (Static) %%% Color map %%%
        map = ColorMap(varargin)
            % Defines a color map that is suitable for both color and b&w
            % printing.
    end
end
