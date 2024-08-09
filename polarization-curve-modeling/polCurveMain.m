%% Polarization Curve Modeling
% Steven A. Policastro, Ph.D. 
% Center for Corrosion Science and Engineering, 
% U.S. Naval Research Laboratory
% 4555 Overlook Avenue SW
% Washington, DC 20375 
% 
%  
% This document serves as a summary and README for the set of physics-based
% and response surface models that are used to calculate simplified
% polarization curves for several alloy systems.  
% 
% Polarization curves can be obtained across a range of temperatures, 
% $Cl^-$ concentrations, pH values, and electrolyte flow velocities.  Note,
% though, that the simplified nature of the polarization curves arises
% becuase these are intended to serve as boundary conditions for FEM of
% complex parts and components exposed to corrosive environments.  The BCs
% on the electrode surfaces in those calculations are restricted to
% single-valued functions.  Thus, observed polarization behavior that
% results in active-passive transitions, for example, will not be captured
% by these models. 
% 
% All functions and classes are written in MATLAB(R).
%
% Contact info: steven.a.policastro.civ@us.navy.mil 
%
% Record of Revisions: 
% Created - October 2021
%
% Last revision: 27-Jul-24
%
% Contents
%
%% Main Function
% This example function serves as the entry point to the code for modeling 
% a polarization curve. It must accomplish the following actions:
%
% * Clear the command window, all figures, and all variables from the 
% workspce.
% * Instantiate an instance of the Constants class.
% * Define the potential region of interest
% * Define the various environmental conditions for generating the
% polarization curves
% * Iterate through the environmental conditions to define and plot the
% defined polarization curves.
% 
% Two examples of the polarization curves that can be generated are shown
% below the following main function:
function polCurveMain
    clc;
    clear all;
    C = Constants;
    vRange = -1.5:0.005:0.5; %VSCE    
    
    alloy = {'ss316'};
    T = [25.0,5.0]+ C.convertCtoK; %C
    pH = [7.0,3.0]; 
    cCl = [0.6,0.02]; %M
    velocity = [5.0,5.0]; %m/s
    envConds = [T;pH;cCl;velocity];

    for i = 1:size(envConds,2)               
        aPolCurve = PolarizationCurveModel(i,alloy,vRange, envConds(:,i));
        PolarizationCurveModel.Plot_Polarization_Data_and_Model(aPolCurve,600+i)
    end   
end
%% Principal Sub-Functions and Classes
% The principal supporting clases and sub-functions that are used in the
% model include the following: 
%
% * Constants - a class that defines a number of conversion values for
% converting units as well as a number of physical parameters for
% electrochemical calculations. 
% * PolarizationCurveModel - a class that creates an object to calculate 
% and plot a polarization curve.
%%
