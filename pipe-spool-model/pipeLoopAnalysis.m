%% Pipe Loop Analysis
% Steven A. Policastro, Ph.D. 
% Center for Corrosion Science and Engineering, 
% U.S. Naval Research Laboratory
% 4555 Overlook Avenue SW
% Washington, DC 20375 
% 
%  
% This document serves as a summary and README for the set of physics-based
% and response surface models that are used to calculate the galvanic
% interactions between 2 pipe spools composed of different materials and
% exposed to a user-specified electrolyte.  The pipe spools are simplified
% to 2-dimensional interactions by assuming axial and radial symmetry. 
% 
% The polarization curves that are used for the electrode boundary
% conditions in this work are modeled using the functions and procedures
% developed in the "Polarization Curve Modeling" project.
% 
% .
% All functions and classes are written in MATLAB(R).
%
% Contact info: steven.a.policastro.civ@us.navy.mil 
%
% _Record of Revisions_:
%
% * Created - january 2023
% * Last revision: 30-Jul-2024 
%
% Contents
% 
%% Main Function
% This example function serves as the entry point to the code for modeling 
% the galvanic interaction between two pipe spools. 
% It must accomplish the following actions:
%
% * Clear the command window, all figures, and all variables from the 
% workspce.
% * Ensure the code can find supporting classes in the respective
% sub-folders
% * Provide paths to the raw data and any other model results (if needed)
% * Instantiate a PipeLoopModel object
% * Define the properties of the pipe spools and the computational cell
% properties
% * Define the various environmental conditions for generating the
% polarization curves
% * Instantiate a GalvanicCorrosion object 
% * Iterate through the environmental conditions to define and plot the
% defined polarization curves.
function pipeLoopAnalysis 
    clc;
    clear all;
    addpath('plotFunctions')
    addpath('supportingClasses')
    addpath('materials')
    %==============================
    % Some definitions
    %==============================
    fn = 'Pipe Loop Data Only';
    fnP1 = 'PipeLoop4ElsycaModel';
    ext = '.csv';
    filename = fullfile(strcat(fn,ext));     
    PLE = PipeLoopExperiments(filename);
    fn2 = fullfile(strcat(fnP1,ext));  
    PLM = PipeLoopModel(fn2);
    PLM.hasMyModel = true;
    PLM.hasElsyca = true;    
    %=============================
    % Plots of potential along the pipes at the last time-stamp in the
    % files
    %============================= 
    endTimePlots = true;     
    cathodeWidth = 3.0; %m
    anodeWidth = 3.0; %m
    widthTotal = cathodeWidth + anodeWidth;     
    deltaDistance = 3.0e-2;%m
    aNodes = round(anodeWidth/deltaDistance); 
    cNodes = round(cathodeWidth/deltaDistance);
    totNodes = aNodes + cNodes;
    electrolyteNodes = 50;
    electrolyteHeight = electrolyteNodes * deltaDistance; 
    edgesBC = {'neumann','neumann','neumann'}; 
    vApp = 0.0; 
    aReactName = 'cuni'; 
    cReactName = 'i625'; 
    env = [0.6, 25.0, 9.0, 5.0];
    gC = galvanicCorrosion(edgesBC,widthTotal,electrolyteHeight,deltaDistance,electrolyteNodes,totNodes,aNodes,cNodes,aReactName,cReactName,vApp,12,env);
    gC.aSim.phi = galvanicCorrosion.JacobiSolver(gC.aSim);
    [PLM.pipeLoopPotentialSimpleModel,PLM.xSimple] = PipeLoopModel.GetData(PLM,gC.aSim);
    if endTimePlots == true
        numPipes = size(PLE.pipeLoopsAnodic,2); 
        endTimeIndex = numel(PLE.time);
        for i = 2:2 %1:numPipes
            cPotsExp = PLE.pipeLoopsCathodic(1,i);
            aPotsExp = PLE.pipeLoopsAnodic(1,i);
            corrosionPotentialsEndTime(i,endTimeIndex,cPotsExp,aPotsExp,PLM)
         end
    end
end
%% Principal Sub-Functions and Classes
% The principal supporting clases and sub-functions that are used in the
% model include the following: 
%
% * PipeLoopModel - a class that contains model predictions and
% experimental results for comparisons.
% * galvanicCorrosion - a class that contains functions for performing the
% potential distribution calculation for the galvanic couple between the
% pipe spools.
%%