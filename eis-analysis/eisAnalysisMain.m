%% Electrochemical Impedance Spectroscopy (EIS) Analysis
% Steven A. Policastro and Rachel M. Anderson 
% Center for Corrosion Science and Engineering, 
% U.S. Naval Research Laboratory
% 4555 Overlook Avenue SW
% Washington, DC 20375 
% 
% This document serves as a summary and README for the classes and
% functions that were developed to analyze EIS data.  This project was 
% developed to facilitate the analysis of electrochemical impedance 
% spectroscopy (EIS) data obtained from potentiostatic EIS experiments.  
%
% *Note*: The code has not been tested to determine if it works with 
% galvanostatic EIS data.
%
% The main function depends on m files in the following sub-directories:
%
% * "Circuit Element Models" 
% * "Equivalent Circuit Models" 
% * "Fitting Algorithms" 
% 
% Data files are expected to be found in the "Data" directory. 
% Fit results are stored in the "Fit" directory. 
% 
% All functions and classes are written in MATLAB(R).
%
% Contact info: steven.a.policastro.civ@us.navy.mil 
%
% Record of Revisions: 
% Created - June 2024
%
% Last revision: 12-July-2024
%
% Contents
%
%% Main Function
%
% This function serves as the entry point for the project to analyze and 
% fit EIS data. It can be customized depending on the physical system the
% EIS data was collected from, but it needs to perform the following
% actions:
%
% * Clear the command window, all figures, and all variables from the 
% workspce.
% * Add the paths to the sub-directories containing the functions for the
% equivalent circuit impedance calculations and other supporting classes.
% * Create a cell array of data filenames
% * Iterate through the filenames to extract the impedance data.
% * Instantiate an instance of the eisFitController class for each
% datafile, pass the data to the constructor as well as the equivalent
% circuit type, and initial guess of the fit parameter values 
% 
% Definitions for the variables in the Main function are provided below:
%
% * datafilenames = cell arrray containinng filenames of the raw Gamry data
% files 
% * selectedEquivalentCircuit = character vector specifying the equivalent 
% circuit to be used to fit the data.  Available circuits are listed in the 
% following sub-section
% * vectorOfInitialParameterEstimates = array of values for the initial
% parameter estimates
% * legendString = character array of the descriptions of the data to be
% plotted on the output plots
% 
% 
%%% Available Equivalent Circuits:
% Specify the equivalent circuit impedance function to use for the
% fitting routine.
% 
% Circuit Name,Name to pass to EISFitController,Number of fit parameters
% needed in $\beta_0$ vector
% Undamaged Coating, UndamagedCoating,3
% Modified Undamaged Coating,ModifiedUndamagedCoating,4
% Randles,Randles,5
% Modified Randles,ModifiedRandles,6
% Nested Randles Coating Defect,NestedRandlesCoatingDefect,5
% Rapid Electrochemical Assessment of Paint (REAP),REAP,6
% Modified REAP,ModifiedREAP,7

function eisAnalysisMain    
    clc;
    clear all;
    addpath( ...
        'Circuit Element Models', ...
        'Equivalent Circuit Models', ...
        'Fitting Algorithms', ...
        'Data', ...
        'Fits' ...
        ) 
    format short
    base_data_dir = 'Data';
    base_fits_dir = 'Fits';    
    ext1 = '.dta';    
    ext2 = '.DTA';
    datafilenames = {'EIS_PureNi_Anodic_24hrOCP__495_Trial3'};
    legendString = 'AM - 504h NSW';

    for fn = datafilenames
        ffn1 = fullfile(base_data_dir,strcat(char(fn),ext1));
        ffn2 = fullfile(base_data_dir,strcat(char(fn),ext2)); 
        selectedEquivalentCircuit = 'ModifiedRandles';
        vectorOfInitialParameterEstimates = [1.0e1,3.0e4,1.0e-4,0.8,7.0e2,1.5];        
        good = 0;
        if isfile(ffn1)
            [~,eisdata] = AnalyzeGamryEISData(ffn1);
            good = 1;
        elseif isfile (ffn2)
            [~,eisdata] = AnalyzeGamryEISData(ffn2);
            good = 1;            
        end        
        if good == 1            
            eisFitController(eisdata,selectedEquivalentCircuit,vectorOfInitialParameterEstimates,base_fits_dir,legendString);
        else
            fprintf("File %s not found.",fn,"/n");            
        end
    end        
end
%% Principal Sub-Functions and Classes
% The principal supporting function that is used by the model is the 
% eisFitController function.  This function initiates and calls the
% important functions for performing the equivalent circuit fit of the EIS
% data and then displays the results.
%%