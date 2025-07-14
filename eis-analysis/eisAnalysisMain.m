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
% Revision: 27-June-2025
% Revision: 11-July-2025
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
% Circuit Name = Name to pass to EISFitController,Number of fit parameters
% needed in $\beta_0$ vector
%
% Available Equivalent Circuits                         Circuit Name                 Number of Free Parameters
% =============================                  ==============================      ========================
% Undamaged Coating                                 UndamagedCoating                             3
% Modified Undamaged Coating                        ModifiedUndamagedCoating                     4
% Randles                                           Randles                                      5
% Modified Randles                                  ModifiedRandles                              6
% Nested Randles Coating Defect                     NestedRandlesCoatingDefect                   5
% Rapid Electrochemical Assessment of Paint (REAP)  REAP                                         6
% Modified REAP                                     ModifiedREAP                                 7
% 
% 
% 
% For explanations of the circuit models for these equivalent circuits,
% please go to https://www.mdpi.com/2079-6412/13/7/1285.

function eisAnalysisMain    
    clc;
    clear all;
    
    format short
    % Add sub-folders where necessary support files are found
    addpath( ...
        'Circuit Element Models', ...
        'Equivalent Circuit Models', ...
        'Fitting Algorithms') 
    
    % ==============================================================
    % File handling
    % ==============================================================    
    % Define array of file extensions used for reading or writing files.
    exts = [".dta",".DTA",".xlsx",".csv"];
    % Define array of top directories where data can loaded and fit results
    % can be stored.  These should remain unchanged.
    dirsForDataAndFits = ["Data","Fits"];
    % Define sub-directories where data is located
    allSubdirectories =  "Coatings EIS";    
    basefilename = "7-3 48 hr";    
    % Determine if a file containing a list of data files exists by using
    % the file extension.  An xlsx file is assumed to contain a list of
    % filenames.
    if isfile(strcat(basefilename,exts(3)))
        % This type of file exists. 
        fnw = strcat(basefilename,exts(3));
        T = readtable(fnw); 
        datafilenames = string(T.file_name);               
    else 
        % File containig a list of data filenames does not exist, so just
        % treat each filename as a datafile
        datafilenames = strcat(basefilename,exts(1));
        lstring = {'Data from epoxy coating (7.3 48 h) on chromated Al'};
    end
    fns = numel(datafilenames);
    outputFitName = strcat(basefilename,exts(4));      
    
    % ==============================================================  
    % Specify the equvalent circuit model to be used for the fit
    % The initial parameter estimate vector must contain the same
    % number of parameter values as are free parameters in the circuit
    % model.  For exampls the Modified Randles circuit has 6 free
    % parameters, so the initial parameter estimate vector needs to
    % contain 6 initial guesses.  A more detailed explanation is provided
    % in the header information.
    % ==============================================================
    selectedEquivalentCircuit = "ModifiedREAP"; %"NestedRandlesCoatingDefect"; %"ModifiedRandles"; %""; %'';

    switch selectedEquivalentCircuit
        case "UndamagedCoating"
            varNames = {"Run","GOF","Rs","Rp", "Cdl"}; 
            vectorOfInitialParameterEstimates = [1.0e1, 1.0e3, 1.0e-4];
        case "ModifiedUndamagedCoating"
            varNames = {"Run","GOF","Rs","Rp","Y0dl","adl"}; 
            vectorOfInitialParameterEstimates = [1.0e1, 1.0e3, 1.0e-4, 0.8];
        case "Randles"
            varNames = {"Run","GOF","Rs","Rp", "Cdl","sigma","B"};
            vectorOfInitialParameterEstimates = [1.0e1, 1.0e3, 1.0e-4, 0.55, 7.0e2]; 
        case "ModifiedRandles"
            varNames = {'Run',"GOF",'Rs','Rp','Y0dl','adl','sigma','B'};
            vectorOfInitialParameterEstimates = [1.0e1,3.0e4,1.0e-4,0.55,7.0e2,1.5];            
        case "ModifiedRandles_SemiInfinite"
            varNames = {"Run","GOF","Rs","Rp", "Y0dl","adl","sigma"};
            vectorOfInitialParameterEstimates = [1.0e1,3.0e4,1.0e-4,0.8,7.0e2,1.5];  
        case "NestedRandlesCoatingDefect"
            varNames = {"Run","GOF","Rs","Rpo","Cc","Rp","Cdl"}; 
            % vectorOfInitialParameterEstimates = [1.0e2, 5.0e10, 1.0e-9, 1.0e4,8.0e-10];    
            vectorOfInitialParameterEstimates = [1.0e2, 2.0e7, 2.0e-10, 9.0e9, 4.0e-10]; 
        case "REAP"
            varNames = {"Run","GOF","Rs","Rpo","Cc","Rp","Y0dl","adl"}; 
            vectorOfInitialParameterEstimates = [1.0e1, 1.0e4, 1.0e-4, 1.0e2,1.0e-6, 0.8];                 
        case "ModifiedREAP"
            varNames = {"Run","GOF","Rs","Rpo", "Y0C","nC","Rp", "Y0DL","nDL"}; 
            vectorOfInitialParameterEstimates = [1.0e2, 2.0e7, 2.0e-10, 0.99, 9.0e9, 4.0e-10, 0.95];
    end
    ffvals = zeros(numel(datafilenames),numel(varNames));
    
    for fn_idx = 1:fns
        ffn1 = fullfile(dirsForDataAndFits(1),allSubdirectories,strcat(datafilenames(fn_idx)));
        legendStrings = {lstring{fn_idx},'Final fit','Initial guess'}; %lstring + "-" + num2str(fn_idx);
        if isfile(ffn1)
            [~,eisDataFromDatafile] = AnalyzeGamryEISData(ffn1);
            % Use this variable to specify which datapoints should be 
            % included in the fit - in case of noise or unclear physical
            % basis to justify more elements in the equivalent circuit.
            % Baseline assumption is that all data points will be included.
            allVals = numel(eisDataFromDatafile.Pt);
            valsToInclude = 1:1:allVals;            
            ffvals(fn_idx,:) = eisFitController(fn_idx,valsToInclude,eisDataFromDatafile,selectedEquivalentCircuit,vectorOfInitialParameterEstimates,dirsForDataAndFits(2),char(legendStrings));
        else
            fprintf("File %s not found. \n",char(ffn1));            
        end
    end  

    area = ones(size(ffvals(:,4))).*14.97; %0.495;
    T1 = table;
    switch selectedEquivalentCircuit
        case "UndamagedCoating"
            T1.Run = lstring';
            T1.GOF = ffvals(:,2);
            T1.Rs = ffvals(:,3);
            T1.Rp = ffvals(:,4);
            T1.Cdl = ffvals(:,5);
            T1.area = area;            
        case "ModifiedUndamagedCoating" 
            T1.Run = lstring'; 
            T1.GOF = ffvals(:,2);
            T1.Rs = ffvals(:,3);
            T1.Rp = ffvals(:,4);
            T1.Y0dl = ffvals(:,5);
            T1.adl = ffvals(:,6);
            T1.area = area; 
        case "Randles"
            T1.Run = lstring'; 
            T1.GOF = ffvals(:,2);
            T1.Rs = ffvals(:,3);
            T1.Rp = ffvals(:,4);
            T1.Cdl = ffvals(:,5);            
            T1.sigma = ffvals(:,6);
            T1.B = ffvals(:,7);
            T1.area = area;               
        case "ModifiedRandles_SemiInfinite"
            T1.Run = lstring'; 
            T1.GOF = ffvals(:,2);
            T1.Rs = ffvals(:,3);
            T1.Rp = ffvals(:,4);
            T1.Y0 = ffvals(:,5);
            T1.alpha = ffvals(:,6);
            T1.sigma = ffvals(:,7);
            T1.area = area;            
        case "ModifiedRandles"
            T1.Run = lstring'; 
            T1.GOF = ffvals(:,2);
            T1.Rs = ffvals(:,3);
            T1.Rp = ffvals(:,4);
            T1.Y0 = ffvals(:,5);
            T1.alpha = ffvals(:,6);
            T1.sigma = ffvals(:,7);
            T1.Bval = ffvals(:,8);
            T1.area = area;  
        case "NestedRandlesCoatingDefect"
            T1.Run = lstring'; 
            T1.GOF = ffvals(:,2);
            T1.Rs = ffvals(:,3);
            T1.Rpo = ffvals(:,4);
            T1.Cc = ffvals(:,5);
            T1.Rp = ffvals(:,6);
            T1.Cdl = ffvals(:,7);            
            T1.area = area;              
        case "REAP"
            T1.Run = lstring';
            T1.GOF = ffvals(:,2);
            T1.Rs = ffvals(:,3);
            T1.Rpo = ffvals(:,4);
            T1.Cc = ffvals(:,5);
            T1.Rp = ffvals(:,6);
            T1.Y0DL = ffvals(:,7);
            T1.alphaDL = ffvals(:,8);
            T1.area = area;              
        case "ModifiedREAP"
            T1.Run = lstring';
            T1.GOF = ffvals(:,2);
            T1.Rs = ffvals(:,3);
            T1.Rpo = ffvals(:,4);
            T1.Y0C = ffvals(:,5);
            T1.alphaC = ffvals(:,6);
            T1.Rp = ffvals(:,7);
            T1.Y0DL = ffvals(:,8);
            T1.alphaDL = ffvals(:,9);
            T1.area = area;                      
    end

    writetable(T1,outputFitName)

end
%% Principal Sub-Functions and Classes
% The principal supporting function that is used by the model is the 
% eisFitController function.  This function initiates and calls the
% important functions for performing the equivalent circuit fit of the EIS
% data and then displays the results.
%%
