%% eisFitController
% Steven A. Policastro, Ph.D. 
% Center for Corrosion Science and Engineering, 
% U.S. Naval Research Laboratory
% 4555 Overlook Avenue SW
% Washington, DC 20375 
% 
%  
% This function calls the specified equivalent circuit model that will be 
% fit to the data and passes the initial parameter set estimate to the 
% fitting function.  It also establishes the constratint limts on the fit 
% parameters for the simplex algorithm and calls the fitting function.  
% Lastly, it stores the results of the fit in output files.  Fit parameter 
% values are output to the command window.
% 
% Record of Code Revisions:
%
% * Created - June 2024
% * Last revision: 12-July-2024
%% Function definition
% 
% Function inputs
% 
% * T = instance of eisdata class
% * fType = name of the equivalent circuit model to be used to analyze the
% data
% * beta0 = array of values for the initial guess for the fit
% * fitDirectory = directory where fit results are stored
% * legendForPlotString = string of descriptions for the data on the plot
%
% Function outputs 
% 
% * Files with fit values and list of fit parameters to the command window
function eisFitController(T,fType,beta0,fitDirectory,legendForPlotString)
%%% Branch execution on equivalent circuit type
% This extended switch statement determines the equivalent circuit fit
% function, the parameters that will be altered to obtain the fit, the
% constraint limits on the parameters, and initializes the output string
    switch char(fType)                    

        case 'UndamagedCoating'
            fun = @matrixImpedance_UndamagedCoatingPlain;
            activeParams = [1,1,1]';
            constraintLimits = [ ...
                1.0e-1, 1.0e3; ...
                1.0e1, 1.0e12; ...
                1.0e-12, 1.0e-1 ...
                ];

            s0 = "Results from the fit algorithm"+ newline;
            s00 = "=============================="+ newline;
            s1 = "Rs = %d" + newline;
            s2 = "Rp = %d" + newline;
            s3 = "C = %d" + newline;
            s7 = s00 + s0 + s00 + s1 + s2 + s3 + s00;                 

        case 'ModifiedUndamagedCoating'
            fun = @matrixImpedance_ModifiedUndamagedCoatingPlain;
            activeParams = [1,1,1,1]';  
            constraintLimits = [ ...
                1.0e-1, 1.0e3; ...
                1.0e-1, 1.0e12; ...
                1.0e-12, 1.0e-1; ...
                0.4 1.0];
            s0 = "Results from the fit algorithm"+ newline;
            s00 = "=============================="+ newline;
            s1 = "Rs = %d" + newline;
            s2 = "Rp = %d" + newline;
            s3 = "Y0 = %d" + newline;
            s4 = "α = %d" + newline;
            s7 = s00 + s0 + s00 + s1 + s2 + s3 + s4 + s00;           

        case 'Randles'
            fun = @matrixImpedance_RandlesPlain;
            activeParams = [1,1,1,1,1]';
            constraintLimits = [  ...
                1.0e-1, 1.0e3; ...
                1.0e1, 1.0e12; ...
                1.0e-12, 1.0e-1; ...
                1.0e-1,1.0e12; ...
                1.0e-2, 1.0e2 ...
                ];     
            s0 = "Results from the fit algorithm"+ newline;
            s00 = "=============================="+ newline;
            s1 = "Rs = %d" + newline;
            s2 = "Rp = %d" + newline;
            s3 = "C = %d" + newline;
            s5 = "σ = %d" + newline;
            s6 = "B = %d" + newline;
            s7 = s00 + s0 + s00 + s1 + s2 + s3 + s5 + s6 + s00;           

        case 'ModifiedRandles'
            fun = @matrixImpedance_ModifiedRandlesPlain;
            activeParams = [1,1,1,1,1,1]';
            constraintLimits = [  ...
                1.0e-1, 1.0e3; ...
                1.0e-3, 1.0e12; ...
                1.0e-12, 1.0e-1; ...
                0.4, 1.0; ...
                1.0e-1,1.0e12; ...
                1.0e-2, 1.0e2 ...
                ];
            s0 = "Results from the fit algorithm"+ newline;
            s00 = "=============================="+ newline;
            s1 = "Rs = %d" + newline;
            s2 = "Rp = %d" + newline;
            s3 = "Y0 = %d" + newline;
            s4 = "α = %d" + newline;
            s5 = "σ = %d" + newline;
            s6 = "B = %d" + newline;
            s7 = s00 + s0 + s00 + s1 + s2 + s3 + s4 + s5 + s6 + s00;

        case 'NestedRandlesCoatingDefect'
            fun = @matrixImpedance_NestedRandlesCoatingDefectPlain;
            activeParams = [1,1,1,1,1]'; 
            constraintLimits = [  ...
                1.0e-1, 1.0e3; ...
                1.0e1, 1.0e12; ...
                1.0e-12, 1.0e-1; ...
                1.0e1, 1.0e12; ...
                1.0e-12, 1.0e-1; ...
                ];
            s0 = "Results from the fit algorithm"+ newline;
            s00 = "=============================="+ newline;
            s1 = "Rs = %d" + newline;
            s2 = "Rc = %d" + newline;
            s3 = "Cc = %d" + newline;
            s4 = "Rp = %d" + newline;
            s5 = "Cdl = %d" + newline;
            s7 = s00 + s0 + s00 + s1 + s2 + s3 + s4 + s5 + s00;           

        case 'REAP'
            fun = @matrixImpedance_REAPPlain;
            activeParams = [1,1,1,1,1,1]';  
            constraintLimits = [  ...
                1.0e-1, 1.0e3; ...
                1.0e-1, 1.0e12; ...
                1.0e-12, 1.0e-1; ...
                1.0e-1,1.0e12; ...
                1.0e-12, 1.0e-1; ...
                0.4, 1.0; ...
                ];
            s0 = "Results from the fit algorithm"+ newline;
            s00 = "=============================="+ newline;
            s1 = "Rs = %d" + newline;
            s2 = "Rpo = %d" + newline;
            s3 = "Cc = %d" + newline;
            s4 = "Rp = %d" + newline;
            s5 = "Y0 = %d" + newline;
            s6 = "α = %d" + newline;
            s7 = s00 + s0 + s00 + s1 + s2 + s3 + s4 + s5 + s6 + s00;           

        case 'ModifiedREAP'
            fun = @matrixImpedance_ModifiedREAPPlain;
            activeParams = [1,1,1,1,1,1,1]';
            constraintLimits = [  ...
                1.0e-1, 1.0e3; ...
                1.0e-1, 1.0e12; ...
                1.0e-12, 1.0e-1; ...
                0.4, 1.0; ...
                1.0e-1,1.0e12; ...
                1.0e-12, 1.0e-1; ...
                0.4, 1.0; ...
                ];

            s0 = "Results from the fit algorithm"+ newline;
            s00 = "=============================="+ newline;
            s1 = "Rs = %d" + newline;
            s2 = "Rpo = %d" + newline;
            s3 = "Y0c = %d" + newline;
            s4 = "αC = %d" + newline;
            s5 = "Rp = %d" + newline;
            s6 = "Y0dl = %d" + newline;
            s8 = "αdl = %d" + newline;
            s7 = s00 + s0 + s00 + s1 + s2 + s3 + s4 + s5 + s6 + s8 + s00;             

        case 'PolynomialTestPlain'
            fun = @PolynomialTestPlain;
            activeParams = [1,1]';
            constraintLimits = [  ...
                1.0e-9, 1.0e3; ...
                1.0e-9, 1.0e12; ...
                ];

    end
%%% Instantiation of simplexFit class
% An instance of the simplexFit class is created with the fit function and
% constraint limits passed to the constructor.
% Some local variables are then defined and then the a check is made to
% determine if the parallelization toolbox is installed.  If it is
% available, a parallel for loop is used.  If not, a sequential for loop is
% used.  In both cases, the simplex fit algorithm is called multiple times
% because the initial simplex is generated using the initial guess and some
% randomly generated points near the initial guess.  The fit parameters
% that return with the lowest mean-square error are output and used to
% estimate the fit impedance that is plotte.d
    fitClass = simplexFit(fun,constraintLimits);   
    fun2 = @fitClass.fitFn;
    check = ver('parallel');

    if  ~isempty(check)
        numIters = 10*numel(activeParams);
        bFits = zeros(numIters,length(beta0));
        mseVals = zeros(numIters,1);        
        parfor iter = 1:numIters
            % ========================
            % Fitting routine called
            [mdlEIS1,mseVals(iter,1),~] = fun2([T.Freq,T.Zreal,T.Zimag],beta0,activeParams);
            % ========================
            bFits(iter,:) = mdlEIS1.coefficients(:);            
        end
    else
        numIters = numel(activeParams) + 1; 
        bFits = zeros(numIters,length(beta0));
        mseVals = zeros(numIters,1);        
        % ========================
        for iter = 1:numIters
            % ========================
            % Fitting routine  called
            [mdlEIS1,mseVals(iter,1),~] = fun2([T.Freq,T.Zreal,T.Zimag],beta0,activeParams);
            % ========================
            bFits(iter,:) = mdlEIS1.coefficients(:);
        end        
    end

    [~,iMSE] = min(mseVals); 
    beta1 = bFits(iMSE,:);
    fprintf(s7,beta1);
    [zMod1,~,~,zPhase1,~,~] = fun(beta1,T.Freq);
    [zMod0,~,~,zPhase0,~,~] = fun(beta0,T.Freq);
    outputString = "Analysis of EIS data completed.";
    disp(outputString)
    outputString = newline;
    disp(outputString)      
    plotBodeEIS([T.Freq,T.Freq],[T.Zmod,zMod1],[T.Zphz,zPhase1],legendForPlotString)
    % ========================
    % Write the fit values and data to an output file
    % ========================
    Aexp = [T.Freq,T.Zmod,T.Zphz, zMod0, zPhase0, zMod1, zPhase1];
    oN = strcat(legendForPlotString,'.csv');
    outputName = fullfile(fitDirectory, oN);
    writematrix(Aexp, outputName)      
end
%------------- END OF CODE --------------