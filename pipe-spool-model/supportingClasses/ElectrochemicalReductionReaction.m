classdef ElectrochemicalReductionReaction
    %ElectrochemicalReductionReaction - Class for modeling reduction
    %reactions
    %
    %The purpose of this class is to create an object that can be used to
    %model and manipulate cathodic polarization data and models
    %
    %==========================================================================
    % Author:   Steve Policastro, Ph.D., Materials Science
    % Center for Corrosion Science and Engineering, U.S. Naval Research
    % Laboratory
    % email address: steven.policastro@nrl.navy.mil  
    % Website: 
    % October 2021; Last revision: 27-Oct-2022
    %==========================================================================

    properties (SetAccess = private)
        c = Constants;
        maxIter = 1000;
        tol = 1.0e-6;
        tau = 1.0e-3;        
    end   
    properties (SetAccess = public)
        name reactionNames        
        Temperature double
        lambda_0 double   
        alpha double
        z double
        cOxidized double
        cReduced double
        diffusionLength double
        diffusionCoefficient double
        eApp double
        E0 double
        EN double
        eta double
        dG_cathodic double
        dG_anodic double
        i0_Cathodic double        
        i0_Anodic double   
        iCathodic double
        iAnodic double
        iAct double
        iLim double
        qRatio double
        i double
        plotSymbol string
    end

    methods
        function obj = ElectrochemicalReductionReaction(nameString, pS, T, c_ox, c_red, dG, alpha, del, z, e0, D)
            %ElectrochemicalReductionReaction Construct an instance of the ElectrochemicalReductionReaction class
            %The purpose of this method is to construct an instance of the
            %ElectrochemicalReductionReaction class.
            obj.name = nameString;
            obj.Temperature = T; % + Constants.convertCtoK
            obj.lambda_0 = (obj.c.kb*T)/obj.c.planck_h;           
            obj.dG_cathodic = dG(1);
            obj.dG_anodic = dG(2);
            obj.alpha = alpha;
            obj.plotSymbol = pS;
            obj.cOxidized = c_ox;
            obj.cReduced = c_red;
            obj.E0 = e0;

            obj.z = z;
            pre_factor = (obj.c.R*obj.Temperature)/(obj.z*obj.c.F); 
            EN_log = log((c_red(1) * c_red(2))/(c_ox(1) * c_ox(2)));
            obj.EN = obj.E0 + (pre_factor * EN_log); %V_SHE     
            obj.diffusionLength = del; %mm ?
            obj.diffusionCoefficient = D; %cm2/s ?
        end    
    
        function obj = setDiffusionLimitedCurrent(obj,delta)
            %Computes iLim
            switch objCatModel.name
                case 'HER'
                    num = objCatModel.z*objCatModel.c.F*objCatModel.diffusionCoefficient*(objCatModel.cOxidized(1)/Constants.M_H2);
                case 'ORR'
                    num = objCatModel.z*objCatModel.c.F*objCatModel.diffusionCoefficient*(objCatModel.cOxidized(1)/Constants.M_O2);
                case 'Fe_Red'
                    num = objCatModel.z*objCatModel.c.F*objCatModel.diffusionCoefficient*(objCatModel.cOxidized(1)/Constants.M_Fe);
            end
            
            obj.diffusionLength = delta;
            obj.iLim = -num/obj.diffusionLength; %A/cm2 
        end
        
        function [iL,diLdDel] = getDiffusionLimitedCurrent(obj)
            if obj.name == reactionNames.None
                iL = 0.0;
                diLdDel = 0.0;
            else
                %Computes iLim, and diLim/dDel
                iL = obj.iLim; %A/cm2   
                diLdDel = iL/obj.diffusionLength;
            end
        end

    end
    methods (Static)
        
        function [iL,diLdDel] = GetDiffusionLimitedCurrent(objCatModel)
            if objCatModel.name == reactionNames.None
                iL = 0.0;
                diLdDel = 0.0;
            else
                %Computes iLim, and diLim/dDel
                switch objCatModel.name
                    case 'HER'
                        num = objCatModel.z*objCatModel.c.F*objCatModel.diffusionCoefficient*(objCatModel.cOxidized(1)/Constants.M_H2);
                    case 'ORR'
                        num = objCatModel.z*objCatModel.c.F*objCatModel.diffusionCoefficient*(objCatModel.cOxidized(1)/Constants.M_O2);
                    case 'Fe_Red'
                        num = objCatModel.z*objCatModel.c.F*objCatModel.diffusionCoefficient*(objCatModel.cOxidized(1)/Constants.M_Fe);
                end
                
                iL = -num/objCatModel.diffusionLength; %A/cm2   
                diLdDel = num/(objCatModel.diffusionLength^2);
            end
        end

        function [iAct,diAdDGc,diAda,diAdz] = GetActivationCurrent(objCatModel, vApp) 

            if objCatModel.name == reactionNames.None
                diAdDGc = 0.0;
                diAda = 0.0;
                diAdz = 0.0;
                iAct = 0.0;
            else
                %Computes iAct, diAct/dDG, diAct/da, and diAct/dz
                eta = vApp - (objCatModel.EN - objCatModel.c.E_SHE_to_SCE);
                R = objCatModel.c.R;
                T = objCatModel.Temperature;
                F = objCatModel.c.F;
                RT = R*T;
                dGC = objCatModel.dG_cathodic;
                dGA = objCatModel.dG_anodic;
                z = objCatModel.z;
                l0 = objCatModel.lambda_0;
                pF = z*F*l0;
                i0C = pF*exp(-dGC/RT);
                i0A = pF*exp(-dGA/RT);
    
                objCatModel.eta = eta;
                objCatModel.i0_Cathodic = i0C;            
                objCatModel.i0_Anodic = i0A;
    
                a = objCatModel.alpha;
                expValC = -(((1.0-a)*z*F).*eta)./RT;
                expValA = ((a*z*F).*eta)./RT;
    
                objCatModel.iCathodic = i0C.*exp(expValC);
                objCatModel.iAnodic = i0A.*exp(expValA);
    
                diAdDGc = objCatModel.iCathodic.*(1.0/RT);
                diAda = ((z*F.*eta)./RT).*(objCatModel.iAnodic - objCatModel.iCathodic);
                
                diAdz = ((i0A/z).*exp(expValA)) + (i0A.*exp(expValA).*((a*F.*eta)./RT)) - ((i0C/z).*exp(expValC)) - (i0C.*exp(expValC).*((-(1-a)*F.*eta)./RT));
    
                iAct = objCatModel.iAnodic - objCatModel.iCathodic;
            end
        end        

        function i = GetKouteckyLevich(objCatModel)
            if objCatModel.name == reactionNames.None
                i = 0.0;
            else
                i = (objCatModel.iLim  .* objCatModel.iAct)./(objCatModel.iAct + objCatModel.iLim ); 
            end
        end
        
        function newB = LM_Fit_Activation2(objCatModel, yData, xData)
            
            B = [objCatModel.dG_cathodic, objCatModel.alpha];
            J = zeros(length(xData), length(B));
            k = length(B);

            y = yData; 

            lambda = 1.0e1; 
            lambdaUp = 3.0;
            lambdaDown = 4.0;
            kMult = 1.0;
%             trackMSE = zeros(objCatModel.maxIter,1);

            for iter = 1:objCatModel.maxIter
                testIterSum = 0;
                
                objCatModel.dG_cathodic = B(1);
                objCatModel.alpha = B(2);   
                [i0_Cathodic,didg,dida,~] = ElectrochemicalReductionReaction.GetActivationCurrent(objCatModel,xData);

                mse0 = CalculateSSD(i0_Cathodic,y,k);

                J(:,1) = didg;
                J(:,2) = dida;

                eK = y - i0_Cathodic;    
                A = J'*J;
                g = J'*eK;
                
                aDiag = diag(A);
                aD = diag(aDiag);

                lambdaOld = lambda;
                marquardtVal = lambdaOld.*aD;
                Aprime = A + marquardtVal;

                delOld = Aprime\g;

                B0 = B;
                Bold(:) = B0(:) + (kMult.*delOld(:));

                objCatModel.dG_cathodic = Bold(1);
                objCatModel.alpha = Bold(2);
                [iOld,~,~,~] = ElectrochemicalReductionReaction.GetActivationCurrent(objCatModel,xData);

                mseOld = CalculateSSD(iOld,y,k);

                if lambda > 1.0e-5
                    lambdaNew = lambda/lambdaDown;
                    marquardtVal = lambdaNew.*aD;
                    Aprime = A + marquardtVal;
                    delNew = Aprime\g;
                    
                    B0 = B;
                    Bnew(:) = B0(:) + (kMult.*delNew(:)); 

                    objCatModel.dG_cathodic = Bnew(1);
                    objCatModel.alpha = Bnew(2); 
                    [iNew,~,~,~] = ElectrochemicalReductionReaction.GetActivationCurrent(objCatModel,xData);
                    mseNew = CalculateSSD(iNew,y,k);
                else
                    mseNew = mseOld;
                    iNew = iOld;
                end

                if mseNew <= mse0 || mseOld <= mse0
                    kMult = 1.0;
                    if mseNew < mseOld
                        lambda = lambdaNew;
                        B = Bnew;
                        
                        for i = 1:length(B)
                            testConvergence = abs(delNew(i))/(objCatModel.tau + abs(Bnew(i)));
                            if testConvergence < objCatModel.tol
                                testIterSum = testIterSum + 1;
                            end
                        end                     
                        if testIterSum == length(B)
%                             mse = mseNew;
                            break;
                        end                          
                    else
                        lambda = lambdaOld;
                        B = Bold;
                        for i = 1:length(B)
                            testConvergence = abs(delOld(i))/(objCatModel.tau + abs(Bold(i)));
                            if testConvergence < objCatModel.tol
                                testIterSum = testIterSum + 1;
                            end
                        end                        
                        if testIterSum == length(B)
%                             mse = mseOld;
                            break;
                        end                          
                    end                   
                else %if mseNew > mse0 && mseOld > mse0                    
                    lambda = lambda * lambdaUp;
                    if mseNew <= mseOld
                        kMult = 1.0e-1*(mse0/mseNew);
%                         trackMSE(iter,1) = mseNew;
                    else
                        kMult = 1.0e-1*(mse0/mseOld);
%                         trackMSE(iter,1) = mseOld;
                    end
                    if lambda >= 1.0e30
%                         mse = trackMSE(iter,1);

%                         figure(200)
%                         hold on
%                         plot(1:iter,trackMSE(1:iter),'-bo')                        
%                         ax = gca;
%                         ax.YScale = 'log';
%                         hold off

                        break;
                    end
                end                    
               
            end

            newB = B;

%             figure(800)
%             hold on
% %             plot(objCatModel.iCathodic,objCatModel.eApp,'--k')
%             plot(abs(yData),xData,'-k','LineWidth',2)
%             plot(abs(i0_Cathodic),xData,'-g','LineWidth',1)
%             plot(abs(iNew),xData,'-b','LineWidth',4)
%             plot(abs(iOld),xData,'-r','LineWidth',4)            
%             ax = gca;
%             ax.XScale = 'log';
%             hold off
%             disp(newB)

        end

        function newB = LM_Fit_Transport(objCatModel, yDataAll, xData, fP)
            
            yData = zeros(size(yDataAll));
            yMean = mean(yDataAll);
            yMax = max(yDataAll);

            switch fP
                case 1
                    yData(:) = yMean;
                case 2
                    yData(:) = yDataAll(end);
                case 3
                    yData(:) = yMax;
                otherwise
                    yData(:) = yMean;
            end
            
            y = -1.0.*yData;
            B = [objCatModel.diffusionLength];
            J = zeros(length(xData), length(B));
            maxIter = 500;
            kMult = 1.0;
            k = length(B);

            lambda = 1.0e2;
            lambdaUp = 3.0;
            lambdaDown = 4.0;

            for iter = 1:maxIter
                testIterSum = 0;
                [iL0,diLdDel] = ElectrochemicalReductionReaction.GetDiffusionLimitedCurrent(objCatModel);
                iL(1:length(yData)) = iL0;
                mse0 = CalculateSSD(iL,yData,k);

                J(:,1) = diLdDel;

                eK = y - iL';
    
                A = J'*J;
                g = J'*eK;
       
                aDiag = diag(A);
                aD = diag(aDiag);       

                lambdaOld = lambda;
                marquardtVal = lambdaOld.*aD;
                Aprime = A + marquardtVal;

                delOld = Aprime\g; 

                B0 = B;
                Bold(:) = B0(:) + (kMult.*delOld(:)); 

                objCatModel.diffusionLength = Bold(1);
                [iOld0,~] = ElectrochemicalReductionReaction.GetDiffusionLimitedCurrent(objCatModel);
                iOld(1:length(yData)) = iOld0;
                mseOld = CalculateSSD(iOld,y,k);

                lambdaNew = lambda/lambdaDown;
                marquardtVal = lambdaNew.*aD;
                Aprime = A + marquardtVal;

                delNew = Aprime\g; 

                B0 = B;
                Bnew(:) = B0(:) + (kMult.*delNew(:)); 

                objCatModel.diffusionLength = Bnew(1);
                [iNew0,~] = ElectrochemicalReductionReaction.GetDiffusionLimitedCurrent(objCatModel);
                iNew(1:length(yData)) = iNew0;
                
                mseNew = CalculateSSD(iNew,y,k);

                if mseNew < mse0 || mseOld < mse0
                    kMult = 1.0;
                    if mseNew < mseOld
                        lambda = lambdaNew;
                        B = Bnew;
                        
                        for i = 1:length(B)
                            testConvergence = abs(delNew(i))/(objCatModel.tau + abs(Bnew(i)));
                            if testConvergence < objCatModel.tol
                                testIterSum = testIterSum + 1;
                            end
                        end                     
                        if testIterSum == length(B)
%                             mse = mseNew;
                            break;
                        end                          
                    else
                        lambda = lambdaOld;
                        B = Bold;
                        for i = 1:length(B)
                            testConvergence = abs(delOld(i))/(objCatModel.tau + abs(Bold(i)));
                            if testConvergence < objCatModel.tol
                                testIterSum = testIterSum + 1;
                            end
                        end                        
                        if testIterSum == length(B)
%                             mse = mseOld;
                            break;
                        end                          
                    end                   
                else %if mseNew > mse0 && mseOld > mse0
                    lambda = lambda * lambdaUp;
                    if mseNew <= mseOld
                        kMult = 1.0e-1*(mse0/mseNew);
                    else
                        kMult = 1.0e-1*(mse0/mseOld);
                    end
                    if lambda >= 1.0e30


%                         figure(200)
%                         hold on
%                         plot(1:iter,trackMSE(1:iter),'-bo')                        
%                         ax = gca;
%                         ax.YScale = 'log';
%                         hold off

                        break;
                    end
                end                
            end

            newB = B;
%             disp(newB)
        end
        
        function iaf = NR_ia(ia0,iL,i)
            N = 100;
            tol = 1.0e-6;
            iaf = 0.0;

            for j = 1:N
                f = (ia0*iL) - (i*ia0) - (i*iL);
                fp = iL-i;
                ia1 = ia0 - (f/fp);

                check = abs(ia1 - ia0);
                if check < tol
                    iaf = ia1;
                    break;
                else
                    ia0 = ia1;
                end
            end
        end
    
        function iLf = NR_iL(iL0,ia,i)
            N = 100;
            tol = 1.0e-6;
            iLf = 0.0;

            for j = 1:N
                f = (iL0*ia) - (i*iL0) - (i*ia);
                fp = ia-i;
                iL1 = iL0 - (f/fp);

                check = abs(iL1 - iL0);
                if check < tol
                    iLf = iL1;
                    break;
                else
                    iL0 = iL1;
                end
            end
        end        
            
    end

end