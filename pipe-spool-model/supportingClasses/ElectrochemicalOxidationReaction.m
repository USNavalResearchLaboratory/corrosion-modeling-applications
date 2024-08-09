classdef ElectrochemicalOxidationReaction
    %ElectrochemicalOxidationReaction - Class for modeling oxidation
    %reactions
    %
    %The purpose of this class is to create an object that can be used to
    %model and manipulate simple anodic polarization data and models
    %
    %==========================================================================
    % Author:   Steve Policastro, Ph.D., Materials Science
    % Center for Corrosion Science and Engineering, U.S. Naval Research
    % Laboratory
    % email address: steven.policastro@nrl.navy.mil  
    % Website: 
    % October 2021; Last revision: 19-Oct-2021
    %==========================================================================

   properties (SetAccess = private)
       C = Constants;
   end   
   properties (SetAccess = public)
        name reactionNames
        Temperature double
        lambda_0 double   
        alpha double
        z double
        cReactants double
        cProducts double
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
        iLim double
        i double
        plotSymbol string
    end

    methods
        function obj = ElectrochemicalOxidationReaction(nameString, pS, T, c_react, c_prod, vApp, dG, alpha)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.name = nameString;
            obj.Temperature = T;
            obj.lambda_0 = (obj.C.kb*T)/obj.C.planck_h;           
            obj.dG_cathodic = dG(1);
            obj.dG_anodic = dG(2);
            obj.alpha = alpha;
            obj.plotSymbol = pS;
            obj.diffusionLength = -1.0;
            
            switch nameString
                case reactionNames.Cr_Ox
                    obj.z = obj.C.z_Cr_ox; 
                    c_Cr_g_cm3 = c_prod(1)*obj.C.M_Cr/1000; %g/cm3          
                    pre_factor = (obj.C.R*obj.Temperature)/(obj.z*obj.C.F); 
                    EN_log = log(c_react(1)/c_Cr_g_cm3);
                    obj.EN = obj.C.e0_Cr_ox + (pre_factor * EN_log);                    
                case reactionNames.Fe_Ox
                    obj.z = obj.C.z_Fe_ox; 
                    c_Fe_g_cm3 = c_prod(1)*obj.C.M_Fe/1000; %g/cm3          
                    pre_factor = (obj.C.R*obj.Temperature)/(obj.z*obj.C.F); 
                    EN_log = log(c_react(1)/c_Fe_g_cm3);
                    obj.EN = obj.C.e0_Fe_ox + (pre_factor * EN_log);                     
                case reactionNames.Ni_Ox
                    obj.z = obj.C.z_Ni_ox; 
                    c_Ni_g_cm3 = c_prod(1)*obj.C.M_Ni/1000; %g/cm3          
                    pre_factor = (obj.C.R*obj.Temperature)/(obj.z*obj.C.F); 
                    EN_log = log(c_react(1)/c_Ni_g_cm3);
                    obj.EN = obj.C.e0_Ni_ox + (pre_factor * EN_log);                    
                otherwise
            end
            obj.eta = vApp - (obj.EN - obj.C.E_SHE_to_SCE);    

            exp_val = -obj.dG_cathodic/(obj.C.R * obj.Temperature);
            exp_term = exp(exp_val);
            obj.i0_Cathodic = (obj.z*obj.C.F*obj.lambda_0) * exp_term;            

            exp_val = -obj.dG_anodic/(obj.C.R * obj.Temperature);
            exp_term = exp(exp_val);
            obj.i0_Anodic = (obj.z*obj.C.F*obj.lambda_0) * exp_term;
        
            obj.iCathodic = -obj.i0_Cathodic.*exp((-(1-obj.alpha)*obj.z*obj.C.F.*obj.eta)./(obj.C.R*obj.Temperature));
            obj.iAnodic = obj.i0_Anodic.*exp((obj.alpha*obj.z*obj.C.F.*obj.eta)./(obj.C.R*obj.Temperature));

            iAct = obj.iCathodic + obj.iAnodic;
            obj.i = iAct;              
        end
        
        function i = CalculateCurrent(obj,vApp)
            obj.eta = vApp - (obj.EN - obj.C.E_SHE_to_SCE);    

            exp_val = -obj.dG_cathodic/(obj.C.R * obj.Temperature);
            exp_term = exp(exp_val);
            obj.i0_Cathodic = (obj.z*obj.C.F*obj.lambda_0) * exp_term;            

            exp_val = -obj.dG_anodic/(obj.C.R * obj.Temperature);
            exp_term = exp(exp_val);
            obj.i0_Anodic = (obj.z*obj.C.F*obj.lambda_0) * exp_term;
        
            obj.iCathodic = -obj.i0_Cathodic.*exp((-(1-obj.alpha)*obj.z*obj.C.F.*obj.eta)./(obj.C.R*obj.Temperature));
            obj.iAnodic = obj.i0_Anodic.*exp((obj.alpha*obj.z*obj.C.F.*obj.eta)./(obj.C.R*obj.Temperature));

            iAct = obj.iCathodic + obj.iAnodic;
            i = iAct;               
        end

    end
    methods (Static)

        function newB = LM_Fit_Activation(objAnModel, yData, xData )
            
            maxIter = 500;
            tol = 1.0e-4;
            tau = 1.0e-3;
            B = [objAnModel.dG_anodic, objAnModel.alpha];
            k = length(B);
            NV = [1,1];
            Z = zeros(length(xData), length(B));
            lambda = 1.0e2;
            lambdaUp = 3.0;
            lambdaDown = 4.0;

            vMax = max(xData);
            vMin = min(xData);
            deltaV = (vMax - vMin)/length(xData);
            vRange = (vMax-deltaV):-deltaV:vMin;
            vRange = vRange';

            for iter = 1:maxIter
                testIterSum = 0;
                deltaB = 0.01.*B;                
                i0 = ActivationCurrentAnodic(objAnModel,vRange,B);
                ssd0 = CalculateSSD(i0,yData,k);

                for i=1:length(B)
                    testB = B;
                    testB(i) = B(i) + deltaB(i);
                    iH = ActivationCurrentAnodic(objAnModel,vRange,testB);
                    testB = B;
                    testB(i) = B(i) - deltaB(i);                
                    iL = ActivationCurrentAnodic(objAnModel,vRange,testB); 
                    Z(:,i) = (iH - iL)./(2*deltaB(i));
                end
    
                yMF = yData - abs(i0);
    
                ZprimeZ = Z'*Z;
                zDiag = diag(ZprimeZ);
                zD = diag(zDiag);

                lambdaOld = lambda;
                marquardtVal = lambdaOld.*zD;
                adjustedBold = linsolve((ZprimeZ + marquardtVal),Z'*yMF);
                B0 = B;
                Bold(:) = B0(:) + adjustedBold(:); 
                iOld = ActivationCurrentAnodic(objAnModel,vRange,Bold);
                ssdOld = CalculateSSD(iOld,yData,k);

                lambdaNew = lambda/lambdaDown;
                marquardtVal = lambdaNew.*zD;
                adjustedBNew = linsolve((ZprimeZ + marquardtVal),Z'*yMF);
                B0 = B;
                Bnew(:) = B0(:) + adjustedBNew(:); 
                iNew = ActivationCurrentAnodic(objAnModel,vRange,Bnew);                
                ssdNew = CalculateSSD(iNew,yData,k);

                if ssdNew < ssd0 && ssdOld < ssd0                    
                    if ssdNew < ssdOld
                        testConvergence = ssdNew/ssd0;
                        if testConvergence < tol
                        end
                        lambda = lambdaNew;
                        B = Bnew;
                        for i = 1:length(B)
                            testConvergence = abs(adjustedBNew(i))/(tau + abs(Bnew(i)));
                            if testConvergence < tol
                                NV(i) = 0;
                                testIterSum = testIterSum + 1;
                            end
                        end
                        if testIterSum == length(B)
                            break;
                        end                          
                    else
                        testConvergence = ssdOld/ssd0;
                        if testConvergence < tol

                        end
                        lambda = lambdaOld;
                        B = Bold;
                        for i = 1:length(B)
                            testConvergence = abs(adjustedBold(i))/(tau + abs(Bold(i)));
                            if testConvergence < tol
                                NV(i) = 0;
                                testIterSum = testIterSum + 1;
                            end
                        end
                        if testIterSum == length(B)
                            break;
                        end                          
                    end
                elseif ssdNew < ssd0 && ssdOld >= ssd0
                    lambda = lambdaNew;
                    B = Bnew;
                    for i = 1:length(B)
                        testConvergence = abs(adjustedBold(i))/(tau + abs(Bold(i)));
                        if testConvergence < tol
                            NV(i) = 0;
                            testIterSum = testIterSum + 1;
                        end
                    end
                    if testIterSum == length(B)
                        break;
                    end                     
                elseif ssdNew >= ssd0 && ssdOld < ssd0
                    lambda = lambdaOld;
                    B = Bold;
                    for i = 1:length(B)
                        testConvergence = abs((B(i) - B0(i))/B0(i));
                        if testConvergence < tol
                            NV(i) = 0;
                            testIterSum = testIterSum + 1;
                        end
                    end
                    if testIterSum == length(B)
                        break;
                    end                     
                else
                    lambda = lambda * lambdaUp;
                end

               
            end
            
%             figure(21)
%             hold on
%             plot(yData,xData,'ko')
%             plot(iOld,vRange,'-b')
%             ax = gca;
%             ax.XScale = 'log';
%             hold off

            newB = B;
%             newI = iOld;
%             disp(newB)

        end
        
    end
end