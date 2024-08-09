classdef ElectrochemicalOxidationReaction < handle
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
        % i0_oxgrowth double
        iCathodic double
        iAnodic double
        iLim double
        i double

    end

    methods        
        function obj = ElectrochemicalOxidationReaction(nameString, c_react, c_prod, T, vApp, metal)
            %ElectrochemicalOxidationReaction Construct an instance of this class
            %   Detailed explanation goes here
            obj.name = nameString;
            obj.Temperature = T;
            obj.lambda_0 = (obj.C.kb*T)/obj.C.planck_h;   

            switch nameString
                case reactionNames.Cr_Ox
                    obj.dG_cathodic = metal.DeltaGMetalOxidation(1);
                    obj.dG_anodic = metal.DeltaGMetalOxidation(2);
                    obj.alpha = metal.BetaMetalOxidation;                   
                case reactionNames.Fe_Ox
                    obj.dG_cathodic = metal.DeltaGMetalOxidation(1);
                    obj.dG_anodic = metal.DeltaGMetalOxidation(2);
                    obj.alpha = metal.BetaMetalOxidation;  
                case reactionNames.Cu_Ox
                    obj.dG_cathodic = metal.DeltaGMetalOxidation(1);
                    obj.dG_anodic = metal.DeltaGMetalOxidation(2);
                    obj.alpha = metal.BetaMetalOxidation;                      
                case reactionNames.Ni_Ox
                    obj.dG_cathodic = metal.DeltaGMetalOxidation(1);
                    obj.dG_anodic = metal.DeltaGMetalOxidation(2);
                    obj.alpha = metal.BetaMetalOxidation;                  
                case reactionNames.Passivation
                    obj.dG_cathodic = metal.DeltaGMetalPassivation(1);
                    obj.dG_anodic = metal.DeltaGMetalPassivation(2);
                    obj.alpha = metal.BetaMetalPassivation;   
                case reactionNames.Pitting
                    obj.dG_cathodic = metal.DeltaGMetalPitting(1);
                    obj.dG_anodic = metal.DeltaGMetalPitting(2);
                    obj.alpha = metal.BetaMetalPitting;
            end

            obj.diffusionLength = -1.0;
            obj.z = metal.OxidationLevelZ;
            pre_factor = (obj.C.R*obj.Temperature)/(obj.z*obj.C.F);

            switch nameString
                case reactionNames.Cr_Ox                     
                    c_g_cm3 = c_prod(1)*obj.C.M_Cr/1000; %g/cm3 
                    EN_log = log(c_react(1)/c_g_cm3);
                    obj.EN = obj.C.e0_Cr_ox + (pre_factor * EN_log);

                case reactionNames.Fe_Ox
                    c_g_cm3 = c_prod(1)*obj.C.M_Fe/1000; %g/cm3                              
                    EN_log = log(c_react(1)/c_g_cm3);
                    obj.EN = obj.C.e0_Fe_ox + (pre_factor * EN_log);                   

                case reactionNames.Cu_Ox
                    c_g_cm3 = c_prod(1)*obj.C.M_Cu/1000; %g/cm3          
                    EN_log = log(c_react(1)/c_g_cm3);
                    obj.EN = obj.C.e0_Cu_ox + (pre_factor * EN_log);     

                case reactionNames.Passivation
                    c_g_cm3 = c_prod(1)*obj.C.M_Cr/1000; %g/cm3          
                    EN_log = log(c_react(1)/c_g_cm3);
                    obj.EN = obj.C.e0_Cr_ox + (pre_factor * EN_log);   

                case reactionNames.Pitting
                    c_g_cm3 = c_prod(1)*obj.C.M_Cr/1000; %g/cm3
                    EN_log = log(c_react(1)/c_g_cm3);
                    obj.EN = obj.C.e0_Cr_ox + (pre_factor * EN_log);                    

                otherwise
                    error('No reaction name found!')
            end   

            eta0 = vApp - (obj.EN - obj.C.E_SHE_to_SCE);
            obj.eta = eta0;    

            exp_val = -obj.dG_cathodic/(obj.C.R * obj.Temperature);
            exp_term = exp(exp_val);
            obj.i0_Cathodic = (obj.z*obj.C.F*obj.lambda_0) * exp_term;            

            exp_val = -obj.dG_anodic/(obj.C.R * obj.Temperature);
            exp_term = exp(exp_val);
            obj.i0_Anodic = (obj.z*obj.C.F*obj.lambda_0) * exp_term;
            
            iActC = -obj.i0_Cathodic.*exp((-(1-obj.alpha)*obj.z*obj.C.F.*obj.eta)./(obj.C.R*obj.Temperature));
            iActA = obj.i0_Anodic.*exp((obj.alpha*obj.z*obj.C.F.*obj.eta)./(obj.C.R*obj.Temperature));

            if nameString == reactionNames.Passivation
   
                hOx = filmThickness(obj,metal,vApp);
                Resistance = metal.ResistivityOfOxide.*hOx;

                % figure(100)
                % hold on
                % plot(vApp,Resistance,'-b+')
                % xlabel('Time (s)')
                % ylabel('Resistance (\Omega)')                    
                % hold off

                iActCorrected = zeros(size(eta0));
                nReps = 50;                    
                tol = 1.0e-6;

                for j = 1:numel(vApp)                    
                    iOld = iActA(j);
                    Rp = Resistance(j);
                    eta00 = obj.eta(j);
                    for k = 1:nReps
                        if Rp > 0
                            % [feta,dfeta] = calculateneweta(obj,iOld,eta0(j),Rp); %etaOld
                            [fi,dfi] = calculatenewCurr(obj,iOld,obj.i0_Anodic,obj.alpha,eta00,obj.z,obj.C.F,obj.C.R,obj.Temperature,Rp);
                            iNew = iOld - (fi/dfi);
                        else
                            iNew = iOld;
                        end
                        err = abs((iNew-iOld)/iOld);
                        if err <= tol
                            iActCorrected(j) = iNew;
                            break;
                        elseif err > tol && k == nReps
                            % iAnCorrected(j) = iNew;
                            if iNew < iActCorrected(j-1)
                                iActCorrected(j) = 1.001*iNew;
                            else
                                iActCorrected(j) = iNew;
                            end
                        else
                            iOld = iNew;
                        end
                    end
                end

                % iActCorrected = iActA.*exp((obj.alpha*obj.z*obj.C.F.*iCorrected)./(obj.C.R*obj.Temperature)); 
                
                % % iSum = iCathodic + iActCorrected;
                % % obj.i = iSum;
                % etaErr = zeros(size(iAct));
                % for i = 1:numel(iAct)
                %     etaErr(i) = (eta0(i) - iCorrected(i))/eta0(i);
                % end

                % figure(200)
                % hold on
                % plot(abs(iActA),vApp,'-ko')
                % plot(abs(iActCorrected),vApp,'-b+')
                % % plot(abs(iActA),etaErr,'bo')
                % ax = gca;
                % ax.XScale = 'log';
                % xlabel('Current density (A/cm^2)')
                % % ylabel('Potential (V_{SCE})')
                % ylabel('Correction term') 
                % hold off   

                iActA = iActCorrected;
            end

            obj.iCathodic = iActC;
            obj.iAnodic = iActA;

            iSum = obj.iCathodic + obj.iAnodic;
            obj.i = iSum;   

        end
        
        function [fi,dfi] = calculatenewCurr(obj,iOld,i0,beta,eta,z,F,R,T,Rp)
            C1 = (beta*z*F)/(R*T);
            C2 = i0 * exp(C1*eta);
            Rcorrect = exp(-C1*Rp*iOld);
            fi = iOld - C2*Rcorrect;
            dfi = 1 + C2*C1*Rp*Rcorrect;
        end

        function hOx = filmThickness(obj,metal,vApp)
            %=================
            % Passive oxide layer formation....?
            %=================
            totalT = ((vApp(end) - vApp(1))*1000)/0.167;
            deltaT = totalT/(numel(vApp));                    
            tOx = 0.0:deltaT:((numel(vApp)-1)*deltaT);
                      
            aOx = obj.alpha; 
            eps2f = ((obj.C.R*obj.Temperature)/(aOx*obj.z*obj.C.F*metal.PassiveFilmThickness))*log(metal.PassiveCurrentDensity/obj.i0_Anodic);
            
            gOx = (aOx*obj.C.F)/(obj.C.R*obj.Temperature);

            kf = metal.OxideMass/(obj.z*obj.C.F*metal.OxideDensity); %cm3/C
            r = 1.0;   

            hOx = zeros(size(vApp));

            multi = 1.40e2; %14.0e10;
            eps2 = multi*eps2f; %2.2086e+07;
            A = eps2 * gOx;
            B = A*(kf/r)*metal.PassiveCurrentDensity*exp(gOx*eps2f*metal.PassiveFilmThickness); % metal.PassiveFilmThickness
        
            for j = 1:numel(vApp)   
                C1 = B*tOx(j);
                C2 = 1 + C1;
                aLayer = (1.0/A)*log(C2);
                hOx(j) = aLayer;
            end            
        end
    end

end