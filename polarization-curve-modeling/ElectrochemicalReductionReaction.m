classdef ElectrochemicalReductionReaction < handle
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
    end

    methods
        % reactionNames.ORR, cReact, cProd, Temp, C.z_orr, C.e0_orr_alk, Dcoeff, ...
        %                 obj.eAppModel, metal
        function obj = ElectrochemicalReductionReaction(nameString, c_ox, c_red, T, z, e0, D, vApp, metal)
            %ElectrochemicalReductionReaction Construct an instance of the ElectrochemicalReductionReaction class
            %The purpose of this method is to construct an instance of the
            %ElectrochemicalReductionReaction class.
            obj.name = nameString;
            switch obj.name
                case reactionNames.ORR
                    obj.dG_cathodic = metal.DeltaGORR(1);
                    obj.dG_anodic = metal.DeltaGORR(2);  
                    obj.alpha = metal.BetaORR;
                    obj.diffusionLength = metal.delORR; %mm?
                case reactionNames.HER
                    obj.dG_cathodic = metal.DeltaGHER(1);
                    obj.dG_anodic = metal.DeltaGHER(2);  
                    obj.alpha = metal.BetaHER;
                    obj.diffusionLength = metal.delHER; %mm?                    
            end
            obj.Temperature = T; % + Constants.convertCtoK
            obj.lambda_0 = (obj.c.kb*T)/obj.c.planck_h;           

            
            obj.cOxidized = c_ox;
            obj.cReduced = c_red;
            obj.E0 = e0;

            obj.z = z;
            pre_factor = (obj.c.R*obj.Temperature)/(obj.z*obj.c.F); 
            EN_log = log((c_red(1) * c_red(2))/(c_ox(1) * c_ox(2)));
            obj.EN = obj.E0 + (pre_factor * EN_log); %V_SHE     
            obj.diffusionCoefficient = D; %cm2/s ?

            if obj.name == reactionNames.None
                obj.iAct = zeros(size(vApp));
            else
                %Computes iAct, diAct/dDG, diAct/da, and diAct/dz
                obj.eta = vApp - (obj.EN - obj.c.E_SHE_to_SCE);
                R = obj.c.R;
                T = obj.Temperature;
                F = obj.c.F;
                RT = R*T;
                dGC = obj.dG_cathodic;
                dGA = obj.dG_anodic;
                obj.z = obj.z;
                l0 = obj.lambda_0;
                pF = obj.z*F*l0;
                i0C = pF*exp(-dGC/RT);
                i0A = pF*exp(-dGA/RT);
    
                obj.eta = obj.eta;
                obj.i0_Cathodic = i0C;            
                obj.i0_Anodic = i0A;
    
                a = obj.alpha;
                expValC = -(((1.0-a)*obj.z*F).*obj.eta)./RT;
                expValA = ((a*obj.z*F).*obj.eta)./RT;
    
                obj.iCathodic = i0C.*exp(expValC);
                obj.iAnodic = i0A.*exp(expValA); 
                obj.iAct = obj.iAnodic - obj.iCathodic;
            end            

            %Computes iLim
            switch obj.name
                case reactionNames.HER
                    num = obj.z*obj.c.F*obj.diffusionCoefficient*(obj.cOxidized(1)/obj.c.M_H2);
                case reactionNames.ORR
                    num = obj.z*obj.c.F*obj.diffusionCoefficient*(obj.cOxidized(1)/obj.c.M_O2);
                case reactionNames.Fe_Red
                    num = obj.z*obj.c.F*obj.diffusionCoefficient*(obj.cOxidized(1)/obj.c.M_Fe);
                case reactionNames.None
                    num = 0.0;
                otherwise
                    error('No reaction name found!')
            end

            obj.iLim = ones(size(obj.eta)).*(-num/obj.diffusionLength); %A/cm2    

            if obj.name == reactionNames.None
                obj.i = zeros(size(obj.eta));
            else
                obj.i = (obj.iLim  .* obj.iAct)./(obj.iAct + obj.iLim ); 
            end    

            % figure(400)
            % hold on
            % plot(abs(obj.i), vApp,'-b+')
            % ax = gca;
            % ax.XScale = 'log';
            % hold off

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