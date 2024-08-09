classdef butlerVolmer < handle
    %butlerVolmer Summary of this class goes here
    %   Detailed explanation goes here

    properties
        C Constants
        reaction char
        % alpha 
        % beta
        i0
        E0
        E0ORR
        E0HER
        E0Pass
        E0ox
        kappa
        T
        pH
        cCl
        velocity

        conc
        diff
        del

        aMetal
        aSoln        
    end

    methods
        function obj = butlerVolmer(r,a,b,i,e,k,t,z,c,D,del)
            %butlerVolmer Construct an instance of this class
            %   Detailed explanation goes here    
            obj.C = Constants();
            if nargin == 0
                obj.reaction ='null';
                obj.alpha = 0.0;
                obj.beta = 0.0;
                obj.i0 = 0.0;
                obj.E0 = 0.0; 
                obj.kappa = 0.0; 
                obj.T = 0.0;
                obj.z = 0; 
                obj.conc = 0.0;
                obj.diff = 0.0;
                obj.del = 1.0;

            elseif nargin == 1
                aReactionName = r;
                obj.reaction = r;
                solutionConductivity = 0.42; %4.2; %0.042; %S m^-1
                switch aReactionName
                    case 'ORR'
                    case 'HER'                        
                    case 'anodic'
                        obj.alpha = 0.05;
                        obj.beta = 0.05;
                        obj.i0 = 1.0; %1.0e-9; %
                        obj.E0 = -0.5; %0.5; %-1.5; %
                        obj.kappa = 10.0; 
                        obj.T = 25.0 + obj.C.convertCtoK;
                        obj.z = 1; 
                        obj.conc = 1.0;
                        obj.diff = 1.0;
                        obj.del = 1.0;
                    case 'cathodic'
                        obj.alpha = 0.05;
                        obj.beta = 0.05;
                        obj.i0 = 1.0; %1.0e-9; %
                        obj.E0 = 0.5; %1.5; %-0.5; %
                        obj.kappa = 10.0; 
                        obj.T = 25.0 + obj.C.convertCtoK;
                        obj.z = 1; 
                        obj.conc = 1.0e-5;
                        obj.diff = 1.0e-9;
                        obj.del = 850e-6;                                                    
                    case 'Zn'
                        obj.alpha = 0.025;
                        obj.beta = 0.05;
                        obj.i0 = 1.0; 
                        obj.E0 = -0.985; 
                        obj.kappa = solutionConductivity;
                        obj.T = 25.0 + obj.C.convertCtoK;
                        obj.z = 1; 
                        obj.conc = 1.0;
                        obj.diff = 1.0;
                        obj.del = 1.0;
                    case 'Cu'
                        obj.alpha = 0.001;
                        obj.beta = 0.05;
                        obj.i0 = 1.0; %1.0e-9; %
                        obj.E0 = -0.845; 
                        obj.kappa = solutionConductivity;
                        obj.T = 25.0 + obj.C.convertCtoK;
                        obj.z = 1; 
                        obj.conc = 1.0e1; %g/m^3
                        obj.diff = 1.0e-8; %m2/s
                        obj.del = 850e-6; %m                              
                end
            elseif nargin == 5 

                aReactionName = r;
                ccl = a;
                T = b;
                pH = i; 
                v = e;
                obj.reaction = r;
                obj.T = T + obj.C.convertCtoK;
                obj.pH = pH;
                obj.cCl = ccl;
                obj.velocity = v;
                obj.aSoln = naclSolutionChemistry(ccl,obj.T);
                
                resistivityScale = 1.0;

                switch aReactionName
                    case 'cuni'
                        obj.aMetal = CuNi(aReactionName,ccl,obj.T,obj.pH,obj.velocity);
                        obj.kappa = (obj.aSoln.rhoNaCl * resistivityScale); %Ohm/m  

                    case 'i625'
                        obj.aMetal = I625(aReactionName,ccl,obj.T,obj.pH,obj.velocity);
                        obj.kappa = (obj.aSoln.rhoNaCl * resistivityScale); %Ohm/m 
                end                
            else
                obj.reaction =r;
                obj.alpha = a;
                obj.beta = b;
                obj.i0 = i;
                obj.E0 = e; 
                obj.kappa = k; 
                obj.T = t + obj.C.convertCtoK;
                obj.z = z; 
                obj.conc = c;
                obj.diff = D;
                obj.del = del;                
            end
        end

        function [iRed,ia,ic] = multiCathodicI625(obj,eApp)

            %======================================
            % ORR
            %======================================
            dG_cathodic = obj.aMetal.DeltaGORR(1);
            dG_anodic = obj.aMetal.DeltaGORR(2);  
            alpha = obj.aMetal.BetaORR;
            diffusionLength = obj.aMetal.delORR;   

            lambda_0 = (obj.C.kb*obj.T)/obj.C.planck_h;

            pre_factor = (obj.C.R*obj.T)/(obj.C.z_orr*obj.C.F);
            [cOH, ~] = Constants.calculatecHandcOH(obj.pH);
            denom = (obj.aSoln.cO2/(obj.C.M_O2*1000))*(obj.aSoln.aW)^2;% V
            num = 1.0*(cOH^4); %/1000)*obj.obj.C.M_OH

            EN_log = log(num/denom); %log(denom/num); %
            obj.E0ORR = obj.C.e0_orr_alk + (pre_factor * EN_log);            
            
            eta = eApp - (obj.E0ORR - obj.C.E_SHE_to_SCE);
            R = obj.C.R;
            F = obj.C.F;
            RT = R*obj.T;
            dGC = dG_cathodic;
            dGA = dG_anodic;
            z = obj.C.z_orr;
            l0 = lambda_0;
            pF = z*F*l0;
            i0C = pF*exp(-dGC/RT);
            i0A = pF*exp(-dGA/RT);
            a = alpha;

            expValC = -(((1.0-a)*z*F).*eta)./RT;
            expValA = ((a*z*F).*eta)./RT;
            iCathodic = i0C.*exp(expValC);
            iAnodic = i0A.*exp(expValA); 
            
            iActORR = iAnodic - iCathodic;
            num = z*obj.C.F*obj.aSoln.dO2*(obj.aSoln.cO2/obj.C.M_O2);
            iLimORR = ones(size(eta)).*(-num/diffusionLength);

            iORR = (iLimORR.*iActORR)./(iActORR + iLimORR);
            %======================================
            % HER
            %======================================
            dG_cathodic = obj.aMetal.DeltaGHER(1);
            dG_anodic = obj.aMetal.DeltaGHER(2);  
            alpha = obj.aMetal.BetaHER;
            diffusionLength = obj.aMetal.delHER; %mm?   

            pre_factor = (obj.C.R*obj.T)/(obj.C.z_her*obj.C.F);

            [cOH, cH] = Constants.calculatecHandcOH(obj.pH);
            denom = 1.0 * ((cOH/1000)*obj.C.M_OH)^2;% V
            num = ((obj.aSoln.aW/1000)*obj.C.M_H2O)^2 * 1.0;

            EN_log = log(denom/num); %log(denom/num); %
            obj.E0HER = obj.C.e0_her_alk + (pre_factor * EN_log);            
            
            eta = eApp - (obj.E0HER - obj.C.E_SHE_to_SCE);            

            dGC = dG_cathodic;
            dGA = dG_anodic;
            z = obj.C.z_her;
            l0 = lambda_0;
            pF = z*F*l0;
            i0C = pF*exp(-dGC/RT);
            i0A = pF*exp(-dGA/RT);
            a = alpha;

            expValC = -(((1.0-a)*z*F).*eta)./RT;
            expValA = ((a*z*F).*eta)./RT;
            iCathodic = i0C.*exp(expValC);
            iAnodic = i0A.*exp(expValA); 
            
            iActHER = iAnodic - iCathodic;
            num = z*obj.C.F*obj.C.D_H2O*(obj.aSoln.aW/obj.C.M_H2O);
            iLimHER = ones(size(eta)).*(-num/diffusionLength);

            iHER = (iLimHER.*iActHER)./(iActHER + iLimHER);     

            %======================================
            % Oxidation
            %======================================
            c_prod = 1.0e-6;
            dG_cathodic = obj.aMetal.DeltaGMetalPassivation(1);
            dG_anodic = obj.aMetal.DeltaGMetalPassivation(2);
            alpha = obj.aMetal.BetaMetalPassivation;              

            c_g_cm3 = c_prod*obj.C.M_Cr/1000; %g/cm3 
            EN_log = log(1.0/c_g_cm3);
            obj.E0ox = obj.C.e0_Cr_ox + (pre_factor * EN_log); 
            z = obj.C.z_Cr_ox;

            eta = eApp - (obj.E0ox - obj.C.E_SHE_to_SCE);
            
            exp_val = -dG_cathodic/(obj.C.R * obj.T);
            exp_term = exp(exp_val);
            i0_Cathodic = (z*obj.C.F*lambda_0) * exp_term;            

            exp_val = -dG_anodic/(obj.C.R * obj.T);
            exp_term = exp(exp_val);
            i0_Anodic = (z*obj.C.F*lambda_0) * exp_term;
            
            iActC = -i0_Cathodic.*exp((-(1-alpha)*z*obj.C.F.*eta)./(obj.C.R*obj.T));
            iActA = i0_Anodic.*exp((alpha*z*obj.C.F.*eta)./(obj.C.R*obj.T));                    
            iOxidation = iActA + iActC;

            ia = 0.0;
            ic = 0.0;

            iRed = iORR + iHER + iOxidation;

            % figure(99)
            % hold on
            % plot(abs(iORR),eApp,':b')
            % plot(abs(iHER),eApp,':r')
            % plot(abs(iOxidation),eApp,':g')
            % plot(abs(iRed),eApp,'--k')
            % ax = gca;
            % ax.XScale = 'log';
            % hold off

            [imin, idxmin] = min(abs(iRed));
            obj.E0 = eApp(idxmin); 
            obj.i0 = abs(imin);
        end
    
        function [iSum,iaOx,icOx] = anodeKineticsCuNi(obj,eApp)

            %======================================
            % Oxidation
            %======================================            
            dG_cathodic = obj.aMetal.DeltaGMetalOxidation(1);
            dG_anodic = obj.aMetal.DeltaGMetalOxidation(2);  
            alpha = obj.aMetal.BetaMetalOxidation;

            lambda_0 = (obj.C.kb*obj.T)/obj.C.planck_h;

            pre_factor = (obj.C.R*obj.T)/(obj.C.z_orr*obj.C.F);
            
            denom = 1.0;% V
            num = 1.0e-6; %/1000)*obj.obj.C.M_OH

            EN_log = log(num/denom); %log(denom/num); %
            obj.E0ox = obj.C.e0_CuNi_ox + (pre_factor * EN_log); 
            eta =  eApp - obj.E0ox;
            
            R = obj.C.R;
            F = obj.C.F;
            RT = R*obj.T;
            dGC = dG_cathodic;
            dGA = dG_anodic;
            z = obj.aMetal.OxidationLevelZ;
            l0 = lambda_0;
            pF = z*F*l0;
            i0C = pF*exp(-dGC/RT);
            i0A = pF*exp(-dGA/RT);
            a = alpha;

            expValC = -(((1.0-a)*z*F).*eta)./RT;
            expValA = ((a*z*F).*eta)./RT;
            icOx = i0C.*exp(expValC);
            iaOx = i0A.*exp(expValA);    
            iOx = iaOx - icOx;

            %======================================
            % ORR
            %======================================
            dG_cathodic = obj.aMetal.DeltaGORR(1);
            dG_anodic = obj.aMetal.DeltaGORR(2);  
            alpha = obj.aMetal.BetaORR;
            diffusionLength = obj.aMetal.delORR;   

            lambda_0 = (obj.C.kb*obj.T)/obj.C.planck_h;

            pre_factor = (obj.C.R*obj.T)/(obj.C.z_orr*obj.C.F);
            [cOH, ~] = Constants.calculatecHandcOH(obj.pH);
            denom = (obj.aSoln.cO2/(obj.C.M_O2*1000))*(obj.aSoln.aW)^2;% V
            num = 1.0*(cOH^4); %/1000)*obj.obj.C.M_OH

            EN_log = log(num/denom); %log(denom/num); %
            obj.E0ORR = obj.C.e0_orr_alk + (pre_factor * EN_log);            
            
            eta = eApp - (obj.E0ORR - obj.C.E_SHE_to_SCE);
            R = obj.C.R;
            F = obj.C.F;
            RT = R*obj.T;
            dGC = dG_cathodic;
            dGA = dG_anodic;
            z = obj.C.z_orr;
            l0 = lambda_0;
            pF = z*F*l0;
            i0C = pF*exp(-dGC/RT);
            i0A = pF*exp(-dGA/RT);
            a = alpha;

            expValC = -(((1.0-a)*z*F).*eta)./RT;
            expValA = ((a*z*F).*eta)./RT;
            iCathodic = i0C.*exp(expValC);
            iAnodic = i0A.*exp(expValA); 
            
            iActORR = iAnodic - iCathodic;
            num = z*obj.C.F*obj.aSoln.dO2*(obj.aSoln.cO2/obj.C.M_O2);
            iLimORR = ones(size(eta)).*(-num/diffusionLength);

            iORR = (iLimORR.*iActORR)./(iActORR + iLimORR);
            %======================================
            % HER
            %======================================
            dG_cathodic = obj.aMetal.DeltaGHER(1);
            dG_anodic = obj.aMetal.DeltaGHER(2);  
            alpha = obj.aMetal.BetaHER;
            diffusionLength = obj.aMetal.delHER; %mm?   

            pre_factor = (obj.C.R*obj.T)/(obj.C.z_her*obj.C.F);

            [cOH, cH] = Constants.calculatecHandcOH(obj.pH);
            denom = 1.0 * ((cOH/1000)*obj.C.M_OH)^2;% V
            num = ((obj.aSoln.aW/1000)*obj.C.M_H2O)^2 * 1.0;

            EN_log = log(denom/num); %log(denom/num); %
            obj.E0HER = obj.C.e0_her_alk + (pre_factor * EN_log);            
            
            eta = eApp - (obj.E0HER - obj.C.E_SHE_to_SCE);            

            dGC = dG_cathodic;
            dGA = dG_anodic;
            z = obj.C.z_her;
            l0 = lambda_0;
            pF = z*F*l0;
            i0C = pF*exp(-dGC/RT);
            i0A = pF*exp(-dGA/RT);
            a = alpha;

            expValC = -(((1.0-a)*z*F).*eta)./RT;
            expValA = ((a*z*F).*eta)./RT;
            iCathodic = i0C.*exp(expValC);
            iAnodic = i0A.*exp(expValA); 
            
            iActHER = iAnodic - iCathodic;
            num = z*obj.C.F*obj.C.D_H2O*(obj.aSoln.aW/obj.C.M_H2O);
            iLimHER = ones(size(eta)).*(-num/diffusionLength);

            iHER = (iLimHER.*iActHER)./(iActHER + iLimHER);   

            iSum = iOx + iORR + iHER;
            % figure(99)
            % hold on
            % plot(abs(iaOx),eApp,':b')
            % plot(abs(icOx),eApp,':r')
            % plot(abs(iOx),eApp,'--k')
            % ax = gca;
            % ax.XScale = 'log';
            % hold off            
            
            [imin, idxmin] = min(abs(iSum));
            obj.E0 = eApp(idxmin);
            obj.i0 = abs(imin);            
        end    
        
    end

    methods (Static)
        function [fA,dfA] = onlyAnodic(EboundGuess,E1Mesh,h,bV)
            A = bV.kappa/(h*bV.i0);
            eta = E1Mesh - EboundGuess;            
            EB = bV.E0 + bV.a * log(A*eta);
            dEB = bV.a/eta;

            fA = EboundGuess - EB;
            dfA = 1 + dEB;            
        end
        
        function [fA,dfA] = onlyCathodic(EboundGuess,E1Mesh,h,bV)
            A = bV.kappa/(h*bV.i0);
            eta = EboundGuess - E1Mesh;
            EB = bV.E0 - bV.b * log(A*eta);
            dEB = -bV.b/eta;

            fA = EB - EboundGuess;
            dfA = dEB - 1;           
        end
        
        function [f,i,ia,ic,Eboundary] = anodePotential(Eguess,E1,h,bV,x)

            switch bV.reaction
                case 'cuni'
                    [i,ia,ic] = butlerVolmer.anodeKineticsCuNi(Eguess,bV);                    
                case default
                    [i,ia,ic] = butlerVolmer.anodeKineticsTafel(Eguess,bV);
            end
            
            p = pi;
            mult = 1.0e3;
            resistance = mult.*((1.0/bV.kappa)/((p*h)^2)).*abs(x);
            check = resistance < 1.0e0;
            resistance(check) = 1.0e0;
            conductivity = 1.0./resistance;

            Eboundary = zeros(size(x));
            Eguesses = ones(size(x)).*Eguess(1);

            switch bV.reaction
                case 'cuni'
                    for jj = 1:numel(Eboundary)
                        [current,~,~] = butlerVolmer.anodeKineticsCuNi(Eguesses(jj),bV);
                        Eboundary(jj) = E1 - ((h/conductivity(jj))*current); %bV.kappa
                    end                    
                case default
                    for jj = 1:numel(Eboundary)
                        [current,~,~] = butlerVolmer.anodeKineticsTafel(Eguesses(jj),bV);
                        Eboundary(jj) = E1 - ((h/conductivity(jj))*current); %bV.kappa
                    end                    
            end

            f = Eboundary - Eguesses;
        end

        function [f,i,ia,ic,Eboundary] = cathodePotential(Eguess,E1,h,bV,x)
            switch bV.reaction
                case 'cathodic'
                    [i,ia,ic] = butlerVolmer.cathodeKineticsTafel(Eguess,bV);
                case 'Cu'
                    [i,ia,ic] = butlerVolmer.fullCathodic(Eguess,bV);
                case 'SS138'
                    [i,ia,ic] = butlerVolmer.fullCathodic(Eguess,bV);
                case 'i625'
                    [i,ia,ic] = butlerVolmer.multiCathodicI625(Eguess,bV);                    
                case default
                    [i,ia,ic] = butlerVolmer.cathodeKineticsTafel(Eguess,bV);
            end
            p = pi;
            mult = 5.0e3;
            % multh = mult*h;              
            resistance = mult.*((1.0/bV.kappa)/((p*h)^2)).*abs(x);
            check = resistance < 1.0;
            resistance(check) = 1.0;
            conductivity = 1.0./resistance;

            Eboundary = zeros(size(x));
            Eguesses = ones(size(x)).*Eguess(1);
          
            switch bV.reaction
                case 'i625'
                    for jj = 1:numel(Eboundary)
                        [current,~,~] = butlerVolmer.multiCathodicI625(Eguesses(jj),bV);
                        Eboundary(jj) = E1 - ((h/conductivity(jj))*current); %bV.kappa
                    end                      
                case default
                    for jj = 1:numel(Eboundary)
                        [current,~,~] = butlerVolmer.fullCathodic(Eguesses(jj),bV);
                        Eboundary(jj) = E1 - ((h/conductivity(jj))*current); %bV.kappa
                    end                    
            end

            f = Eboundary - Eguesses;
        end        

        function [i,ia,ic] = anodeKineticsTafel(Eapp,bV)
            eta =  Eapp - bV.E0;
            
            eTermA = exp(eta./bV.alpha);
            ia = bV.i0 .* eTermA;

            eTermC = exp(-eta./bV.beta);
            ic = bV.i0 .* eTermC;      

            i = ia - ic;            
        end

        function [i,ia,ic] = cathodeKineticsTafel(Eapp,bV)
            eta =  Eapp - bV.E0;

            eTermC = exp(-eta./bV.beta);
            ic = bV.i0 .* eTermC;

            eTermA = exp(eta./bV.alpha);            
            ia = bV.i0 .* eTermA;

            i = -(ic - ia);         
        end        

        function [i,ia,iact] = fullCathodic(Eapp,bV)
            F = bV.C.F;            
            eta =  Eapp - bV.E0;

            eTermC = exp(-eta./bV.beta);
            iact = bV.i0 .* eTermC;
            iLim = ((bV.z*F*bV.diff*bV.conc)/bV.del);
            ic = (iact.*iLim)./(iact + iLim);

            eTermA = exp(eta./bV.alpha);            
            ia = bV.i0 .* eTermA;

            i = -(ic - ia);
        end   

    end
end
