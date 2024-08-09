classdef SS316 < corrodingMetal & handle
    %SS316 Summary of this class goes here
    %   Detailed explanation goes here

    properties
    end

    methods

        function obj = SS316(name,ccl,T,pH)
            %corrodingMetal Construct an instance of this class
            %   Detailed explanation goes here
            obj.Name = name;
            obj.cCl = ccl;
            obj.T = T;
            obj.pH = pH;

            obj.MetalMass = 51.9961; %Molar mass of Cr, g/mol
    
            obj.OxidationLevelZ = 3; %Cr -> Cr3+ + 3e-
    
            % obj.DoesPit int32
            obj.PitPotential = -0.2; %VSCE
            obj.DeltaGMetalPitting = obj.CalculateDeltaG('Pitting',obj.cCl,obj.pH,obj.T);
            obj.BetaMetalPitting = 0.9999;
    
            obj.DeltaGMetalPassivation = obj.CalculateDeltaG('Passivation',obj.cCl,obj.pH,obj.T);
            obj.BetaMetalPassivation = 0.6; %0.21;
    
            obj.DeltaGORR = obj.CalculateDeltaG('ORR',obj.cCl,obj.pH,obj.T); 
            obj.BetaORR = 0.89;
            obj.delORR = 0.085; %cm 
    
            obj.DeltaGHER = obj.CalculateDeltaG('HER',obj.cCl,obj.pH,obj.T);
            obj.BetaHER = 0.8;
            obj.delHER = 0.15; %cm 
    
            obj.OxideMass = 151.99; %g/mol
            obj.OxideDensity = 5.22; %g/cm3 
            obj.ResistivityOfOxide = 5000.0e9; %Ohm/cm
            obj.PassiveCurrentDensity = 1.0e-3; %A/cm2
            obj.PassiveFilmThickness = 2.5e-7; %cm 
        end
   
        function actEnergies = CalculateDeltaG(obj,whichBarrier,ccl,ph,T)
            pHmax = 13;
            pHmin = 1;
            switch whichBarrier
                case 'ORR'
                    %Expression to predict dG_cathodic without pH
                    %dependence
                    data = readmatrix('SS316ORRCoeffs.csv');
                    p00 = data(1);
                    p10 = data(2);
                    p01 = data(3);
                    p20 = data(4);
                    p11 = data(5);
                    p02 = data(6);           
                    dgCathodic_nopH = p00 + p10*ccl + p01*T + p20*ccl^2 + p11*ccl*T + p02*T^2;
                    
                    dGCmax = 1.1*dgCathodic_nopH;
                    dGCmin = 0.9*dgCathodic_nopH;
                    m = (dGCmin-dGCmax)/(pHmax-pHmin);
                    dgCathodic = m*(ph-pHmax) + dGCmin;
                    
                    dgAnodic = 800.0e4;

                case 'HER'
                    %Expression to predict dG_cathodic without pH
                    %dependence
                    data = readmatrix('SS316HERCoeffs.csv');
                    p00 = data(1);
                    p10 = data(2);
                    p01 = data(3);
                    p20 = data(4);
                    p11 = data(5);
                    p02 = data(6);  
                    dgCathodic_nopH = p00 + p10*ccl + p01*T + p20*ccl^2 + p11*ccl*T + p02*T^2;
                    
                    dGCmax = 1.1*dgCathodic_nopH;
                    dGCmin = 0.9*dgCathodic_nopH;
                    m = (dGCmin-dGCmax)/(pHmax-pHmin);
                    dgCathodic = m*(ph-pHmax) + dGCmin;
                    
                    dgAnodic = 1000.0e4;              

                case 'Oxidation'
                    dgAnodic = 0.0;
                    dgCathodic = 0.0;

                case 'Passivation'
                    %Expression to predict dG_anodic without pH
                    %dependence
                    data = readmatrix('SS316PassCoeffs.csv');
                    p00 = data(1);
                    p10 = data(2);
                    p01 = data(3);
                    p20 = data(4);
                    p11 = data(5);
                    p02 = data(6);  
                    dgAnodic_nopH = p00 + p10*ccl + p01*T + p20*ccl^2 + p11*ccl*T + p02*T^2;
                    
                    dGAmax = 1.1*dgAnodic_nopH;
                    dGAmin = 0.9*dgAnodic_nopH;
                    m = (dGAmin-dGAmax)/(pHmax-pHmin);
                    dgAnodic = m*(ph-pHmax) + dGAmin;
                    
                    dgCathodic = 100.0e4; %15.0e4;    

                case 'Pitting'
                    %Expression to predict dG_cathodic without pH
                    %dependence
                    data = readmatrix('SS316PitCoeffs.csv');
                    p00 = data(1);
                    p10 = data(2);
                    p01 = data(3);
                    p20 = data(4);
                    p11 = data(5);
                    p02 = data(6);  
                    dgAnodic_nopH = p00 + p10*ccl + p01*T + p20*ccl^2 + p11*ccl*T + p02*T^2;
                    
                    dGAmax = 1.1*dgAnodic_nopH;
                    dGAmin = 0.9*dgAnodic_nopH;
                    m = (dGAmax-dGAmin)/(pHmax-pHmin);
                    dgAnodic = m*(ph-pHmin) + dGAmin;
                    
                    dgCathodic = 20.0e4;                        
            end

            actEnergies = [dgCathodic,dgAnodic]; 

        end
    
    end
end