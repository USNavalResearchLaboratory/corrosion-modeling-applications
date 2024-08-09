classdef HY100 < corrodingMetal & handle
    %HY100 Summary of this class goes here
    %   Detailed explanation goes here

    methods

        function obj = HY100(name,ccl,T,pH)
            %corrodingMetal Construct an instance of this class
            %   Detailed explanation goes here
            obj.Name = name;
            obj.cCl = ccl;
            obj.T = T;
            obj.pH = pH;

            obj.MetalMass = 55.845; %Molar mass of Fe, g/mol
    
            obj.OxidationLevelZ = 2; %Fe -> Fe2+ + 2e-
    
            % obj.DoesPit int32
            obj.PitPotential = -0.2; %VSCE
            obj.DeltaGMetalPitting = obj.CalculateDeltaG('Pitting',obj.cCl,obj.pH,obj.T);
            obj.BetaMetalPitting = 0.9999;
    
            obj.DeltaGMetalOxidation = obj.CalculateDeltaG('Oxidation',obj.cCl,obj.pH,obj.T);
            obj.BetaMetalOxidation = 0.3;
    
            obj.DeltaGORR = obj.CalculateDeltaG('ORR',obj.cCl,obj.pH,obj.T); 
            obj.BetaORR = 0.89;
            obj.delORR = 0.085; %cm 
    
            obj.DeltaGHER = obj.CalculateDeltaG('HER',obj.cCl,obj.pH,obj.T);
            obj.BetaHER = 0.72;
            obj.delHER = 0.15; %cm 
    
            obj.OxideMass = 151.99; %g/mol
            obj.OxideDensity = 5.22; %g/cm3 
            obj.ResistivityOfOxide = 5000.0e9; %Ohm/cm
            obj.PassiveCurrentDensity = 1.0e-6; %A/cm2
            obj.PassiveFilmThickness = 2.5e-7; %cm 
        end
   
        function actEnergies = CalculateDeltaG(obj,whichBarrier,ccl,ph,T)
            switch whichBarrier
                case 'ORR'
                    %Expression to predict dG_cathodic without pH
                    %dependence
                    data = readmatrix('HY100ORRCoeffs.csv');
                    p00 = data(1);
                    p10 = data(2);
                    p01 = data(3);
                    p20 = data(4);
                    p11 = data(5);
                    p02 = data(6);          
                    dgCathodic_nopH = p00 + p10*ccl + p01*T + p20*ccl^2 + p11*ccl*T + p02*T^2;
                    
                    dGCmax = 1.1*dgCathodic_nopH;
                    dGCmin = 0.9*dgCathodic_nopH;
                    m = (dGCmin-dGCmax)/(13-1);
                    dgCathodic = m*(ph-13) + dGCmin;
                    
                    dgAnodic = 800.0e4;

                case 'HER'
                    %Expression to predict dG_cathodic without pH
                    %dependence
                    data = readmatrix('HY100HERCoeffs.csv');
                    p00 = data(1);
                    p10 = data(2);
                    p01 = data(3);
                    p20 = data(4);
                    p11 = data(5);
                    p02 = data(6);    
                    dgCathodic_nopH = p00 + p10*ccl + p01*T + p20*ccl^2 + p11*ccl*T + p02*T^2;
                    
                    dGCmax = 1.1*dgCathodic_nopH;
                    dGCmin = 0.9*dgCathodic_nopH;
                    m = (dGCmin-dGCmax)/(13-1);
                    dgCathodic = m*(ph-13) + dGCmin;
                    
                    dgAnodic = 1000.0e4;              

                case 'Oxidation'
                    %Expression to predict dG_anodic without pH
                    %dependence
                    data = readmatrix('HY100FeOxCoeffs.csv');
                    p00 = data(1);
                    p10 = data(2);
                    p01 = data(3);
                    p20 = data(4);
                    p11 = data(5);
                    p02 = data(6);
                    dgAnodic_nopH = p00 + p10*ccl + p01*T + p20*ccl^2 + p11*ccl*T + p02*T^2;
                    
                    dGAmax = 1.1*dgAnodic_nopH;
                    dGAmin = 0.9*dgAnodic_nopH;
                    m = (dGAmin-dGAmax)/(13-1);
                    dgAnodic = m*(ph-13) + dGAmin;
                    
                    dgCathodic = 80.0e4;    

                case 'Pitting'
                    %Expression to predict dG_cathodic without pH
                    %dependence
                    data = readmatrix('HY100PitCoeffs.csv');
                    p00 = data(1);
                    p10 = data(2);
                    p01 = data(3);
                    p20 = data(4);
                    p11 = data(5);
                    p02 = data(6);
                    dgAnodic_nopH = p00 + p10*ccl + p01*T + p20*ccl^2 + p11*ccl*T + p02*T^2;
                    
                    dGAmax = 1.1*dgAnodic_nopH;
                    dGAmin = 0.9*dgAnodic_nopH;
                    m = (dGAmax-dGAmin)/(13-1);
                    dgAnodic = m*(ph-1) + dGAmin;
                    
                    dgCathodic = 20.0e4;                        
            end

            actEnergies = [dgCathodic,dgAnodic]; 

        end
    
    end
end