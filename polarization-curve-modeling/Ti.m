classdef Ti < corrodingMetal & handle
    %Ti Summary of this class goes here
    %   Detailed explanation goes here


    methods

        function obj = Ti(name,ccl,T,pH)
            %Ti Construct an instance of this class
            %   Detailed explanation goes here
            obj.Name = name;
            obj.cCl = ccl;
            obj.T = T;
            obj.pH = pH;

            obj.MetalMass = 47.88; %Molar mass of Ti, g/mol
    
            obj.OxidationLevelZ = 3; %Ti -> Ti3+ + 3e-
   
            obj.DeltaGMetalPassivation = obj.CalculateDeltaG('Passivation',obj.cCl,obj.pH,obj.T);
            obj.BetaMetalPassivation = 0.3;
    
            obj.DeltaGORR = obj.CalculateDeltaG('ORR',obj.cCl,obj.pH,obj.T); 
            obj.BetaORR = 0.65;
            obj.delORR = 0.085; %cm 
    
            obj.DeltaGHER = obj.CalculateDeltaG('HER',obj.cCl,obj.pH,obj.T);
            obj.BetaHER = 0.75;
            obj.delHER = 0.15; %cm 
    
            obj.OxideMass = 143.76; %g/mol
            obj.OxideDensity = 4.49; %g/cm3 
            obj.ResistivityOfOxide = 50000.0e9; %Ohm/cm
            obj.PassiveCurrentDensity = 1.0e-6; %A/cm2
            obj.PassiveFilmThickness = 2.5e-7; %cm 
        end
   
        function actEnergies = CalculateDeltaG(obj,whichBarrier,ccl,ph,T)
            switch whichBarrier
                case 'ORR'
                    %Expression to predict dG_cathodic without pH
                    %dependence
                    data = readmatrix('TiORRCoeffs.csv');
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
                    data = readmatrix('TiHERCoeffs.csv');
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
                    dgAnodic = 0.0;
                    dgCathodic = 0.0;

                case 'Passivation'
                    %Expression to predict dG_anodic without pH
                    %dependence
                    data = readmatrix('TiPassCoeffs.csv');
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
                      
            end

            actEnergies = [dgCathodic,dgAnodic]; 

        end
    
    end
end