classdef CuNi < corrodingMetal & handle
    %CuNi Summary of this class goes here
    %   Detailed explanation goes here

    methods

        function obj = CuNi(name,ccl,T,pH,liqV)
            %CuNi Construct an instance of this class
            %   Detailed explanation goes here
            obj.Name = name;
            obj.cCl = ccl;
            obj.T = T;
            obj.pH = pH;

            obj.MetalMass = 63.546; %Molar mass of Cu, g/mol
    
            obj.OxidationLevelZ = 1; %Cu -> Cu1+ + 1e-
    
            obj.DeltaGMetalOxidation = obj.CalculateDeltaG('Oxidation',obj.cCl,obj.pH,obj.T);
            obj.BetaMetalOxidation = 0.7; %0.4;
    
            obj.DeltaGORR = obj.CalculateDeltaG('ORR',obj.cCl,obj.pH,obj.T); 
            obj.BetaORR = .72; %0.65;
            liqV0 = 7.5; %m/s
            obj.delORR = 0.085*(1.0-(liqV/liqV0)); %cm 
    
            obj.DeltaGHER = obj.CalculateDeltaG('HER',obj.cCl,obj.pH,obj.T);
            obj.BetaHER = 0.6;
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
                    data = readmatrix('cuniORRCoeffs.csv');
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
                    data = readmatrix('cuniHERCoeffs.csv');
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
                    data = readmatrix('cuniCuOxCoeffs.csv');
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
                    
                    dgCathodic = 25.0e4; %13.2e4;                         
            end

            actEnergies = [dgCathodic,dgAnodic]; 

        end
    
    end
end