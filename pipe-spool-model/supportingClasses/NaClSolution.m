classdef NaClSolution
    %O2 contains functions for calculating dissolved oxygen properties
    %   This class contains properties for the dissolved oxygen
    %   concentration, cO2, in an NaCl solution and the diffusion
    %   coefficient, dO2, for the dissolved oxygen.  
    %   The setContants property imports the Constants class.
    %==========================================================================
    % Author:   Steven A. Policastro, Ph.D., Materials Science
    % Center for Corrosion Science and Engineering, U.S. Naval Research
    % Laboratory
    % email address: steven.policastro@nrl.navy.mil  
    % Website: 
    % October 2022; Last revision: 12-Oct-2022
    %==========================================================================
    properties
        setConstants Constants
        dO2 double
        cO2 double
        rhoNaCl double
        aW double
    end

    methods
        function obj = NaClSolution(cCl,T)
            %NaClSolution Constructs an instance of this class
            %   This contructor requires a temperature and Cl-
            %   concentration in (mol/L) as inputs so that values for the
            %   O2 diffusivity, dissolved O2 concentration, and solution
            %   conductivity can be calculated
            obj.setConstants = Constants;
            obj.cO2 = obj.calcConcO2(T,cCl); %mol/cm3
            obj.dO2 = obj.calcDiffO2(T,cCl); %cm2/s
            obj.rhoNaCl = 1.0./obj.calcSolnCond(T,cCl); % Ohm * m
            obj.aW = obj.waterActivity(cCl);
        end
    end
    methods (Access = private)
                
        function cO2 = calcConcO2(obj,T,cCl)
                %calcConcO2 Determines the dissolved oxygen concentration
                %in an NaCl solution.
                %   This function requires a temperature value (C) and Cl- 
                %   concentration (mol/L) to calculate the concentration of
                %   dissolved oxygen in the solution.         
                molecular_mass_Cl = obj.setConstants.M_Cl * obj.setConstants.convertGtoKg;
                if T >= obj.setConstants.convertCtoK
                    temperature_k = T;
                else
                    temperature_k = T + obj.setConstants.convertCtoK;
                end
                 
                Cl_molality_kg = molecular_mass_Cl .* cCl;
                Cl_molality_mg = Cl_molality_kg .* obj.setConstants.convertKgtoMg;
        
                a1 = 31820.0;
                b1 = -229.9;
                c1 = -19.12;
                d1 = 0.3081;
    
                a2 = -1409.0;
                b2 = 10.4;
                c2 = 0.8628;
                d2 = -0.0005235;
                
                d3 = 0.07464;
        
                acentric_factor_O2 = 0.022;
                num1 = (a1 * acentric_factor_O2) + a2;
                num2 = (b1 * acentric_factor_O2) + b2;
                denom1 = (c1 * acentric_factor_O2) + c2;
                denom2 = 1.0 + (denom1 .* temperature_k);
        
                Ln_H_s_0 = (num1 + (num2 .* temperature_k)) ./ denom2;
        
                num3 = d1 + (d2 .* temperature_k);
                denom3 = 1.0 + (d3 .* temperature_k);
        
                salinity_factor = 0.001;
                salinity = salinity_factor .* Cl_molality_mg;
                exp_term = (num3 ./ denom3) .* salinity;
        
                Ln_H_s = Ln_H_s_0 + exp_term;
        
                K_H = exp(Ln_H_s);
                x1 = obj.setConstants.cO2 ./ K_H;
                
                molecular_mass_O2 = obj.setConstants.M_O2 * obj.setConstants.convertGtoKg; %convert_g_to_kg;
                x1_g_L = x1 .* (molecular_mass_O2 / obj.setConstants.convertGtoKg); %convert_g_to_kg;
    
                x1_g_cm3 = x1_g_L./obj.setConstants.convertLtoCm3;
                cO2 = x1_g_cm3; %g/cm3
        end
        
        function DO2 = calcDiffO2(obj,T,cCl)
            %calcDiffO2 Determines the diffusivity of dissolved oxygen in
            %an NaCl solution.
            %   This function requires a temperature value (C) and Cl- 
            %   concentration (mol/L) to calculate the diffusivity of 
            %   dissolved oxygen in an NaCl solution.   
            lenT = length(T);
            lenC = length(cCl);
            if lenT >= lenC
                DO2 = zeros(size(T));
                N = lenT;
            else
                DO2 = zeros(size(cCl));
                N = lenC;
            end

            for j = 1:N
                if length(T) > 1
                    if T(j) >= obj.setConstants.convertCtoK
                        TK = T(j);
                    else
                        TK = T(j) + obj.setConstants.convertCtoK;
                    end                    
%                     TK = T(j)+ obj.setConstants.convertCtoK; 
                else
                    if T >= obj.setConstants.convertCtoK
                        TK = T(j);
                    else
                        TK = T(j) + obj.setConstants.convertCtoK;
                    end                      
%                     TK = T + obj.setConstants.convertCtoK;
                end
    
                params = [0.193015581, -0.000936823, -3738.145703; ...
                    0.586220598, -0.001982362, -0.003767555; ...
                    -2058331786, 7380780.538, -725742.0949; ...
                    -12341118, 7397.380585, -1024619.196; ...
                    -0.082481761, 8.05605E-06, -0.005230993; ...
                    -13685.50552, 11.9799009, -0.05822883];
            
                numModelParameters = size(params,1);
            
                b = zeros(numModelParameters,1);
                for i = 1:numModelParameters
                    b(i,1) = Constants.LinearLinear(params(i,:),TK);
                end
            
                DO2(j) = StokesModel2(obj,b,TK,cCl(j)); %cm2/s
            end
        end
        
        function D = StokesModel2(obj,b,TK,cCl)
            %StokesModel2 Determines the diffusivity of dissolved oxygen in
            %an NaCl solution using a Stokes model.
            %   This function requires a temperature value (C), Cl- 
            %   concentration (mol/L), and a set of b-parameters to 
            %   calculate the diffusivity of dissolved oxygen in an NaCl 
            %   solution.        
            phi = 2.6; 
        
            eta0 = b(5).*exp(b(6)./TK);
            B = b(3) + b(4).*(TK - Constants.convertCtoK); %
            A = b(2); 
            eta = eta0.*(1.0 + A.*sqrt(cCl) + B.*cCl);
        
            D = b(1).* ((sqrt(phi*obj.setConstants.M_H2O).*TK)./((obj.setConstants.VO2.*eta).^0.6));
        
        end

        function k = calcSolnCond(obj,TC,cCl)
        % solutionConductivity - Calculates the conductivity of a solution
        %
        % Function to calculate the conductivity of an NaCl solution. 
        % 
        % Wadsworth, J.C. 
        % The Statistical Description of Precision Conductivity Data for Aqueous Sodium Chloride. 
        % J Solution Chem 41, 715–729 (2012). 
        % https://doi.org/10.1007/s10953-012-9823-6
        %
        % Syntax:  k = solutionConductivity(T,c)
        %
        % Inputs: 
        % c = solution concentration (M) 
        % T = solution temperature (C)
        %  
        %
        % Outputs: Conductivity (S/m)
        %
        %
        % Other m-files required: none
        % Subfunctions: none
        % MAT-files required: none
        %
        % See also: 
        %
        %==========================================================================
        % Author:   Steve Policastro, Ph.D., Materials Science
        % Center for Corrosion Science and Engineering, U.S. Naval Research
        % Laboratory
        % email address: steven.policastro@nrl.navy.mil  
        % Website: 
        % June 2022; Last revision: 24 June 2022
        %==========================================================================
            if TC >= obj.setConstants.convertCtoK
                TC = TC - obj.setConstants.convertCtoK;
            end         
            b0 = -0.014; %uS/cm

            b10 = 66591.0;
            b11 = 2172.2;
            b12 = 9.1584;
            Lambda0 = b10 + b11.*TC + b12.*(TC.^2);

            b13 = 37515.0;
            b14 = -3471.9;
            b15 = 69.11;
            b16 = -1.0777;
            S = b13 + (b14.*TC) + (b15.*(TC.^2)) + b16.*(TC.^3);

            b17 = -23.47;
            E = b17.*(TC.^2);

            b18 = 46091;
            b19 = 8760;
            b20 = -352.06;
            b21 = 3.8403;
            J1 = b18 + (b19.*TC) + (b20.*(TC.^2)) + b21.*(TC.^3);
        
            b22 = -77300;
            b23 = -10646;
            b24 = 481.02;
            b25 = -4.9759;
            J2 = b22 + (b23.*TC) + (b24.*(TC.^2)) + b25.*(TC.^3);

            b26 = 98097;
            b27 = 5539.6;
            b28 = -242.12;
            b29 = 2.6452;
            J3 = b26 + (b27.*TC) + (b28.*(TC.^2)) + b29.*(TC.^3);  

            b30 = -68419;
            b31 = -1014.3;
            b32 = 43.97;
            b33 = -0.4871;
            J4 = b30 + (b31.*TC) + (b32.*(TC.^2)) + b33.*(TC.^3);  

            b34 = 22654;
            J5 = b34;

            b35 = -2799.6;
            J6 = b35;
            
            k1 = b0 + (Lambda0.*cCl) - (S.*(cCl.^(3/2))) + ...
                (E.*(cCl.^2).*log(cCl)) + (J1.*(cCl.^2)) + ...
                (J2.*(cCl.^(5/2))) + (J3.*(cCl.^(3))) + ...
                (J4.*(cCl.^(7/2))) + (J5.*(cCl.^(4))) + ...
                (J6.*(cCl.^(9/2))); %uS/cm

            k = k1.*(1.0e-6/0.01); %S/m
            
        end
        
        function aW = waterActivity(obj,clConc)
            % waterActivity - Calculates the water activity
            %
            % Function to calculate the activity of water in solution for a given NaCl
            % concentration
            %
            % Syntax:  aW = waterActivity(clConc)
            %
            % Inputs: 
            % clConc = solution concentration (M)
            %  
            %
            % Outputs: Water activity (M)
            %
            %
            % Other m-files required: none
            % Subfunctions: none
            % MAT-files required: none
            %
            % See also: 
            %
            %==========================================================================
            % Author:   Steve Policastro, Ph.D., Materials Science
            % Center for Corrosion Science and Engineering, U.S. Naval Research
            % Laboratory
            % email address: steven.policastro@nrl.navy.mil  
            % Website: 
            % June 2022; Last revision: 24 June 2022
            %==========================================================================
            mNaCl = obj.setConstants.M_NaCl * obj.setConstants.convertGtoKg; %kg/mol
            mH2O = obj.setConstants.M_H2O * obj.setConstants.convertGtoKg; %kg/mol
            molarity_of_water = 55.55; % mol/L
            test_solution_volume = 1.0; % L
        
            mass_NaCl_per_vol = mNaCl .* clConc; %kg/L
            mass_H2O_per_vol = mH2O .* obj.setConstants.cH2O; %kg/L
            total_mass_per_vol = mass_NaCl_per_vol + mass_H2O_per_vol;
        
            mass_percent_NaCl_in_sol = (mass_NaCl_per_vol ./ total_mass_per_vol) .* 100;
        
            d1 = 1.0001;% '0.0034287
            d2 = -0.0064603;% '0.99611
            density_NaCl_sol = d1 ./ (1.0 + (d2 .* mass_percent_NaCl_in_sol)) .* (obj.setConstants.convertGtoKg * obj.setConstants.convertLtoCm3); %g/mL - > kg/L
            mass_solution = density_NaCl_sol .* test_solution_volume; % kg, g/mL * 1000 mL/L
        
            mass_solvent = mass_solution - (mass_NaCl_per_vol .* test_solution_volume); %kg
            %mass_solvent = mass_solvent / 1000 'kg
            mol_Cl = clConc .* test_solution_volume;
        
            conc_cl_molality = mol_Cl ./ mass_solvent;
        
            c1 = 1.0001;
            c2 = -0.065634;
            c3 = -0.033533;
        
            %a_W = c1 + (c2 * conc_cl_molality) + (c3 * conc_cl_molality ^ 2)
            activity_coefficient = (c1 + (c2 .* conc_cl_molality)) ./ (1.0 + (c3 .* conc_cl_molality));
        
            aW = molarity_of_water .* activity_coefficient;    
        end    
    end
end