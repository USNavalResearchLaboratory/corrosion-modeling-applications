%% Constants and Conversion Factors
% Steven A. Policastro, Ph.D. 
% Center for Corrosion Science and Engineering, 
% U.S. Naval Research Laboratory
% 4555 Overlook Avenue SW
% Washington, DC 20375 
% 
%  
% This class contains variables storing physical constant values for use in
% electrochemical calculations.
% 
% _Record of Revisions_:
%
% * Created - October 2022
% * Last revision: 12-Oct-2022
classdef Constants
    % Defines physical constants and conversion factors for use in
    % electrochemical simulations.
    properties (Constant = true)
        R = 8.314; %Ideal gas constant, J/(mol K)
        F = 96485.3; %Faraday constant, coul/mol    
        kb = 1.38e-23;%Boltzmann's constant, m2 kg/(s2 K)
        planck_h = 6.626e-34; %Planck's constant, m2 kg / s
        eps0 = 8.85418782e-12; %Permittivity of free space m-3 kg-1 s4 A2
        e = 1.6e-19; % Charge on electron, C
        E_SHE_to_SCE = 0.244; %Convert standard hydrogen eelctrode potential to saturated calomel potential  
        convertCtoK = 273.15; %Convert celsius to Kelvin
        convertGtoKg = 1.0 / 1000.0; %Convert grams to kilograms
        convertKgtoMg = 1000.0 * 1000.0; %Convert kilograms to milligrams
        convertLtoCm3 = 1000; %Convert liters to cm3 (or mL)
        convertMtoCm = 100; %Convert meters to cm
        convertCm2toM2 = 1.0e-4; %Convert cm2 to m2

        cH2O = 55.55; %Concentration of water, mol/L
        cO2 = 0.209476; %Concentration of oxygen in air, %

        M_H2 = 2.016; %Molar mass of H2, g/mol
        M_OH = 17.008; %Molar mass of OH-, g/mol
        M_O2 = 32.0; %Molar mass of O2, g/mol
        M_H2O = 18.01528; %Molar mass of H2O, g/mol
        M_Cl = 35.5; %Molar mass of Cl-, g/mol
        M_NaCl = 58.4; %Molar mass of NaCl, g/mol
    
        M_Cr = 51.9961; %Molar mass of Cr, g/mol
        MCr2O3 = 151.99; %g/mol
        densityCr2O3 = 5.22; %g/cm3 
        rhoCr2O3_0 = 5000.0e9; %Ohm/cm
        iPassCr2O3 = 1.0e-6; %A/cm2
        tCr2O3film = 2.5e-7; %cm

        M_Fe = 55.845; %Molar mass of Fe, g/mol
        M_Ni = 58.6934; %Molar mass of Ni, g/mol
        M_Cu = 63.546;%Molar mass of Cu, g/mol
    
        D_H = 9.311e-5; %Diffusivity of H3O-, cm2/sec
        D_H2O = 2.299e-5; %Diffusivity of H2O, cm2/sec  
        D_Fe = 2.5e-12; %Diffusivity of Fe in oxide, cm2/sec  
		
		epsH2O = 80.1; %Dielectric constant of water
		epsPolyurethane = 6.19; %Dielectric constant of polyurethane
		epsEpoxy  = 3.6; %Dielectric constant of epoxy
		
        e0_orr_acid = 1.223; %Thermodynamic electrode potential of ORR in acidic solution, V_SHE
        e0_orr_2e_alk = -0.065;  %Thermodynamic electrode potential of ORR in acidic solution, V_SHE
        e0_orr_alk = 0.401;  %Thermodynamic electrode potential of ORR in neutral and alkaline solutions, V_SHE
        e0_her_alk = -0.83; %Thermodynamic electrode potential of HER in neutral and alkaline solutions, V_SHE
        e0_her_acid = 0.0; %Thermodynamic electrode potential of HER in acidic solution , V_SHE
        e0_me_ox = 0.0; %Thermodynamic electrode potential of generic metal, V_SHE
        e0_Cr_ox = -0.74; %Thermodynamic electrode potential of Cr oxidation, V_SHE
        e0_Fe_ox = -0.501; %-0.41;  %Thermodynamic electrode potential of Fe oxidation, V_SHE
        e0_Ni_ox = -0.23; %Thermodynamic electrode potential of Ni oxidation, V_SHE  
        e0_Cu_ox = 0.52; %Thermodynamic electrode potential of Ni oxidation, V_SHE
    
        z_orr = 4; %Number of electrons exchanged in the ORR reaction        
        z_her = 2; %Number of electrons exchanged in the HER reaction
        z_Cr_ox = 3; %Number of electrons exchanged in the Cr oxidation reaction
        z_Fe_ox = 2; %Number of electrons exchanged in the Fe oxidation reaction
        z_Cu_ox = 1; %Number of electrons exchanged in the Fe oxidation reaction
        z_Fe_red = 1; %Number of electrons exchanged in an Fe reduction reaction
        z_Ni_ox = 2; %Number of electrons exchanged in the Ni oxidation reaction

        VO2 = 22.414; %Molar volume of 1 mol of O2, L/mol
        VNaCl = 16.6; %Molar volume of 1 mol NaCl, L/mol
        % ==========================
    end

    methods
        function obj = Constants() 
            %Constants Construct an instance of this class
            %   Constructor of the constants class that takes no arguments
        end
    end

    methods (Static = true)
        function [cH,cOH] = calculatecHandcOH(pH)
            %calculatecHandcOH calculate cH and cOH
            %   Function that returns the concentration of H+ and OH-
            %   ions in solution, based on the input pH value.
            cH = 10.0^-(pH); %mol/L
            cOH = 10.0^-(14.0-pH); %mol/L
        end  

        function y = LinearLinear(b,x)
            %LinearLinear returns the solution of a linear-linear rational
            %   Function which, for a set of 3 parameters, b(1), b(2), and,
            %   b(3), and a domain, x, returns a range, y, for a
            %   linear-linear rational function.
            num = b(1) + b(2).*x;
            denom = 1.0 + b(3).*x;
            y = num./denom;
        end
        
    end

end