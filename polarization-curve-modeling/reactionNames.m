classdef reactionNames < uint8
    %reactionNames defines an enumeration to contain the available
    %reactions that can be modeled
    %   Enumerates a list of reaction names used by the
    %   electrochemicalReaction class to determine which reaction to model
    %==========================================================================
    % Author:   Steven A. Policastro, Ph.D., Materials Science
    % Center for Corrosion Science and Engineering, U.S. Naval Research
    % Laboratory
    % email address: steven.policastro@nrl.navy.mil  
    % Website: 
    % October 2022; Last revision: 12-Oct-2022
    %==========================================================================
    % =====
    % ORR - 4e- Alkaline
    % O2 + 2H2O + 4e- -> 4OH- 
    % =====
    % ORR - 4e- acid
    % =====
    % O2 + 4H+ + 4e- -> 2H2O    
    % =====
    % ORR - 2e- alkaline
    % =====
    % O2 + H2O + 2e- -> HO2- + OH-
    % HO2- + H2O + 2d- -> 3OH-
    % =====       
    % =====    
    % HER
    % 2H+ + 2e- -> H2
    % 2H2O + 2e- -> H2 + 2OH-
    % =====      
   enumeration
      HER (1)
      ORR (2)
      Fe_Ox (3)
      Fe_Red (4)
      Cr_Ox (5)
      Cr_Red (6)
      Ni_Ox (7)
      Ni_Red (8)
      Cu_Ox (9)
      Cu_Red (10)
      Al_Ox (11)
      Al_Red (12)
      Passivation (13)
      Scale (14)
      Pitting (15)
      None (16)
   end
end
 