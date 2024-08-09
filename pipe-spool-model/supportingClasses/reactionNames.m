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
   enumeration
      HER (1)
      ORR (2)
      Fe_Ox (3)
      Fe_Red (4)
      Cr_Ox (5)
      Cr_Red (6)
      Ni_Ox (7)
      Ni_Red (8)
      None (9)
   end
end
 