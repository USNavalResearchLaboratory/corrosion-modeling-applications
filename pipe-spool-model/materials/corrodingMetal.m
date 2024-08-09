classdef corrodingMetal < handle
    %corrodingMetal Summary of this class goes here
    %   Detailed explanation goes here
   properties (SetAccess = private)
       C = Constants;
   end      
    properties
        %activationEnergiesRedox
        Name char
        MetalMass double

        OxidationLevelZ

        DoesPit int32
        PitPotential double
        DeltaGMetalPitting
        BetaMetalPitting double  

        DeltaGMetalPassivation 
        BetaMetalPassivation double

        DeltaGMetalOxidation 
        BetaMetalOxidation double

        DeltaGORR 
        BetaORR double
        delORR double

        DeltaGHER 
        BetaHER double      
        delHER double

        OxideMass double
        OxideDensity double
        ResistivityOfOxide double
        PassiveCurrentDensity double
        PassiveFilmThickness double

        cCl double
        T double
        pH double
    end

end