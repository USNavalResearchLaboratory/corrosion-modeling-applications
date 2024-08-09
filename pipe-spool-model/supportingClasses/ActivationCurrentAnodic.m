function i = ActivationCurrentAnodic(objAnModel, vApp, b) 
%ActivationCurrentAnodic - Calculate the oxidation reaction activation-controlled
%current density
%
% The purpose of this function is to calculate the activation-controlled
% current density for a given electrochemical oxidation reaction.
%
% Syntax:  i = ActivationCurrentAnodic(objCatModel, vApp, b)
%
% Inputs:
%    objAnModel - object of type ElectrochemicalOxidationReaction
%    vApp - Vector of applied potentials
%    b - 2-element vector containing values for the energy barrier for the
%    anodic reaction and the reaction symmetry parameter, beta
%
% Outputs:
%    i - Activation-controlled current density
%
%    Example
%    i0 = ActivationCurrentAnodic(objAnModel,vRange,B);
%
% Other m-files required: ElectrochemicalOxidationReaction
% Subfunctions: none
% MAT-files required: none
%
% See also: ActivationCurrentAnodic, TransportLimitedCurrent
%
%==========================================================================
% Author:   Steve Policastro, Ph.D., Materials Science
% Center for Corrosion Science and Engineering, U.S. Naval Research
% Laboratory
% email address: steven.policastro@nrl.navy.mil  
% Website: 
% October 2021; Last revision: 19-Oct-2021
%==========================================================================
    L1 = length(b);
    if L1 == 2
        exp_val = -b(1)/(objAnModel.C.R * objAnModel.Temperature);
        exp_term = exp(exp_val);
        i0_Anodic = (objAnModel.z*objAnModel.C.F*objAnModel.lambda_0) * exp_term; 
    
        an_eta = vApp - (objAnModel.EN - objAnModel.C.E_SHE_to_SCE);
        RT = objAnModel.C.R*objAnModel.Temperature;
        preFactor1 = b(2)*objAnModel.z*objAnModel.C.F;
        preFactor2 = preFactor1/RT;
        i = i0_Anodic.*exp(preFactor2.*an_eta);        
    else
       ME = MException('NotEnoughParameters','ActivationCurrentAnodic requires 2 parameters but received %s ',L1);
        throw(ME)         
    end

end