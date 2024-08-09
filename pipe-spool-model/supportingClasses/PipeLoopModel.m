%% Pipe Loop Model
% Steven A. Policastro, Ph.D. 
% Center for Corrosion Science and Engineering, 
% U.S. Naval Research Laboratory
% 4555 Overlook Avenue SW
% Washington, DC 20375 
% 
% 
% This class serves as a holding center for gathering experimental data
% and model outcomes from files.
%% Class definition
classdef PipeLoopModel < handle
    % PipeLoopModel - Class for handling experimental data and model
    % results.
%%% Public properties
    properties
        convertInToCm = 2.54;
        convertCmToM = 0.01;
        nPos
        xElsyca
        xSimple
        x2D
        y2D
        hasElsyca 
        hasMyModel 
        pipeLoopPotentialElsyca
        pipeLoopCurrentElsyca
        pipeLoopPotentialSimpleModel
        pipeLoopPotentialSimpleModel2
        pipeLoopPotentialSimpleModel3
    end

    methods
%%% Class Constructor 
% This mehtod receives a filenname that contains other model outcomes that
% can be loaded for comparison with the simple pipe spool model.
    function obj = PipeLoopModel(fn2)
            %PipeLoopModel - Class constructor.
            aCalc = readtable(fn2,'NumHeaderLines',1); 
            obj.nPos = size(aCalc,1);
            halfWay = round(obj.nPos/2);
            obj.xElsyca = aCalc.Var1- aCalc.Var1(halfWay); %  (obj.nPos:-1:1)
            obj.pipeLoopCurrentElsyca = aCalc.Var2./(100*100); %cm2
            obj.pipeLoopPotentialElsyca = aCalc.Var3; 
            obj.hasElsyca = false;
            obj.hasMyModel = false;
        end
    end
%%% Public static methods    
% This method, for determining the potential distribution from another
% model file and statistics on the dimensions of the pipes, do not require 
% an instantiated class object to exist
methods (Static)
        function [pot,dist] = GetData(~,aSim)            
            pot = aSim.phi(:,10); 
            totdist = aSim.L;
            npp = numel(pot);
            dx = totdist/(npp-1);
            halfDist = totdist/2;
            dist = -halfDist:dx:halfDist;            
        end
    end
end