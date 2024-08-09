classdef aPipeLoop
    %aPipeLoop Summary of this class goes here
    %   Detailed explanation goes here
    properties(SetAccess = private)
        convFtToIn = 12.0;
        convInToCm = 2.54;
        convCmToM = 0.01;   
        pi = 3.14159265;
    end
    properties
        ID
        alloyType
        length {mustBeNumeric}
        innerDiameter {mustBeNumeric}
        areaExposure {mustBeNumeric}
        numberOfReferenceElectrodes {mustBeInteger}
        referenceElectrodeLocationsAll
        referenceElectrodeLocationsLL
        Places {mustBeNumeric}
        Potentials {mustBeNumeric}
        flowRate {mustBeNumeric}
    end

    methods
        function obj = aPipeLoop(i,material,L,innerD,locations,fR,refLL,pots)
            %aPipeLoop Construct an instance of this class
            %   Detailed explanation goes here
            obj.ID = i;
            obj.alloyType = material;
            obj.length = L*obj.convFtToIn*obj.convInToCm;
            obj.innerDiameter = innerD *obj.convInToCm;
            obj.areaExposure = obj.calculatePipeArea();
            [obj.numberOfReferenceElectrodes,obj.referenceElectrodeLocationsAll] = obj.setRefCellLocs(locations);
            obj.flowRate = fR;
            if numel(refLL) > 0
                obj.referenceElectrodeLocationsLL = refLL;
            end
            obj.Potentials = pots;
            obj.Potentials(obj.Potentials==0) = nan;
            obj.Places = zeros(size(pots));
            
            for i = 1:obj.numberOfReferenceElectrodes
                obj.Places(:,i) = obj.referenceElectrodeLocationsAll(i);
            end            
        end
    end
%     methods (Static)
%         function assignPotentials(aPipe,pots)
%             aPipe.Potentials = pots;
%             aPipe.Potentials(aPipe.Potentials==0) = nan;
%             aPipe.Places = zeros(size(pots));
%             for i = 1:aPipe.numberOfReferenceElectrodes
%                 aPipe.Places(:,i) = aPipe.referenceElectrodeLocationsAll(i);
%             end                
%         end        
%     end
    methods (Access = private)
        function area = calculatePipeArea(obj)
            %calculatePipeArea Summary of this method goes here
            %   Detailed explanation goes here
            r = obj.innerDiameter/2.0;            
            area = obj.length + (obj.pi * r^2);            
%             if strcmp(obj.alloyType, 'I625HX')
%                 area = 5760 * obj.convInToCm*obj.convInToCm;
%             else
%                 r = obj.innerDiameter/2.0;
%                 pi = 3.14159265;
%                 area = obj.length + (pi * r^2);
%             end

        end

        function [n,places] = setRefCellLocs(obj,locations)
            n = numel(locations);
            places = zeros(size(locations));
            for i = 1:n
                places(i) = locations(i)*obj.convInToCm *obj.convCmToM;
            end
        end
    end
end