%% Analyze Gamry POTDYN Data Function
% Steven A. Policastro, Ph.D. 
% Center for Corrosion Science and Engineering, 
% U.S. Naval Research Laboratory
% 4555 Overlook Avenue SW
% Washington, DC 20375 
% 
% 
% This function loads a Gamry DTA file obtained from a potentiodynamic
% polarization experiment. It extracts the open-circuit potential data and
% the polarization data from their respective areas of the file.  It will
% also extract exposed area value from AREA entry in the file.  The area
% value and the OCP and polarization data are returned from the function as
% a variable and a MATLAB table, respectively.
%
% Record of Software Revisions: 
% Created - June 2024
%
% Last revision: 20-July-2024
% 
% Contents
% 
%% Function definition
function [area,potdynTable] = AnalyzeGamryPOTDYNData(fn)
% AnalyzeGamryPOTDYNData - Extracts data from POTDYN data files
%%% Define Necessary Parameters
% These variables for analysis of the datafile  
    areaMarker = 'AREA'; % Find the recorded area value
    ocvMarker = 'OCVCURVE'; % Find the start of the columns of OCP data
    curveMarker = 'CURVE'; % Find the start of the columns of POTDYN data
    expAbortMarker = 'EXPERIMENTABORTED'; % Find out if the experiment ran to completion
    
    ocpdatavals = -1;
    startocpheader = 0;
    startocpunits = 0;
    startocpdata = 0;

    startcycpolheader = 0;
    startcycpolunits = 0;
    startcycpoldata = 0;
    dataline = 1;    

%%% File Opening
% Try to open filename and error if not successful  
    [fid,msg] = fopen(fn,'rt');
    assert(fid>=3,msg)
%%% File Scanning and Data Extraction
% Scan through the DTA file line-by-line.  Use the Marker variables to
% detect the different sections of the datafile so that we know when to
% extract certain lines.
    while ~feof(fid)
        st = fgetl(fid);    
        % First, find the exposure area stored in the data file.  It may
        % not be correct, but let's grab it anyway...
        if contains(st,areaMarker)
            stread = split(st);
            area = str2double(stread{3,1});    
        end 
        % Now, find the start of the OCP data, if it's present...
        foundocpcurve = 0;
        if contains(st,ocvMarker)
            foundocpcurve = 1;
            stread = split(st);
            numocplines = str2double(stread{3,1});
            startocpheader = 1;
        end
        if startocpheader == 1
            st = fgetl(fid);
            ocpvars = char(split(st));
            startocpunits = 1;
            startocpheader = 0;
            ocpdatavals = zeros(numocplines,size(ocpvars,1));
        end
        if startocpunits == 1
            st = fgetl(fid);
            units = char(split(st)); % extracts the units, but we don't do anything with them...
            startocpdata = 1;
            startocpunits = 0;
            st = fgetl(fid);
        end
        if startocpdata == 1
            if dataline >= numocplines
                startocpdata = 0;
            else
                datas = str2double(split(st));
                ocpdatavals(dataline,:) = datas;
                dataline = dataline + 1;
            end
        end
        % Now, find the start of the potdyn data
        if contains(st,curveMarker)
            stread = split(st);
            numpotdynlines = str2double(stread{3,1});
            startcycpolheader = 1;
            dataline = 1;
        end
        if startcycpolheader == 1
            st = fgetl(fid);
            vars = char(split(st));
            startcycpolunits = 1;
            startcycpolheader = 0;
            potdyndatavals = zeros(numpotdynlines,size(vars,1));
        end
        if startcycpolunits == 1
            st = fgetl(fid);
            units = char(split(st)); % extracts the units, but we don't do anything with them...
            startcycpoldata = 1;
            startcycpolunits = 0;
            st = fgetl(fid);
        end
        if startcycpoldata == 1
            if dataline <= numpotdynlines && ~contains(st,expAbortMarker)
                datas = str2double(split(st));
                potdyndatavals(dataline,:) = datas;
                dataline = dataline + 1;
            end

        end
    end

    fclose(fid);
%%% Create Output Tables
% Combine data vectors into the different data tables that will be returned
% by the function.
    if numel(ocpdatavals) > 1
        ocpvarnames = strtrim(string(ocpvars));
        datavals1 = ocpdatavals(:,2:size(ocpdatavals,2));
        ocpvarnames = ocpvarnames(2:numel(ocpvarnames));
        ocpTable = array2table(datavals1,'VariableNames',ocpvarnames);
    end

    potdynvarnames = strtrim(string(vars));    
    datavals2 = potdyndatavals(:,2:size(potdyndatavals,2));
    potdynvarnames = potdynvarnames(2:numel(potdynvarnames));
    potdynTable = array2table(datavals2,'VariableNames',potdynvarnames);            
end