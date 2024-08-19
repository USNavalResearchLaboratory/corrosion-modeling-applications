function [area,cycpolTable] = AnalyzeGamryCYCPOLData_2(fn)
    %AnalyzeGamryCYCPOLData_2 - Extract data from CYCPOL data files
    %   This function opens and extracts the relevant corrosion data and
    %   experimental parameters from a Gamry cyclic polarization experiment
    %   data file.
    lyric1 = 'AREA'; % Find the recorded area value
    lyric2 = 'OCVCURVE'; % Find the start of the columns of OCP data
    lyric3 = 'CURVE'; % Find the start of the columns of CYCPOL data
    lyric4 = 'EXPERIMENTABORTED'; % Find out if the experiment ran to completion
    ocpdatavals = -1;
    
    [fid,msg] = fopen(fn,'rt');
    assert(fid>=3,msg)
    
    startocpheader = 0;
    startocpunits = 0;
    startocpdata = 0;

    startcycpolheader = 0;
    startcycpolunits = 0;
    startcycpoldata = 0;
    dataline = 1;

    while ~feof(fid)
        st = fgetl(fid);    
        % st = split(ast);

        % First, find the exposure area stored in the data file.  It may
        % not be correct, but let's grab it anyway...
        if contains(st,lyric1)
            stread = split(st);
            area = str2double(stread{3,1});    
        end 

        % Now, find the start of the OCP data, if it's present...
        foundocpcurve = 0;
        if contains(st,lyric2)
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

        % Now, find the start of the cycpol data
        if contains(st,lyric3)
            stread = split(st);
            numcycpollines = str2double(stread{3,1});
            startcycpolheader = 1;
            dataline = 1;
        end
        if startcycpolheader == 1
            st = fgetl(fid);
            vars = char(split(st));
            startcycpolunits = 1;
            startcycpolheader = 0;
            cycpoldatavals = zeros(numcycpollines,size(vars,1));
        end
        if startcycpolunits == 1
            st = fgetl(fid);
            units = char(split(st)); % extracts the units, but we don't do anything with them...
            startcycpoldata = 1;
            startcycpolunits = 0;
            st = fgetl(fid);
        end
        if startcycpoldata == 1
            if dataline <= numcycpollines && ~contains(st,lyric4)
                datas = str2double(split(st));
                cycpoldatavals(dataline,:) = datas;
                dataline = dataline + 1;
            end

        end
    end

    fclose(fid);

    if numel(ocpdatavals) > 1
        ocpvarnames = strtrim(string(ocpvars));
        datavals1 = ocpdatavals(:,2:size(ocpdatavals,2));
        ocpvarnames = ocpvarnames(2:numel(ocpvarnames));
        ocpTable = array2table(datavals1,'VariableNames',ocpvarnames);
    end

    cycpolvarnames = strtrim(string(vars));    
    datavals2 = cycpoldatavals(:,2:size(cycpoldatavals,2));
    cycpolvarnames = cycpolvarnames(2:numel(cycpolvarnames));
    cycpolTable = array2table(datavals2,'VariableNames',cycpolvarnames);            
end