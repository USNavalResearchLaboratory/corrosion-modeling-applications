function [area,eisTable] = AnalyzeGamryEISData(fn)
    %AnalyzeGamryEISData Extract data from EISPOT data files
    %   This function opens and extracts the relevant corrosion data and
    %   experimental parameters from a Gamry potentiostatic electrochemical
    %   impedance spectroscopy (EIS) data file.
    
    lyric1 = 'AREA'; % Find the recorded area value
    lyric2 = 'ZCURVE'; % Find the start of the columns of data
    
    [fid,msg] = fopen(fn,'rt');
    assert(fid>=3,msg)
    startheader = 0;
    startunits = 0;
    startdata = 0;
    dataline = 1;

    while ~feof(fid)
        st = fgetl(fid);    
        if contains(st,lyric1)
            stread = split(st);
            % disp(stread{1,2})
            area = str2double(stread{3,1});    
        end        
        if contains(st,lyric2)
            startheader = 1;        
        end
        if startheader == 1
            st = fgetl(fid);
            vars = char(split(st));
            startunits = 1;
            startheader = 0;
            datavals = ones(500,size(vars,1)).*-1;
        end
        if startunits == 1
            st = fgetl(fid);
            units = char(split(st));
            startdata = 1;
            startunits = 0;
            st = fgetl(fid);
        end
        if startdata == 1
            datas = str2double(split(st));
            datavals(dataline,:) = datas;
            dataline = dataline + 1;
        end
    end
    fclose(fid);
    varnames = strtrim(string(vars));
    
    dv = datavals(:,2) > -1;
    datavals2 = datavals(dv,2:size(datavals,2));
    varnames = varnames(2:numel(varnames));
    eisTable = array2table(datavals2,'VariableNames',varnames);            
end
            