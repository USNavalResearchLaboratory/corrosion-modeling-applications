function corpotTable = AnalyzeGamryCORPOTData(fn)
    %AnalyzeGamryCORPOTData  Extract data from CORPOT data files
    %   This function opens and extracts the relevant corrosion data and
    %   experimental parameters from a Gamry open circuit experiment data
    %   file.

    lyric2 = 'CURVE'; % Find the start of the columns of data
    
    [fid,msg] = fopen(fn,'rt');
    assert(fid>=3,msg)
    startheader = 0;
    startunits = 0;
    startdata = 0;
    dataline = 1;

    while ~feof(fid)
        st = fgetl(fid);           
        if contains(st,lyric2)
            stread = split(st);
            % disp(stread{1,2})
            numlines = str2double(stread{3,1});
            startheader = 1;        
        end
        if startheader == 1
            st = fgetl(fid);
            vars = char(split(st));
            startunits = 1;
            startheader = 0;
            datavals = zeros(numlines,size(vars,1));
        end
        if startunits == 1
            st = fgetl(fid);
            units = char(split(st)); % extracts the units, but we don't do anything with them...
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
    
    datavals2 = datavals(:,2:size(datavals,2));
    varnames = varnames(2:numel(varnames));
    corpotTable = array2table(datavals2,'VariableNames',varnames);
    
end
        