function exportFilenamesToExcel(folderPath, fileExtension, excelFileName)
% exportFilenamesToExcel - Scans a folder for files with a given extension
% and writes their names into an Excel workbook.
%
% Usage:
%   exportFilenamesToExcel('C:\MyFolder', '.txt', 'file_list.xlsx')

    % Validate inputs
    if ~isfolder(folderPath)
        error('The specified folder does not exist.');
    end

    if ~startsWith(fileExtension, '.')
        fileExtension = strcat('.', fileExtension);
    end

    % Get list of files with the given extension
    allfiles = dir(folderPath);
    listing = struct2table(allfiles); 
    disp(listing.name)
    onesWanted = endsWith(listing.name,fileExtension);
    fileList = listing.name(onesWanted);
    filenames = string(fileList);

    % If no files found, notify and exit
    if isempty(filenames)
        warning('No files with extension %s found in %s.', fileExtension, folderPath);
        return;
    end

    % Write to Excel
    try
        T = array2table(filenames); %,'VariableNames','filenames'
        writetable(T, excelFileName, 'Sheet', 1, 'Range', 'A1');
        fprintf('Successfully wrote %d filenames to %s\n', numel(filenames), excelFileName);
    catch ME
        error('Failed to write to Excel: %s', ME.message);
    end
end
