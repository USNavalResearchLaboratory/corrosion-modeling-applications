function CreateExcelFileOfFilenamesForFitting
    clc;
    clear all;

    excelFileName = "files_dta1.xlsx";    
    folderPath = "C:\Users\steve\OneDrive\EIS Analysis\EIS Analyst MATLAB - Rev 1\Data\JHUAPL AM EIS 2\C58";
    fileExtension = "DTA";
    exportFilenamesToExcel(folderPath, fileExtension, excelFileName)
end