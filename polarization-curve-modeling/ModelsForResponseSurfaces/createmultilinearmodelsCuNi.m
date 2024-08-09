function createmultilinearmodelsCuNi
    clc;
    clear all;
    
    ccl = [0.0 3.1 6.2 0.0 3.1 6.2 0.0 3.1 6.2];
    T = [5.0 5.0 5.0 25.0 25.0 25.0 45.0 45.0 45.0] + 273.15;
    pH = [1 7 13];

    metalName = 'cuni';

    %ORR
    % fn = strcat(metalName,'ORRCoeffs.csv');
    % dg_5 = [14.7 15.2 20.2].*1.0e4;
    % dg_25 = [13.0 15.2 17.7].*1.0e4;
    % dg_45 = [12.2 13.7 16.7].*1.0e4;
    
    %HER  
    % fn = strcat(metalName,'HERCoeffs.csv');
    % dg_5 = [9.7 11.2 15.2].*1.0e4;
    % dg_25 = [9.0 10.2 12.7].*1.0e4;
    % dg_45 = [7.2 8.7 11.7].*1.0e4;

    %CuOx
    fn = strcat(metalName,'CuOxCoeffs.csv');
    adjuster = 3.2;
    dg_5 = ([10.5 12.0 16.0]+adjuster).*1.0e4;
    dg_25 = ([9.8 11.0 13.5]+adjuster).*1.0e4;
    dg_45 = ([8.0 9.5 12.5]+adjuster).*1.0e4;

    dg_All = [dg_5, dg_25, dg_45];

    tbl = table(ccl',T',dg_All');
    
    sf = fit([tbl.Var1,tbl.Var2],tbl.Var3,'poly22');
    figure(2)
    plot(sf,[ccl',T'],dg_All')

    disp(sf)
    MyCoeffs = coeffvalues(sf);
    disp(MyCoeffs)
    writematrix(MyCoeffs,fn)
 
end