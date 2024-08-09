function createmultilinearmodelsHY100
    clc;
    clear all;
    
    ccl = [0.0 3.1 6.2 0.0 3.1 6.2 0.0 3.1 6.2];
    T = [5.0 5.0 5.0 25.0 25.0 25.0 45.0 45.0 45.0] + 273.15;
    pH = [1 7 13];

    metalName = 'HY100';
    %ORR
    % fn = strcat(metalName,'ORRCoeffs.csv');
    % dg_5 = [13.6 11.1 19.1].*1.0e4;
    % dg_25 = [11.9 14.1 16.6].*1.0e4;
    % dg_45 = [11.1 12.6 15.6].*1.0e4;
    
    %HER  <===============
    % fn = strcat(metalName,'HERCoeffs.csv');
    % dg_5 = [9.6 11.1 15.1].*1.0e4;
    % dg_25 = [8.9 10.1 12.6].*1.0e4;
    % dg_45 = [7.2 8.6 11.6].*1.0e4;

    %FeOx
    % fn = strcat(metalName,'FeOxCoeffs.csv');
    % dg_5 = [14.6 16.1 20.1].*1.0e4;
    % dg_25 = [13.9 15.1 17.6].*1.0e4;
    % dg_45 = [12.1 13.6 16.6].*1.0e4;

    %Pitting
    fn = strcat(metalName,'PitCoeffs.csv');
    dg_5 = [34.1 20.1 18.1].*1.0e4;
    dg_25 = [21.6 19.1 17.9].*1.0e4;
    dg_45 = [20.6 17.5 16.1].*1.0e4;

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