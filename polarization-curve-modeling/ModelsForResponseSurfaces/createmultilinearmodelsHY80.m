function createmultilinearmodelsHY80
    clc;
    clear all;
    
    ccl = [0.0 3.1 6.2 0.0 3.1 6.2 0.0 3.1 6.2];
    T = [5.0 5.0 5.0 25.0 25.0 25.0 45.0 45.0 45.0] + 273.15;
    pH = [1 7 13];

    metalName = 'HY80';

    %ORR
    fn = strcat(metalName,'ORRCoeffs.csv');
    dg_5 = [13.5 11.0 19.0].*1.0e4;
    dg_25 = [11.8 14.0 16.5].*1.0e4;
    dg_45 = [11.0 12.5 15.5].*1.0e4;
    
    %HER  
    % fn = strcat(metalName,'HERCoeffs.csv');
    % dg_5 = [9.5 11.0 15.0].*1.0e4;
    % dg_25 = [8.8 10.0 12.5].*1.0e4;
    % dg_45 = [7.0 8.5 11.5].*1.0e4;

    %FeOx
    % fn = strcat(metalName,'FeOxCoeffs.csv');
    % dg_5 = [14.5 16.0 20.0].*1.0e4;
    % dg_25 = [13.8 15.0 17.5].*1.0e4;
    % dg_45 = [12.0 13.5 16.5].*1.0e4;

    %Pitting
    % fn = strcat(metalName,'PitCoeffs.csv');
    % dg_5 = [34.0 20.0 18.5].*1.0e4;
    % dg_25 = [21.5 19.0 17.8].*1.0e4;
    % dg_45 = [20.5 17.5 16.0].*1.0e4;

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