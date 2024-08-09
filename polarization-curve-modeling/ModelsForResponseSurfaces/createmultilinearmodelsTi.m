function createmultilinearmodelsTi
    clc;
    clear all;
    
    ccl = [0.0 3.1 6.2 0.0 3.1 6.2 0.0 3.1 6.2];
    T = [5.0 5.0 5.0 25.0 25.0 25.0 45.0 45.0 45.0] + 273.15;
    pH = [1 7 13];

    metalName = 'Ti';
    %ORR
    % fn = strcat(metalName,'ORRCoeffs.csv');
    % dg_5 = [15.2 16.7 20.7].*1.0e4;
    % dg_25 = [14.5 15.7 18.2].*1.0e4;
    % dg_45 = [12.7 14.2 17.2].*1.0e4;
    
    %HER
    % fn = strcat(metalName,'HERCoeffs.csv');
    % dg_5 = [12.1 13.6 17.6].*1.0e4;
    % dg_25 = [11.4 12.6 15.1].*1.0e4;
    % dg_45 = [9.6 12.1 14.1].*1.0e4;

    %Passivation
    fn = strcat(metalName,'PassCoeffs.csv');
    dg_5 = [15.7 17.2 21.2].*1.0e4;
    dg_25 = [15.0 16.2 18.7].*1.0e4;
    dg_45 = [13.2 14.7 17.7].*1.0e4;

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