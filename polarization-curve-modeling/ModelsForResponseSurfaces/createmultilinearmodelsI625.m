function createmultilinearmodelsI625
    clc;
    clear all;
    
    ccl = [0.0 3.1 6.2 0.0 3.1 6.2 0.0 3.1 6.2];
    T = [5.0 5.0 5.0 25.0 25.0 25.0 45.0 45.0 45.0] + 273.15;
    pH = [1 7 13];

    metalName = 'I625';
    %ORR
    % fn = strcat(metalName,'ORRCoeffs.csv');
    % dg_5 = [14.2 15.7 19.7].*1.0e4;
    % dg_25 = [13.5 14.7 17.2].*1.0e4;
    % dg_45 = [11.7 13.2 16.2].*1.0e4;
    
    %HER
    % fn = strcat(metalName,'HERCoeffs.csv');
    % dg_5 = [11.1 12.6 16.6].*1.0e4;
    % dg_25 = [10.4 11.6 14.1].*1.0e4;
    % dg_45 = [8.6 10.1 13.1].*1.0e4;

    %Passivation
    fn = strcat(metalName,'PassCoeffs.csv');
    adjuster = 4.0;    
    dg_5 = ([16.7 18.2 22.2]+adjuster).*1.0e4;
    dg_25 = ([16.0 17.2 19.7]+adjuster).*1.0e4;
    dg_45 = ([14.2 15.7 18.7]+adjuster).*1.0e4;

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