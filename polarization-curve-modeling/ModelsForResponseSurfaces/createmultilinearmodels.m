function createmultilinearmodels
    clc;
    clear all;
    
    ccl = [0.0 3.1 6.2 0.0 3.1 6.2 0.0 3.1 6.2];
    T = [5.0 5.0 5.0 25.0 25.0 25.0 45.0 45.0 45.0] + 273.15;
    pH = [1 7 13];

    metalName = 'SS316';

    %ORR
    % fn = strcat(metalName,'ORRCoeffs.csv');
    % adjuster = 5.5;
    % dg_5 = ([10.5 12.0 16.0] + adjuster).*1.0e4;
    % dg_25 = ([9.8 11.0 13.5] + adjuster).*1.0e4;
    % dg_45 = ([8.0 9.5 12.5] + adjuster).*1.0e4;
    
    %HER
    % fn = strcat(metalName,'HERCoeffs.csv');
    % adjuster = 5.0;
    % dg_5 = ([7.5 9.0 13.0] + adjuster).*1.0e4;
    % dg_25 = ([6.8 8.0 10.5] + adjuster).*1.0e4;
    % dg_45 = ([5.0 6.5 9.5] + adjuster).*1.0e4;

    %Passivation
    fn = strcat(metalName,'PassCoeffs.csv');
    adjuster = 9.0;
    dg_5 = ([12.5 14.0 18.0] + adjuster).*1.0e4;
    dg_25 = ([11.8 13.0 15.5] + adjuster).*1.0e4;
    dg_45 = ([10.0 11.5 14.5] + adjuster).*1.0e4;

    %Pitting
    % fn = strcat(metalName,'PitCoeffs.csv');
    % adjuster = 2.0;
    % dg_5 = ([52.0 48.0 46.5] + adjuster).*1.0e4;
    % dg_25 = ([49.5 47.0 45.8] + adjuster).*1.0e4;
    % dg_45 = ([48.5 45.5 44.0] + adjuster).*1.0e4;

    dg_All = [dg_5, dg_25, dg_45];

    tbl = table(ccl',T',dg_All');
    
    sf = fit([tbl.Var1,tbl.Var2],tbl.Var3,'poly22');

    figure(2)
    % hold on
    plot(sf,[ccl',T'],dg_All')
    xlabel('[Cl^-] (M)')
    ylabel('T (K)')
    zlabel('\DeltaG (J)')
    % hold off

    disp(sf)
    MyCoeffs = coeffvalues(sf);
    disp(MyCoeffs)
    writematrix(MyCoeffs,fn)    
 
end
