function findeps2
    clc;
    clear all;

    eps0 = 1.0e-5;
    nReps = 50;
    tol = 1.0e-2; %0.4; %
    vals = zeros(1,nReps+1);
    error = zeros(size(vals));
    vals(1) = eps0;

    for j = 1:nReps
        [num,denom] = figuref(eps0); 
        frac = num/denom;
        eps1 = eps0 - frac;
        error(j) = abs((eps1 - eps0)/eps0);
        if error(j) <= tol
            disp(eps1)
            disp(j)
            break
        else
            vals(j+1) = eps1;
            eps0 = eps1;
        end        
    end

    nEps = vals > 0;
    n = 1:1:(nReps+1);
    mEps = vals < 0;

    figure(300)
    hold on
    plot(n(nEps),abs(vals(nEps)),'-bo')
    plot(n(mEps),abs(vals(mEps)),'-go')
    xlim([0 nReps+1])
    ax = gca;
    ax.YScale = 'log';
    hold off

    figure(301)
    hold on
    plot(n(nEps),abs(error(nEps)),'-r+')
    xlim([0 nReps+1])
    ax = gca;
    ax.YScale = 'log';
    hold off    
    
end
function [f,df] = figuref(eps)
    hf = 2.5e-7; %cm  
    tf = 4.4676e3; %s
    R  = 8.314;
    Temperature = 25 + 273.15;
    z = 3;
    F = 96485;

    i0growth = 5.3523e-7; %1.0e-6; %A/cm2    
    aOx = 0.5;
    eps2f = 4.317e4; 

    gOx = (aOx*F)/(R*Temperature);
    MCr2O3 = 151.99; %g/mol
    rhoCr2O3 = 5.22; %g/cm3    
    kf = MCr2O3/(z*F*rhoCr2O3); %cm3/C
    r = 1.0;     

    C1f = gOx*(kf/r)*i0growth*exp(gOx*eps2f*hf);
    f = ((exp(gOx*eps*hf))/tf) - C1f*eps;
    df = ((gOx*hf/tf)*exp(gOx*eps*hf)) - C1f;
end