function sC = calcSC(Deff,EISExpt)
    termA = EISExpt.R*EISExpt.TemperatureImmersion;
    termB = sqrt(2) * ((EISExpt.zORR*EISExpt.F)^2) * EISExpt.Area;
    term1 = termA/termB;
    term2 = 1.0/(EISExpt.cO2 * sqrt(Deff));
    term3 = 1.0/(EISExpt.cOH *sqrt(EISExpt.DOH));
    sC = term1 * (term2 + term3);      
end