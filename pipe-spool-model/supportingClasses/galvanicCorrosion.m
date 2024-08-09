%% Galvanic Corrosion Simulation
% Steven A. Policastro, Ph.D. 
% Center for Corrosion Science and Engineering, 
% U.S. Naval Research Laboratory
% 4555 Overlook Avenue SW
% Washington, DC 20375 
% 
% 
% This class defines the computational cell and contains the static
% functions that perform the Symmetric Over-Relaxation (SOR) calculation to
% obtain the potential distribution by solving the Laplace equation in 2
% dimensions.
%
%% Class definition
classdef galvanicCorrosion
    %galvanicCorrosion Calculates the galvanic interaction between 2 metals
    %   This class contains properties and functions to define material
    %   properties, electrochemical reactions, and a solver for the laplace
    %   equation
%%% Public properties
% Contains only a single property of type galvCorrSim that calculates the
% polarization curves for the different materials in the galvanic couple.
    properties
        aSim galvCorrSim
    end

    methods
%%% Class Constructor 
% This method receives multiple parameters that are used to define the
% computational cell and calculatet the polarization curves using the
% methods developed in the Polarization Curve Model project.
        function obj = galvanicCorrosion(topBC,length,height,dL,eNodes,numNodes,aNodes,cNodes,aReactName,cReactName,Vapp,nF,env)
            %galvanicCorrosion Construct an instance of this class           
            obj.aSim = galvCorrSim(topBC,length,height,dL,eNodes,numNodes,aNodes,cNodes,aReactName,cReactName,Vapp,nF,env); %
            galvCorrSim.aTafelPolCurve(obj.aSim);
        end
    end
%%% Public static methods    
% These methods peform the following tasks:
%
methods (Static)
%%%
%
% * JacobiSolver - Uses the Jacbi method to perform the SOR algorithm to 
% obtain the potential distribution in 2D.
        function aPhi = JacobiSolver(aSim)
            %JacobiSolver - Uses the Jacbi method to perform the SOR
            %algorithm to obtain the potential distribution in 2D
            maxIter = aSim.NY*aSim.NX; 
            maxIter2 = 15;
            tol = 1.0e-6; 
            converged2 = false;

            N = aSim.NX*aSim.NY;
            E = zeros(1,N);
            
            L = 1;
            aPhi = zeros(aSim.NY,aSim.NX);

            idx = 1:1:N;   

            topPts = zeros(1,N); 
            leftPts = zeros(1,N); 
            rightPts = zeros(1,N); 
            botPts = zeros(1,N); 
            rowAboveBotPts = zeros(1,N); 

            topPts(1:aSim.NX) = idx(1:aSim.NX);
            botPts((N-aSim.NX)+1:N) = idx((N-aSim.NX)+1:N);
            rowAboveBotPts(1:aSim.NX) = idx((N-aSim.NX)+1:N); %idx((N-(2*aSim.NX))+1:(N-aSim.NX));
            leftPts(aSim.NX+1:aSim.NX:(N-(2*aSim.NX))+1) = idx(aSim.NX+1:aSim.NX:(N-(2*aSim.NX))+1);
            rightPts((2*aSim.NX):aSim.NX:N-aSim.NX) = idx((2*aSim.NX):aSim.NX:N-aSim.NX);

            intPts1 = ~(idx&topPts); 
            intPts2 = ~(idx&botPts); 
            intPts3 = ~(idx&leftPts); 
            intPts4 = ~(idx&rightPts);
            intPtsL = intPts1&intPts2&intPts3&intPts4;
            intPts = idx(intPtsL);

            topPts = nonzeros(topPts);
            botPts = nonzeros(botPts);
            leftPts = nonzeros(leftPts);
            rightPts = nonzeros(rightPts);
            rowAboveBotPts = nonzeros(rowAboveBotPts);

            % Bottom electrode
            NB = numel(botPts);
            NT = numel(topPts);
            % halfway = floor(NB/2);
            NBC = 1:aSim.NXc;
            NBA = (aSim.NXc+1):(aSim.NXa+aSim.NXc);
            EBL = zeros(size(NBC));
            EBR = zeros(size(NBA));

            Elp = zeros(size(leftPts));
            Erp = zeros(size(rightPts));
            Etop = zeros(size(topPts));
            Eint = zeros(size(intPts));

            % ===============
            % Base
            % =============== 
            EbotsAboveC(NBC) = aSim.eCorr;
            EbotsAboveA(NBA-101) = aSim.eCorr;            
            EBL(NBC) = galvanicCorrosion.newtonRaphson2(EbotsAboveC,'cathodic',aSim, aSim.dy, aSim.xpos);          
            EBR(NBA) = galvanicCorrosion.newtonRaphson2(EbotsAboveA,'anodic',aSim, aSim.dy, aSim.xpos);            

            while L <= maxIter2 && converged2 == false
                converged3 = false;
                converged1 = false;
                k = 1;

                while k <= maxIter && converged1 == false  

                    if converged3 == false
                        Ebot = E(botPts);
                        for ii = 1:10
                            Ebot(1) = (EBL(1) + E(botPts(1)-aSim.NX) + 2*Ebot(2))/4.0;                    
                            for i = 2:numel(Ebot)-1
                                Eleft = Ebot(i-1);
                                Eright = Ebot(i+1);
                                if (i <= NBC(end))
                                    Ebottom = EBL(i);
                                else
                                    Ebottom = EBR(i);
                                end
                                Etopper = E(botPts(i)-aSim.NX);                                                
                                Ebot(i) = (Eleft + Eright + Ebottom + Etopper)./4.0;
                            end
    
                            Ebot(end) = (EBR(end) + E(botPts(end)-aSim.NX) + 2*Ebot(numel(Ebot)-1))/4.0;
                        end

                        E(botPts) = Ebot;
                        converged3 = true;
                        % figure(101)
                        % hold on
                        % plot(Ebot,'-b')
                        % hold off
                        % disp('Here!')
                    end
    
                    % ===============
                    % Sides
                    % ===============                        
                    % Left                   
                    switch aSim.bcL
                        case 'dirichlet'
                            if numel(aSim.Eapp) > 0
                                Elp = aSim.Eapp(2);
                            else
                                Elp = 0.0;
                            end
                        case 'neumann'
                            for i = 1:numel(leftPts)
                                ptNum = leftPts(i);
                                Eright = E(ptNum + 1);
                                Etopper = E(ptNum-aSim.NX);
                                Ebottom = E(ptNum+aSim.NX);
                                Elp(i) = ((2.0*Eright) + Etopper + Ebottom)./4.0;
                            end
                    end
                    E(leftPts) = Elp;
    
                    % Right                 
                    switch aSim.bcR
                        case 'dirichlet'
                            if numel(aSim.Eapp) > 0
                                Erp = aSim.Eapp(3);
                            else
                                Erp = 0.0;
                            end                        
                        case 'neumann'
                            for i = 1:numel(rightPts)
                                ptNum = rightPts(i);
                                Eleft = E(ptNum - 1);
                                Etopper = E(ptNum-aSim.NX);
                                Ebottom = E(ptNum+aSim.NX);                            
                                Erp(i) = ((2.0*Eleft) + Etopper + Ebottom)./4.0;
                            end
                    end
                    E(rightPts) = Erp;
    
                    % ===============
                    % Top
                    % ===============  
                    switch aSim.bcT
                        case 'dirichlet'
                            if numel(aSim.Eapp) > 0
                                Etop = aSim.Eapp(1);
                            else
                                Etop = aSim.Eapp; 
                            end
                        case 'neumann'                        
                            Etop(1) = (2*E(2) + 2*E(1+aSim.NX))/4.0;
                            Etop(NT) = (2*E(NT-1) + 2*E(NT+aSim.NX))/4.0;
                            for i = 2:NT-1
                                ptNum = topPts(i);
                                Eleft = E(ptNum - 1);
                                Eright = E(ptNum + 1);
                                Ebottom = E(ptNum+aSim.NX);
                                Etop(i) = (Eleft + Eright + (2.0*Ebottom))./4.0;                       
                            end
                    end                
                    E(topPts) = Etop;                
                    
                    % ===============
                    % Interior
                    % ===============                  
                    for i = 1:numel(Eint)
                        ptNum = intPts(i);
                        Eleft = E(ptNum-1);
                        Eright = E(ptNum+1);
                        Ebottom = E(ptNum+aSim.NX);
                        Etopper = E(ptNum-aSim.NX);
                        Eint(i) = (Eleft + Eright + Ebottom + Etopper)./4.0;
                    end
                    % EintPlot = reshape(Eint, [aSim.NX-2,aSim.NY-2]);
                    % figure(98)
                    % surf(EintPlot)
    
                    E(intPts) = Eint;
    
                    % ===============
                    % Convergence check!
                    % ===============
                    eps = abs(E(intPts) - ((E(intPts-1) + E(intPts+1) + E(intPts-aSim.NX) + E(intPts+aSim.NX))/4.0));
                    convCheck = eps(eps<=tol);
    
                    if numel(convCheck) == numel(E(intPts))
                        converged1 = true;
                        % aPhi = reshape(E, [aSim.NX,aSim.NY]);

                        % figure(99)
                        % surf(aPhi)                        
                        break;
                    end                
                    k = k + 1;

                end 
                if converged1 == true
                    aPhi = reshape(E, [aSim.NX,aSim.NY]);
                    % figure(98)
                    % surf(aPhi)
                    converged2 = true;
                    % if numel(epsL) == aSim.NXc && numel(epsR) == aSim.NXa
                    %     converged2 = true;
                    % end                    
                end

                L = L + 1;

            end
            
            % figure(99)
            % surf(aPhi)  

            if converged2 == true
                aPhi = reshape(E, [aSim.NX,aSim.NY]);
                % figure(99)
                % surf(aPhi)                                              
            end
        end  
%%%
% * newtonRaphson2 - Uses Newton's method to determine the electrode
% potential at the electrode boundary condition for a given potential at
% the node just outside the electrode and the Nernst potential inside the
% electrode.
        function EBOut = newtonRaphson2(E1,rtype,galvSim,h,xpos)
            %newtonRaphson2 - Uses Newton's method to obtain the electrode
            %boundary condition potential
            numIters = 30; 
            tolNR = 1.0e-3;
            aDelta = 1.0e-2;
            aDelta1 = 1.0e-4;
            EBOut = E1; 

            switch rtype
                case 'anodic' 
                    R = galvSim.bVreactions(1).C.R;
                    T = galvSim.bVreactions(1).T;
                    RT = R*T;
                    F = galvSim.bVreactions(1).C.F;
                    a = galvSim.bVreactions(1).aMetal.BetaMetalOxidation;
                    z = galvSim.bVreactions(1).aMetal.OxidationLevelZ;                        

                    Ec = galvSim.bVreactions(1).E0;
                    beta = RT/(a*z*F);
                    EBold = (1.0-aDelta)*Ec;
                    invkappa = 1.0/galvSim.bVreactions(1).kappa;                  
                    basedistx = galvSim.NXc * h;
                    for i = 1:numel(E1)     
                        Rpass = 0.0; 
                        R2 = invkappa/h;
                        R1 = R2 + Rpass;
                        for jdx = 1:numIters
                            if EBold > E1(i)
                                deltaE = EBold-E1(i);
                                iC = (deltaE/R1)*aDelta1; 
                                i0 = galvSim.bVreactions(1).i0;
                                logterm = log(iC/i0);
                                correction = Ec - beta*logterm;
                                fE = correction - EBold;
                                dFdE = -beta/deltaE - 1;
                                adj = fE/dFdE;
                                EBnew = EBold - adj;     
                            else
                                deltaE = E1(i)-EBold;  
                                iC = (deltaE/R1)*aDelta1;
                                i0 = galvSim.bVreactions(1).i0;
                                logterm = log(iC/i0);
                                correction = Ec - beta*logterm;
                                fE = correction - EBold;
                                dFdE = beta/deltaE - 1;
                                adj = fE/dFdE;
                                EBnew = EBold - adj;                                 
                            end
                            tolCheck = abs((EBnew - EBold)/EBold);            
                            if tolCheck <= tolNR
                                EBOut(i) = EBnew;                               
                                break;
                            else
                                EBold = EBnew;
                            end
                        end
                    end                                           
                case 'cathodic'
                    R = galvSim.bVreactions(2).C.R;
                    T = galvSim.bVreactions(2).T;
                    RT = R*T;
                    F = galvSim.bVreactions(2).C.F;
                    a = galvSim.bVreactions(2).aMetal.BetaORR;
                    z = galvSim.bVreactions(2).C.z_orr;                            
                    Ec = galvSim.bVreactions(2).E0;
                    beta = RT/(a*z*F);
                    EBold = (1.0-aDelta)*Ec;
                    invkappa = 1.0/galvSim.bVreactions(2).kappa;                       
                    for i = 1:numel(E1)
                        Rpass = 0.0; 
                        R2 = invkappa/h; 
                        R1 = R2 + Rpass; 
                        for jdx = 1:numIters                                               
                            if EBold > E1(i)                                
                                deltaE = EBold-E1(i);
                                iC = (deltaE/R1)*aDelta1;
                                i0 = galvSim.bVreactions(2).i0;
                                logterm = log(iC/i0);
                                correction = Ec - beta*logterm;
                                fE = correction - EBold;
                                dFdE = -beta/deltaE - 1;
                                adj = fE/dFdE;
                                EBnew = EBold - adj;                                 
                            else
                                deltaE = E1(i)-EBold;  
                                iC = (deltaE/R1)*aDelta1;
                                i0 = galvSim.bVreactions(2).i0;
                                logterm = log(iC/i0);
                                correction = Ec - beta*logterm;
                                fE = correction - EBold;
                                dFdE = beta/deltaE - 1;
                                adj = fE/dFdE;
                                EBnew = EBold - adj;                                 
                            end
                            tolCheck = abs((EBnew - EBold)/EBold);            
                            if tolCheck <= tolNR
                                EBOut(i) = EBnew;
                                break;
                            else
                                EBold = EBnew;
                            end
                        end
                    end                           
            end

        end        
%%%
% * plotPhis - Outputs contour plots of the equipotential lines throughout
% the depth of the electrolyte as well as cross-sections of the potential
% field at 2 different depths.        
        function plotPhis(aSim)
            %======================================================================
            tick_label_size = 16;
            axis_label_size = 18;
            title_label_size = 20;
            axis_line_width = 3;
            font_weight = 'bold';
            plot_line_width = 4;
            plot_line_width_2 = 2;
            marker_size = 8;
            colorVec = {[0 0 1],[0 0 0.95],[0 0 0.9]}; 
            colorVec1 = {[1 0 0],[0.95 0 0],[0.9 0 0]};   
            markerVec = {'bo','bo','bo','rs','rs','rs','k^','k^','k^','g<','g<','g<','c>','c>','c>'};
            markerVec1 = {'bo','rs','k^','g<','c>','o','s','^','<','>','o','s','^','<','>'};
            linVec = {'-b','--b','-.b','-r','--r','-.r','-k','--k', ...
                '-.k','none','none','none','none','none'};      
            someColors = {'#006400','#00BFFF','#00FF7F','#B0C4DE','#FF7F50','#6495ED','#4682B4','#0000FF','#FFD700'};
            %======================================================================       

            lengthConvert = 1.0e3; %m -> mm
            [X,Y] = meshgrid(aSim.xpos,aSim.ypos);
            fo = double(aSim.nFig);

            %=====================================================
            % Plot the 2D representation of the field
            %=====================================================            
            figure(fo);
            hold on            
            c = colorbar;
            c.Label.String = 'Potential (V_{SCE})';
            contour(X.*lengthConvert,Y.*lengthConvert,aSim.phi + aSim.potentialOffset,20,'LineWidth', plot_line_width); %[-0.188 -0.146, -0.104, -0.064, 0.0, 0.064, 0.104, 0.146, 0.188]
            ylim([0.0 aSim.H].*lengthConvert)
            xlabel('x (mm)', 'FontSize', axis_label_size,'FontWeight',font_weight)
            ylabel('y (mm)', 'FontSize', axis_label_size,'FontWeight',font_weight)
            ax = gca;
            ax.ZGrid = 'off';
            ax.XGrid = 'off';
            ax.YGrid = 'off';
            ax.FontName = 'Times New Roman';
            ax.FontSize = tick_label_size;
            ax.FontWeight = font_weight;
            ax.XAxis.LineWidth = axis_line_width;
            ax.YAxis.LineWidth = axis_line_width;    
            box on
            hold off  
            %===================================================== 
            xPot = aSim.phi(1,:) + aSim.potentialOffset;
            mid = floor(aSim.NY);
            xPot2 = aSim.phi(mid,:) + aSim.potentialOffset;
            ylocLabels = (aSim.bV(2).E0 + aSim.potentialOffset)-0.01;            
            %=====================================================
            % Plot cross-sections of the potential field
            %=====================================================
            figure(fo+1)
            hold on
            plot(aSim.xpos.*lengthConvert,xPot,'-.g')
            plot(aSim.xpos.*lengthConvert,xPot2,'-r')
            xline(0.0,'--k');            
            text(-5.0,ylocLabels,aSim.bV(2).reaction)
            text(10.0,ylocLabels,aSim.bV(1).reaction)
            xlabel('Position (mm)')
            ylabel('Potential (V_{SCE})')            
            ylim([aSim.bV(1).E0 + aSim.potentialOffset, aSim.bV(2).E0 + aSim.potentialOffset])            
            box on
            % legend('depth = 2 mm')
            % legend boxoff
            axis square
            hold off             
        end
        
    end    
end
