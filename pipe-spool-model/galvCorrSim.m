classdef galvCorrSim
    %galvCorrSim Summary of this class goes here
    %   Detailed explanation goes here

    properties        
        xpos double
        ypos double
        L double
        H double
        dx double
        dy double

        bcT char
        bcL char
        bcR char
        
        NX double
        NXa double
        NXc double
        NY double

        conductivity double

        bVreactions butlerVolmer

        Eapp double
        corrosionCurrentAnodic double
        corrosionCurrentCathodic double
        corrosionCurrentTotal double
        eCorr double
        phi double
        nFig int16
    end

    methods
        function obj = galvCorrSim(BCs,length,height,deltaL,eNodes,numNodes,aNodes,cNodes,reactAnode,reactCathode,Vapp,nF,env) %
            %galvCorrSim Construct an instance of this class
              % Detailed explanation goes here
              if nargin == 0
                obj.NX = 50;
                obj.NXa = 25;
                obj.NXc = 25;                
                obj.NY = 50;
    
                obj.L = 1.0;
    
                obj.bcL = 'neumann';
                obj.bcR = 'neumann';
                obj.bcT = 'neumann';
                obj.H = 10.0*obj.L;
    
                obj.dx = obj.L/(numNodes-1);
                obj.dy = obj.H/(numNodes-1);
                leftEnd = -obj.L/2.0;
                rightEnd = obj.L/2.0;
                obj.xpos = leftEnd:obj.dx:rightEnd;
                obj.ypos = 0.0:obj.dy:obj.H;
    
                obj.bVreactions(1) = 'anodic';
                obj.bVreactions(2) = 'cathodic';
    
                obj.Eapp = 0.0;
                obj.nFig = 132; 

              elseif nargin == 1
                  
              else
                obj.NX = numNodes;
                obj.NXa = aNodes;
                obj.NXc = cNodes;
                obj.NY = eNodes;
    
                obj.L = length;
                obj.H = height;
                obj.bcT = char(BCs{1});
                obj.bcL = char(BCs{2});
                obj.bcR = char(BCs{3});
    
                obj.dx = deltaL; %obj.L/(obj.NX-1);
                obj.dy = deltaL; %obj.H/(obj.NY-1); %deltaL; %

                leftEnd = -obj.NXc*obj.dx; %obj.L/2.0;
                rightEnd = obj.NXa*obj.dx; %obj.L/2.0;
                
                obj.xpos = leftEnd:obj.dx:rightEnd;
                obj.NXc = cNodes+1;
                obj.NX = numel(obj.xpos);

                obj.ypos = 0.0:obj.dy:obj.H;
                obj.NY = numel(obj.ypos);

                obj.bVreactions(1) = butlerVolmer(reactAnode,env(1),env(2),env(3),env(4));
                obj.bVreactions(2) = butlerVolmer(reactCathode,env(1),env(2),env(3),env(4));
    
                erange = -1.5:0.01:1.5;
                obj.Eapp = Vapp;
                [obj.corrosionCurrentTotal,obj.corrosionCurrentAnodic,obj.corrosionCurrentCathodic] = obj.GetTotalCurrent(erange);
                [mini,idxmin] = min(abs(obj.corrosionCurrentTotal));
                obj.eCorr = erange(idxmin);

                obj.nFig = nF;                  
              end
        end
        
        function [iSum,ia,ic] = GetTotalCurrent(obj,erange)
            ia = obj.bVreactions(1).anodeKineticsCuNi(erange);
            ic = obj.bVreactions(2).multiCathodicI625(erange);
            iSum = ia + ic;
        end
    
    end

    methods (Static)
        function aTafelPolCurve(aSim)
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
            vapp = -1.5:0.01:1.5;

            %=====================================================
            % Plots the modeled polarizatio curves
            %=====================================================
            figure(100)
            hold on

            % plot(abs(bv1.),vapp,'-b','LineWidth', plot_line_width)
            % plot(abs(i2),vapp,'-r','LineWidth', plot_line_width)

            plot(abs(aSim.corrosionCurrentAnodic),vapp,':r','LineWidth', plot_line_width-1)
            % plot(abs(i1a),vapp,':b','LineWidth', plot_line_width-1)                       

            plot(abs(aSim.corrosionCurrentCathodic),vapp,':b','LineWidth', plot_line_width-1)
            % plot(abs(i2a),vapp,':r','LineWidth', plot_line_width-1)   
            
            plot(abs(aSim.corrosionCurrentTotal),vapp,':k','LineWidth', plot_line_width-2)

            xlabel('i (A/cm^2)', 'FontSize', axis_label_size,'FontWeight',font_weight)
            ylabel('V_{Ag/AgCl}', 'FontSize', axis_label_size,'FontWeight',font_weight)
            ax = gca;
            ax.XScale = 'log';
            ax.FontName = 'Times New Roman';
            ax.FontSize = tick_label_size;
            ax.FontWeight = font_weight;
            ax.XAxis.LineWidth = axis_line_width;
            ax.YAxis.LineWidth = axis_line_width;            
            % 
            ylim([-1.5 0.5])
            xlim([10^-10 10^0])
            axis square
            box on
            legend(aSim.bVreactions(1).reaction,aSim.bVreactions(2).reaction,'location','best')
            legend boxoff
            hold off
            %=====================================================
            % figure(101)
            % hold on
            % plot(vapp,e1,'-b')
            % plot(vapp,e2,'-r')
            % hold off

        end        
    end
end