%% Polarization Curve Class
% Steven A. Policastro, Ph.D. 
% Center for Corrosion Science and Engineering, 
% U.S. Naval Research Laboratory
% 4555 Overlook Avenue SW
% Washington, DC 20375 
% 
% 
% This class serves as the template to create an object that can calculate 
% and plot a polarization curve. 
%
% Alloys specifically included in this version of the code:
%
% * HY 80 
% * HY 100 
% * SS 316 
% * Ti
% * CuNi (Wrought)
% * I625 
% 
% *Note* This version of the code contains estimates for reaction
% properties for the alloys specified above.  If additional alloys 
% are to be considered, their model reaction properties must be
% entered below, or this version of the code must be changed so the
% values can be obtained elsewhere.
% _Record of Revisions_:
%
% * Created - October 2021
% * Last revision: 30-Jul-2024 
% 
% Contents
% 
%% Class definition
classdef PolarizationCurveModel
% PolarizationCurveModel - Class for instantiating, defining, and plotting 
% a model polarization curve for a specified alloy in specified exposure
% conditions.  
%%% Internal class properties 
% These properties are not accessible outside the class.  
   properties (SetAccess = private)
        tick_label_size = 16;
        axis_label_size = 18;
        title_label_size = 20;
        axis_line_width = 3;
        font_weight = 'bold';
        plot_line_width = 3;
        plot_line_width_2 = 2;
        outOfRangeVal = -1000;
        min2Deriv = 0.001; %0.004;
        offsetPeak = 900;
        cutOffPeak1Left = 0.02;
        cutOffPeak1Right = 0.5;
        minFilteredCurrentForORRPeak = 1.0e-6;
        minFileredCurrentForFeRedPeak = 1.0e-8;
        minActivationHERCurrent = 5.0e-5;
        lineTypes = {'--','-.',':','-'};
   end   
%%% Public class properties
% These properties are accessible by outside code, once an instance of this
% class has been initialized.  These properties cover the solution
% properties, environmental conditions, and currents and potential values
% for the polarization curve calculations.
   properties
        ItemNo int16
        Name char
        UNSname char
        solution naclSolutionChemistry
        Temp double
        cCl double
        pH double
        flowvelocity double          
        exposedArea double
        allCatRxnsModel ElectrochemicalReductionReaction
        inclCatRxnsModel reactionNames
        allAnRxnsModel ElectrochemicalOxidationReaction
        inclAnRxnsModel reactionNames        
        eAppModel double
        iTotModel double
        numCatReactions int32
        numAnReactions int32
        catReactsDict
        anReactsDict
        plotSymbol string        
   end
   methods (Access = public)
%%% Object constructor method       
    function obj = PolarizationCurveModel(no, fName, pdata, envcon) 
        %PolarizationCurveModel - Constructor method for this class.
        %
        %The purpose of this method is to construct an instance of the
        %PolarizationCurveModel class.  The constructor
        %requires:
        % =============================================
        % Name                              = char array
        % Applied potentials                = numeric array (V_{SCE})
        % Temperature                       = numeric value (^oC)
        % Chloride concentration            = numeric value ([M])
        % pH                                = numeric value
        % =============================================
        % The constructor returns an instance of the class
        % =============================================
        plotCurrentContributions = false;
        C = Constants;  
        obj.ItemNo = no;
        obj.Name = fName;          
        obj.Temp = envcon(1);
        obj.pH = envcon(2);
        obj.cCl = envcon(3);
        obj.flowvelocity = envcon(4);
        % =============================================
        % Setup the applied potential vector
        % =============================================       
        obj.eAppModel = pdata;    
        % =============================================
        % Calculate solution properties affected by chloride
        % concentration and temperature
        % =============================================   
        [cOH, cH] = Constants.calculatecHandcOH(obj.pH);
        obj.solution = naclSolutionChemistry(obj.cCl,obj.Temp);        
        % =============================================
        % Create the electrochemical reactions
        % ============================================= 
        % Create instances of corroding metal and define reaction
        % properties 
        % ============================================= 
        switch obj.Name
            case 'ss316'
                metal = SS316(obj.Name,obj.cCl,obj.Temp,obj.pH);                    
                unsmetal = 'UNS S31600';
                obj.UNSname = unsmetal;
                acolor = 'b';                
                % =====
                % ORR
                % =====
                cReact = [obj.solution.cO2,((obj.solution.aW/(1000))*C.M_H2O)^2]; %g/cm3
                cProd = [1.0, ((cOH/1000)*C.M_OH)^4]; %g/cm3               
                Dcoeff = obj.solution.dO2;
                obj.allCatRxnsModel(1) = ElectrochemicalReductionReaction(reactionNames.ORR, cReact, cProd, obj.Temp, C.z_orr, C.e0_orr_alk, Dcoeff, ...
                    obj.eAppModel, metal);  
                % =====    
                % HER 
                % =====
                cReact = [((obj.solution.aW/(1000))*C.M_H2O)^2, 1.0]; %g/cm3 55.55;
                cProd = [1.0, ((cOH/1000)*C.M_OH)^2]; %g/cm3 10.0^-(14.0-obj.pH); %mol/L
                Dcoeff = C.D_H2O;
                obj.allCatRxnsModel(2) = ElectrochemicalReductionReaction(reactionNames.HER, cReact, cProd, obj.Temp, C.z_her, C.e0_her_alk, Dcoeff, ...
                    obj.eAppModel,metal);     
                % =====    
                % Passivation
                % =====                     
                cReact = 1.0;
                cProd = 1.0e-6;
                obj.allAnRxnsModel(1) = ElectrochemicalOxidationReaction(reactionNames.Passivation, cReact, cProd, obj.Temp, obj.eAppModel, metal);
                % =====    
                % Pitting 
                % =====                     
                cReact = 1.0;
                cProd = 1.0e-6;
                obj.allAnRxnsModel(2) = ElectrochemicalOxidationReaction(reactionNames.Pitting, cReact, cProd, obj.Temp, obj.eAppModel, metal);
            case 'hy80'
                metal = HY80(obj.Name,obj.cCl,obj.Temp,obj.pH);
                unsmetal = 'UNS K31820';
                obj.UNSname = unsmetal;
                acolor = 'r';                                  
                % =====
                % ORR
                % =====
                cReact = [obj.solution.cO2,((obj.solution.aW/(1000))*C.M_H2O)^2]; %g/cm3
                cProd = [1.0, ((cOH/1000)*C.M_OH)^4]; %g/cm3               
                Dcoeff = obj.solution.dO2;
                obj.allCatRxnsModel(1) = ElectrochemicalReductionReaction(reactionNames.ORR, cReact, cProd, obj.Temp, C.z_orr, C.e0_orr_alk, Dcoeff, ...
                    obj.eAppModel, metal);   
                % =====    
                % HER 
                % =====
                cReact = [((obj.solution.aW/(1000))*C.M_H2O)^2, 1.0]; %g/cm3 55.55;
                cProd = [1.0, ((cOH/1000)*C.M_OH)^2]; %g/cm3 10.0^-(14.0-obj.pH); %mol/L
                Dcoeff = C.D_H2O;
                obj.allCatRxnsModel(2) = ElectrochemicalReductionReaction(reactionNames.HER, cReact, cProd, obj.Temp, C.z_her, C.e0_her_alk, Dcoeff, ...
                    obj.eAppModel,metal);     
                % =====    
                % Fe oxidation
                % =====                     
                cReact = 1.0;
                cProd = 1.0e-6;
                obj.allAnRxnsModel(1) = ElectrochemicalOxidationReaction(reactionNames.Fe_Ox, cReact, cProd, obj.Temp, obj.eAppModel, metal);
                % =====    
                % Pitting 
                % =====                     
                cReact = 1.0;
                cProd = 1.0e-6;
                obj.allAnRxnsModel(2) = ElectrochemicalOxidationReaction(reactionNames.Pitting, cReact, cProd, obj.Temp, obj.eAppModel, metal);
            case 'hy100'
                metal = HY100(obj.Name,obj.cCl,obj.Temp,obj.pH);
                unsmetal = 'UNS K32045';
                obj.UNSname = unsmetal;
                acolor = 'm';                               
                % =====
                % ORR
                % =====
                cReact = [obj.solution.cO2,((obj.solution.aW/(1000))*C.M_H2O)^2]; %g/cm3
                cProd = [1.0, ((cOH/1000)*C.M_OH)^4]; %g/cm3               
                Dcoeff = obj.solution.dO2;
                obj.allCatRxnsModel(1) = ElectrochemicalReductionReaction(reactionNames.ORR, cReact, cProd, obj.Temp, C.z_orr, C.e0_orr_alk, Dcoeff, ...
                    obj.eAppModel, metal);   
                % =====    
                % HER 
                % =====
                cReact = [((obj.solution.aW/(1000))*C.M_H2O)^2, 1.0]; %g/cm3 55.55;
                cProd = [1.0, ((cOH/1000)*C.M_OH)^2]; %g/cm3 10.0^-(14.0-obj.pH); %mol/L
                Dcoeff = C.D_H2O;
                obj.allCatRxnsModel(2) = ElectrochemicalReductionReaction(reactionNames.HER, cReact, cProd, obj.Temp, C.z_her, C.e0_her_alk, Dcoeff, ...
                    obj.eAppModel,metal);      
                % =====    
                % Fe oxidation
                % =====                     
                cReact = 1.0;
                cProd = 1.0e-6;
                obj.allAnRxnsModel(1) = ElectrochemicalOxidationReaction(reactionNames.Fe_Ox, cReact, cProd, obj.Temp, obj.eAppModel, metal);
                % =====    
                % Pitting 
                % =====                     
                cReact = 1.0;
                cProd = 1.0e-6;
                obj.allAnRxnsModel(2) = ElectrochemicalOxidationReaction(reactionNames.Pitting, cReact, cProd, obj.Temp, obj.eAppModel, metal); 
            case 'i625'
                metal = I625(obj.Name,obj.cCl,obj.Temp,obj.pH);
                unsmetal = 'UNS N06625';
                obj.UNSname = unsmetal;
                acolor = 'k';
                % =====
                % ORR
                % =====
                cReact = [obj.solution.cO2,((obj.solution.aW/(1000))*C.M_H2O)^2]; %g/cm3
                cProd = [1.0, ((cOH/1000)*C.M_OH)^4]; %g/cm3               
                Dcoeff = obj.solution.dO2;
                obj.allCatRxnsModel(1) = ElectrochemicalReductionReaction(reactionNames.ORR, cReact, cProd, obj.Temp, C.z_orr, C.e0_orr_alk, Dcoeff, ...
                    obj.eAppModel, metal);  
                % =====    
                % HER 
                % =====
                cReact = [((obj.solution.aW/(1000))*C.M_H2O)^2, 1.0]; %g/cm3 55.55;
                cProd = [1.0, ((cOH/1000)*C.M_OH)^2]; %g/cm3 10.0^-(14.0-obj.pH); %mol/L
                Dcoeff = C.D_H2O;
                obj.allCatRxnsModel(2) = ElectrochemicalReductionReaction(reactionNames.HER, cReact, cProd, obj.Temp, C.z_her, C.e0_her_alk, Dcoeff, ...
                    obj.eAppModel,metal);      
                % =====    
                % Passivation
                % =====                     
                cReact = 1.0;
                cProd = 1.0e-6;
                obj.allAnRxnsModel(1) = ElectrochemicalOxidationReaction(reactionNames.Passivation, cReact, cProd, obj.Temp, obj.eAppModel, metal);
            case 'ti'
                metal = Ti(obj.Name,obj.cCl,obj.Temp,obj.pH);
                unsmetal = 'UNS R50700';
                obj.UNSname = unsmetal;
                acolor = 'c';
                % =====
                % ORR
                % =====
                cReact = [obj.solution.cO2,((obj.solution.aW/(1000))*C.M_H2O)^2]; %g/cm3
                cProd = [1.0, ((cOH/1000)*C.M_OH)^4]; %g/cm3               
                Dcoeff = obj.solution.dO2;
                obj.allCatRxnsModel(1) = ElectrochemicalReductionReaction(reactionNames.ORR, cReact, cProd, obj.Temp, C.z_orr, C.e0_orr_alk, Dcoeff, ...
                    obj.eAppModel, metal);   
                % =====    
                % HER 
                % =====
                cReact = [((obj.solution.aW/(1000))*C.M_H2O)^2, 1.0]; %g/cm3 55.55;
                cProd = [1.0, ((cOH/1000)*C.M_OH)^2]; %g/cm3 10.0^-(14.0-obj.pH); %mol/L
                Dcoeff = C.D_H2O;
                obj.allCatRxnsModel(2) = ElectrochemicalReductionReaction(reactionNames.HER, cReact, cProd, obj.Temp, C.z_her, C.e0_her_alk, Dcoeff, ...
                    obj.eAppModel,metal);      
                % =====    
                % Passivation
                % =====                     
                cReact = 1.0;
                cProd = 1.0e-6;
                obj.allAnRxnsModel(1) = ElectrochemicalOxidationReaction(reactionNames.Passivation, cReact, cProd, obj.Temp, obj.eAppModel, metal);
            case 'cuni'
                metal = CuNi(obj.Name,obj.cCl,obj.Temp,obj.pH,obj.flowvelocity);
                unsmetal = 'UNS C71500';
                obj.UNSname = unsmetal;                
                acolor = 'g';
                % =====
                % ORR
                % =====
                cReact = [obj.solution.cO2,((obj.solution.aW/(1000))*C.M_H2O)^2]; %g/cm3
                cProd = [1.0, ((cOH/1000)*C.M_OH)^4]; %g/cm3               
                Dcoeff = obj.solution.dO2;
                obj.allCatRxnsModel(1) = ElectrochemicalReductionReaction(reactionNames.ORR, cReact, cProd, obj.Temp, C.z_orr, C.e0_orr_alk, Dcoeff, ...
                    obj.eAppModel, metal);   
                % =====    
                % HER 
                % =====
                cReact = [((obj.solution.aW/(1000))*C.M_H2O)^2, 1.0]; %g/cm3 55.55;
                cProd = [1.0, ((cOH/1000)*C.M_OH)^2]; %g/cm3 10.0^-(14.0-obj.pH); %mol/L
                Dcoeff = C.D_H2O;
                obj.allCatRxnsModel(2) = ElectrochemicalReductionReaction(reactionNames.HER, cReact, cProd, obj.Temp, C.z_her, C.e0_her_alk, Dcoeff, ...
                    obj.eAppModel,metal);      
                % =====    
                % Cu oxidation
                % =====                     
                cReact = 1.0;
                cProd = 1.0e-1;
                obj.allAnRxnsModel(1) = ElectrochemicalOxidationReaction(reactionNames.Cu_Ox, cReact, cProd, obj.Temp, obj.eAppModel, metal);                  
        end
        obj.plotSymbol = strcat(char(obj.lineTypes{obj.ItemNo}),acolor);
        if plotCurrentContributions == true
            obj.PlotReactionCurrents();
        end
        obj.iTotModel = zeros(size(obj.eAppModel));
        for j = 1:numel(obj.allCatRxnsModel)
            obj.iTotModel = obj.iTotModel + obj.allCatRxnsModel(j).i;
        end
        for j = 1:numel(obj.allAnRxnsModel)
            obj.iTotModel = obj.iTotModel + obj.allAnRxnsModel(j).i;
        end
        outXML = true;
        if outXML
            dNode = 'PolarizationCurve';
            instNode.name = 'NRL';
            instNode.city = 'Washington';
            instNode.state = 'DC';
            instNode.country = 'USA';
            
            metalNode.code = unsmetal;
            metalNode.name = obj.Name;
            metalNode.surf = '600 grit SiC';
            metalNode.area = 1.0/(100*100);
                        
            electNode.clconc = obj.cCl;
            electNode.temp = obj.Temp;
            electNode.pH = obj.pH;
            electNode.o2 = obj.solution.cO2;
            electNode.sigma = 1.0/obj.solution.rhoNaCl;
            electNode.flow = 0.25;
            electNode.s2 = 0.08;

            metalNode.N = numel(obj.eAppModel);
            dataNode(metalNode.N).ianode = 0.0;
            dataNode(metalNode.N).iorr = 0.0;
            
            dataNode(metalNode.N).iher = 0.0;
            dataNode(metalNode.N).itot = 0.0; 
            dataNode(metalNode.N).v = 0.0;        
            for j = 1:metalNode.N
                if numel(obj.allAnRxnsModel) > 1
                    dataNode(j).ianode = abs(obj.allAnRxnsModel(1).i(j)) + abs(obj.allAnRxnsModel(2).i(j));
                else
                    dataNode(j).ianode = abs(obj.allAnRxnsModel(1).i(j));
                end                    
                dataNode(j).iorr = obj.allCatRxnsModel(1).i(j); % + obj.allCatRxnsModel(2).i(j)
                dataNode(j).iher = obj.allCatRxnsModel(2).i(j);
                dataNode(j).itot = dataNode(j).iorr + dataNode(j).iher;        
                dataNode(j).v = obj.eAppModel(j);
            end
            fn = strcat(obj.Name,'.xml');
            outputxml(dNode,instNode,metalNode,electNode,dataNode,fn)                
        end            
    end       
%%% Plotting method during debugging when addiing new materials   
% Internal plotting method for checking the currents arising from
% individual reactions.  This helps in debugging because currents that are
% non-physical can be difficult to detect if they are much smaller than
% other contributions.
    function PlotReactionCurrents(obj)
        figure(500)
        hold on
        plot(abs(obj.allAnRxnsModel(1).i),obj.eAppModel,'-k','LineWidth', obj.plot_line_width-2)
        plot(abs(obj.allAnRxnsModel(2).i),obj.eAppModel,'-.b','LineWidth', obj.plot_line_width-2)
        plot(abs(obj.allCatRxnsModel(1).i),obj.eAppModel,'--r','LineWidth', obj.plot_line_width-2)
        plot(abs(obj.allCatRxnsModel(2).i),obj.eAppModel,'-.g','LineWidth', obj.plot_line_width-2)            
        axis square
        box on
        ylim([-1.4,0.0])
        xlim([1.0e-13,0.01])
        xlabel('Current density (A/cm^2)', 'FontSize', obj.axis_label_size,'FontWeight',obj.font_weight)
        ylabel('Potential (V_{SCE})', 'FontSize', obj.axis_label_size,'FontWeight',obj.font_weight)
        ax = gca;
        ax.XScale = 'log';
        ax.FontSize = obj.tick_label_size;
        ax.FontWeight = obj.font_weight;
        ax.LineWidth = obj.axis_line_width;
        ax.XTick = [1.0e-13,1.0e-12,1.0e-11,1.0e-10,1.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-4,1.0e-3,0.01,0.1];
        ax.YTick = -1.4:0.2:0.1;       
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        legend('passive','pit','ORR','HER')
        legend boxoff
        hold off            
       end
   end
%%% Public static methods
% These methods, for plotting figures and file outputs, do not require 
% an instantiated class object to exist.
    methods (Static)
        
        function Plot_Polarization_Data_and_Model(obj,fig_num) %,aD
            sName = strcat(char(obj.UNSname), ', T =',num2str(obj.Temp),'K, [Cl^-] = ',num2str(obj.cCl),'M, pH = ',num2str(obj.pH));
            Plot_name = strcat(sName,'_',num2str(fig_num));
            h1 = figure(fig_num);
            set(h1,'Position', [10 20 800 800])
            hold on
            % Plot the average polarization data as a thin black dashed
            % line
            plot(abs(obj.iTotModel(1:length(obj.eAppModel))),obj.eAppModel,obj.plotSymbol,'LineWidth', obj.plot_line_width-1)
            axis square
            box on
            ylim([-1.5,0.5])
            xlim([1.0e-11,1.0e2])
            xlabel('Current density (A/cm^2)', 'FontSize', obj.axis_label_size,'FontWeight',obj.font_weight)
            ylabel('Potential (V_{SCE})', 'FontSize', obj.axis_label_size,'FontWeight',obj.font_weight)            
            ax = gca;
            ax.XScale = 'log';
            ax.FontSize = obj.tick_label_size;
            ax.FontWeight = obj.font_weight;
            ax.LineWidth = obj.axis_line_width;
            ax.XTick = [1.0e-11,1.0e-9,1.0e-7,1.0e-5,1.0e-3,1.0e-1,1.0e1];     
            ax.YTick = -1.5:0.2:0.5;   
            ax.XMinorTick = 'on';
            ax.YMinorTick = 'on';            
            legendString = {sName};
            legend(legendString,'Location','best')
            legend boxoff
            exportgraphics(ax,strcat(char(Plot_name),'.png'),'Resolution',300)
            hold off 
        end        
        
        function OutputFitValues(fileName,TC,cCl,fVals)
            if isfile(fileName)
                delete(char(fileName));
            end   
            writecell({'T','Cl-','dG_Cathodic','dG_Anodic','alpha','Diffusion_Length'},char(fileName),'Range','A1:F1')
            writematrix([TC,cCl],char(fileName),'Range','A2:B2')
            writematrix(fVals,char(fileName),'Range','C2:F5')
        end
    end

end
%------------- END OF CODE --------------
