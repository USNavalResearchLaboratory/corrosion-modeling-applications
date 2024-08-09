function corrosionPotentialsLeaklessVsSilverChloride(fNum,catNums, anNums, time, cathode, anode, lgdString,typePlot,titleString)
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
    linVecCat = {'-c','-b','--b','--c'}; 
    linVecAn = {'-k','-.k','--k',':k'};

    linVecCat2 = {'co','bo','b^','c^'}; 
    linVecAn2 = {'ko','k^','ks','k>'};    
    newcolorsCuNi625 = {'#cce5ff','#99ccff','#66b2ff','#3399ff','#0080ff','#0066cc','#004c99','#003366', ...
        '#ffcccc','#ffcc99','#ffb266','#ff9933','#ff8000','#cc6600','#994c00','#663300','#331900'}; %'#cce5ff',  

    newcolorsCuNiTi = {'#e5ccff','#cc99ff','#b266ff','#9933ff','#7f00ff','#6600cc','#4c0099','#330066', ...
        '#ffcccc','#ffcc99','#ffb266','#ff9933','#ff8000','#cc6600','#994c00','#663300','#331900'}; %'#cce5ff',      
    %====================================================================== 
    figure(fNum)
    hold on
    title(titleString, 'FontSize', title_label_size,'FontWeight',font_weight)
    
    switch typePlot
        case 'DirectCompare'
            yyaxis left
            plot(time,cathode.Potentials(:,catNums(1)),linVecCat2{2},'LineWidth', plot_line_width)
            plot(time,cathode.Potentials(:,catNums(2)),linVecCat2{3},'LineWidth', plot_line_width)
            plot(time,cathode.Potentials(:,catNums(3)),linVecCat2{1},'LineWidth', plot_line_width)
            plot(time,anode.Potentials(:,anNums(1)),linVecAn2{1},'LineWidth', plot_line_width)           
            ylabel('Potential (V_{Ag/AgCl})', 'FontSize', axis_label_size,'FontWeight',font_weight)
            ylim([-0.5 0.3])

            yyaxis right
            deltaV = cathode.Potentials(:,catNums(3)) - cathode.Potentials(:,catNums(1));
            plot(time,deltaV.*1000,'-r','LineWidth', plot_line_width)
            ylabel('\Delta V (mV)', 'FontSize', axis_label_size,'FontWeight',font_weight)
            ylim([0.0 200.0])

            ax = gca;
            ax.FontName = 'Times New Roman';
            ax.FontSize = tick_label_size;
            ax.FontWeight = font_weight;
            ax.XAxis.LineWidth = axis_line_width;
            ax.YAxis(1).LineWidth = axis_line_width;            
            ax.YAxis(2).LineWidth = axis_line_width;
        case 'NoCompare'
            for i = 1:numel(catNums)-1
                plot(time,cathode.Potentials(:,catNums(i)),linVecCat{i},'LineWidth', plot_line_width)
            end
        
            for i = 1:numel(anNums)
                plot(time,anode.Potentials(:,anNums(i)),linVecAn{i},'LineWidth', plot_line_width)
            end    
        
            plot(time,cathode.Potentials(:,catNums(numel(catNums))),'-g','LineWidth', plot_line_width) 
            ylabel('Potential (V_{Ag/AgCl})', 'FontSize', axis_label_size,'FontWeight',font_weight)
            ylim([-0.5 0.3])

            ax = gca;
            ax.FontName = 'Times New Roman';
            ax.FontSize = tick_label_size;
            ax.FontWeight = font_weight;
            ax.XAxis.LineWidth = axis_line_width;
            ax.YAxis.LineWidth = axis_line_width;              
    end

    xlabel('Measurement time (d)', 'FontSize', axis_label_size,'FontWeight',font_weight)        
    legend(lgdString)

    box on
    hold off
    saveas(gcf,titleString,'png')
end