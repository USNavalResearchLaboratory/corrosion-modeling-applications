function corrosionCurrentsPlot(figNum,time,curr,legendString,titleString)
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
    figure(figNum)
    hold on
    title(titleString, 'FontSize', title_label_size,'FontWeight',font_weight)
    numCurr = size(curr);
    
    for i = 1:numCurr(2)
        nums = ones(size(curr(:,i))).*i;
        plot3(time,nums,curr(:,i),'LineWidth', plot_line_width) %linVec{i},
    end
    colororder(someColors) 
    axis square
    box on

    xlabel('Measurement time (d)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    ylabel('Pipe loop', 'FontSize', axis_label_size,'FontWeight',font_weight)
    zlabel('Current (A/cm^2)', 'FontSize', axis_label_size,'FontWeight',font_weight)

    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = tick_label_size;
    ax.FontWeight = font_weight;
    ax.XAxis.LineWidth = axis_line_width;
    ax.YAxis.LineWidth = axis_line_width;
    
    zlim([0.0 50.0].*1.0e-6)

    legend(legendString)
    view([200 20])
    hold off

    saveas(gcf,titleString,'png')
end