function plotBodeEIS(fignum,xDataForPlotting_1,yDataForPlotting_1,yDataForPlotting_2,legendstring)
    % =====================================================================
    % Plotting variables
    tick_label_size = 20;
    axis_label_size = 24;
    title_label_size = 20;
    plot_line_width = 3;
    axis_line_width = 3;
    marker_size = 4;
    font_weight = 'bold';
    lines = {'-b','-k','--r'};
    chars1 = {'bo','r+','g*','k.','cx','ms','yd','b^','rv','g<','k>'};
    % =====================================================================
    lowVal = -1000;

    figure(fignum)
    
    t = tiledlayout(1,2);
    % Zmod vs freq
    axZmod = nexttile(t);
    xplots = xDataForPlotting_1(:,1);
    yplots = yDataForPlotting_1(:,1); 
    plot(axZmod,xplots,yplots,char(chars1{1}),'MarkerSize', marker_size,'LineWidth',  plot_line_width)
    hold(axZmod,'on')
    for i = 2:size(xDataForPlotting_1,2) %2:
        xplots = xDataForPlotting_1(:,i);
        yplots = yDataForPlotting_1(:,i);                               
        plot(axZmod,xplots,yplots,char(lines{i}),'LineWidth',  plot_line_width)
    end
    xlabel(axZmod,'Frequency (Hz)')
    ylabel(axZmod,'Z_{mod} (\Omega)')

    xlim(axZmod,[1.0e-3 1.0e6])
    ylim(axZmod,[1.0e2 1.0e11])
    xticks(axZmod,[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5])
    xticklabels(axZmod,{'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2','10^3','10^4','10^5'})          

    box(axZmod,'on')
    

    axZmod.FontName = 'times';
    axZmod.XScale = 'log';
    axZmod.YScale = 'log';
    axZmod.FontSize =  tick_label_size;
    axZmod.FontWeight =  font_weight;
    axZmod.LineWidth  =  axis_line_width;
    axZmod.XMinorTick = "on";
    axZmod.YMinorTick = "on";  

    legend(axZmod,legendstring,'Location','best')
    legend(axZmod,'boxoff') 

    axis(axZmod,'square')     
    hold(axZmod,'off')

    % Phase vs freq
    axZphz = nexttile(t);
    xplots = xDataForPlotting_1(:,1);
    yplots = yDataForPlotting_2(:,1);        
    plot(axZphz,xplots,yplots,char(chars1{1}),'MarkerSize', marker_size,'LineWidth',  plot_line_width)    
    hold(axZphz,'on')
    
    for i = 2:size(xDataForPlotting_1,2) %1:
        xplots = xDataForPlotting_1(:,i);
        yplots = yDataForPlotting_2(:,i);                
        plot(axZphz,xplots,yplots,char(lines{i}),'LineWidth', plot_line_width)        
    end

    xlabel(axZphz,'Frequency (Hz)')
    ylabel(axZphz,'Phase (^o)')

    xlim(axZphz,[1.0e-3 1.0e6])
    ylim(axZphz,[-90 0])

    xticks(axZphz,[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5])
    xticklabels(axZphz,{'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2','10^3','10^4','10^5'})

    box(axZphz,'on')
    % ax = gca;
    axZphz.FontName = 'times';
    axZphz.XScale = 'log';            
    axZphz.FontSize =  tick_label_size;
    axZphz.FontWeight =  font_weight;
    axZphz.LineWidth  =  axis_line_width;
    axZphz.XMinorTick = "on";
    axZphz.YMinorTick = "on";  

    axis(axZphz,'square')     
    hold(axZphz,'off')  
end
