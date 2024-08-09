function corrosionPotentialsWaterfallPlot(fNum,typePlot,time,pipeLoop,titleString)
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
    
    newcolorsCuNi625Full = {'#cce5ff','#99ccff','#66b2ff','#3399ff','#0080ff','#0066cc','#004c99','#003366', ...
        '#ffcccc','#ffcc99','#ffb266','#ff9933','#ff8000','#cc6600','#994c00','#663300','#331900'}; %'#cce5ff',  

    newcolorsCuNi625HX = {'#cce5ff','#99ccff','#66b2ff','#3399ff','#0080ff','#0066cc','#004c99', ...
        '#ffcccc','#ffcc99','#ffb266','#ff9933','#ff8000','#cc6600','#994c00','#663300','#331900'}; %'#003366','#cce5ff',  
    
    newcolorsCuNi625Short = {'#cce5ff','#99ccff','#66b2ff','#3399ff', ...
        '#ffcccc','#ffcc99','#ffb266','#ff9933','#ff8000','#cc6600','#994c00','#663300','#331900'}; %,'#0080ff','#0066cc','#004c99','#003366','#cce5ff', 
    
    newcolorsCuNiTi = {'#e5ccff','#cc99ff','#b266ff','#9933ff','#7f00ff','#6600cc','#4c0099','#330066', ...
        '#ffcccc','#ffcc99','#ffb266','#ff9933','#ff8000','#cc6600','#994c00','#663300','#331900'}; %'#cce5ff', 

    newcolorsCuNiTiShort = {'#e5ccff','#cc99ff','#b266ff','#9933ff', ...
        '#ffcccc','#ffcc99','#ffb266','#ff9933','#ff8000','#cc6600','#994c00','#663300','#331900'}; %'#7f00ff','#6600cc','#4c0099','#330066','#cce5ff',     
    %====================================================================== 
    figure(fNum)
    hold on
    title(titleString, 'FontSize', title_label_size,'FontWeight',font_weight)
    for j = 1:pipeLoop.numLocsAnode
        plot3(time,pipeLoop.anodePlaces(:,j),pipeLoop.anodePotentials(:,j),'LineWidth', plot_line_width) %'MarkerSize',marker_size-2,
    end    
    for i = 1:pipeLoop.numLocsCathode
        plot3(time,pipeLoop.cathodePlaces(:,i),pipeLoop.cathodePotentials(:,i),'LineWidth', plot_line_width) %
    end
  
    switch typePlot
        case 'I625'
            colororder(newcolorsCuNi625Full) 
        case 'I625HX'
            colororder(newcolorsCuNi625HX)             
        case 'I625Short'
            colororder(newcolorsCuNi625Short)    
        case 'Ti'
            colororder(newcolorsCuNiTi)             
        case 'TiShort'
            colororder(newcolorsCuNiTiShort) 
    end
          

    xlabel('Measurement time (d)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    ylabel('Pipe position (m)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    zlabel('Potential (V_{Ag/AgCl})', 'FontSize', axis_label_size,'FontWeight',font_weight)

%     zlim([-0.2 0.3])
    zlim([-0.3 0.2])

    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = tick_label_size;
    ax.FontWeight = font_weight;
    ax.XAxis.LineWidth = axis_line_width;
    ax.YAxis.LineWidth = axis_line_width;
    ax.ZAxis.LineWidth = axis_line_width;

    box on    
    view([200 20])
    hold off    
    saveas(gcf,titleString,'png')
end