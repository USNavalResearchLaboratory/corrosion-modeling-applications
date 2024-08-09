function corrosionPotentialsEndTime(i,idxT,pipeLoopCathodic,pipeLoopAnodic, Bmodels)
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
    linVec = {':bo',':r^','-.b','-r','--r','-.r','-k','--k', ...
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
    titleString = strcat(pipeLoopCathodic.alloyType,{' '},'to',{' '},pipeLoopAnodic.alloyType,{' '},'Pipe',{' '}); %,num2str(i),'-endTimes'
    fNum = (100*i) + 3;
    f = figure(fNum);
    f.Position = [10 10 900 600];
    hold on
    title(titleString, 'FontSize', title_label_size,'FontWeight',font_weight)

    xpos = pipeLoopCathodic.Places(idxT,1:pipeLoopCathodic.numberOfReferenceElectrodes);
    ypos = pipeLoopCathodic.Potentials(idxT,1:pipeLoopCathodic.numberOfReferenceElectrodes);
    plot(xpos,ypos,linVec{1},'LineWidth', plot_line_width) %

    xpos = pipeLoopAnodic.Places(idxT,1:pipeLoopAnodic.numberOfReferenceElectrodes);
    ypos = pipeLoopAnodic.Potentials(idxT,1:pipeLoopAnodic.numberOfReferenceElectrodes);    
    plot(xpos,ypos,linVec{2},'LineWidth', plot_line_width)   

    N = numel(Bmodels.xElsyca);
    halfWay = round(N/2);
    if Bmodels.hasElsyca == true
        plot(Bmodels.xElsyca(1:halfWay),Bmodels.pipeLoopPotentialElsyca(1:halfWay),'--b','LineWidth', plot_line_width)
        plot(Bmodels.xElsyca(halfWay+1:N),Bmodels.pipeLoopPotentialElsyca(halfWay+1:N),'--r','LineWidth', plot_line_width)
    end

    N = numel(Bmodels.xSimple);
    halfWay = round(N/2);    
    % xAdj = Bmodels.xSimple(halfWay);
    if Bmodels.hasMyModel == true
        plot(Bmodels.xSimple(1:halfWay),Bmodels.pipeLoopPotentialSimpleModel(1:halfWay),'-b','LineWidth', plot_line_width)
        plot(Bmodels.xSimple(halfWay:N),Bmodels.pipeLoopPotentialSimpleModel(halfWay:N),'-r','LineWidth', plot_line_width)    
    end
    
    legend('Data - I625','Data - CuNi','I625-CP Master','CuNi-CP Master','I625-PC Model','CuNi-PC Model')
    legend boxoff

    xlabel('Pipe position (m)', 'FontSize', axis_label_size,'FontWeight',font_weight)
    ylabel('Potential (V_{Ag/AgCl})', 'FontSize', axis_label_size,'FontWeight',font_weight)

    xlim([-3.0 3.0])

    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = tick_label_size;
    ax.FontWeight = font_weight;
    ax.XAxis.LineWidth = axis_line_width;
    ax.YAxis.LineWidth = axis_line_width;

    box on    

    axis square
    hold off    

    saveas(gcf,char(titleString),'png')
end