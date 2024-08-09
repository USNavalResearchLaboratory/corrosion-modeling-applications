function outputxml(dNode, instNode,metalNode,electNode,dataNode,fn)
    
    docNode = com.mathworks.xml.XMLUtils.createDocument(dNode);
    %======================================================================
    institution_type_node = docNode.createElement('InstituteLocation');
    docNode.getDocumentElement.appendChild(institution_type_node);

    name_node = docNode.createElement('Name');
    name_text = docNode.createTextNode(instNode.name);
    name_node.appendChild(name_text);
    institution_type_node.appendChild(name_node);

    city_node = docNode.createElement('City');
    city_text = docNode.createTextNode(instNode.city);
    city_node.appendChild(city_text);    
    institution_type_node.appendChild(city_node);

    state_node = docNode.createElement('State');
    state_text = docNode.createTextNode(instNode.state);
    state_node.appendChild(state_text);    
    institution_type_node.appendChild(state_node);

    country_node = docNode.createElement('Country');
    country_text = docNode.createTextNode(instNode.country);
    country_node.appendChild(country_text);    
    institution_type_node.appendChild(country_node);    

    %======================================================================
    %======================================================================
    material_data_node = docNode.createElement('MaterialData');
    docNode.getDocumentElement.appendChild(material_data_node);

    uns_node = docNode.createElement('UNSCode');
    uns_text = docNode.createTextNode(metalNode.code);
    uns_node.appendChild(uns_text);
    material_data_node.appendChild(uns_node);

    cName_node = docNode.createElement('CommonName');
    cName_text = docNode.createTextNode(metalNode.name);
    cName_node.appendChild(cName_text);
    material_data_node.appendChild(cName_node);

    surf_node = docNode.createElement('SurfacePrep');
    surf_text = docNode.createTextNode(metalNode.surf);
    surf_node.appendChild(surf_text);
    material_data_node.appendChild(surf_node);    

    area_node = docNode.createElement('ExpArea');
    area_node.setTextContent(num2str(metalNode.area))
    area_node.setAttribute('units','m2'); 
    material_data_node.appendChild(area_node);    

    npoints_node = docNode.createElement('NdataPoints');
    npoints_node.setTextContent(num2str(metalNode.N))
    material_data_node.appendChild(npoints_node);       
    %======================================================================
    %======================================================================
    electrolyte_data_node = docNode.createElement('ElectrolyteData');
    docNode.getDocumentElement.appendChild(electrolyte_data_node);

    cl_node = docNode.createElement('ClConc');
    cl_node.setTextContent(num2str(electNode.clconc))
    cl_node.setAttribute('units','M'); 
    electrolyte_data_node.appendChild(cl_node);

    pH_node = docNode.createElement('pH');
    pH_node.setTextContent(num2str(electNode.pH))
    pH_node.setAttribute('units','unitless'); 
    electrolyte_data_node.appendChild(pH_node);

    o2_node = docNode.createElement('O2conc');
    o2_node.setTextContent(num2str(electNode.o2))
    o2_node.setAttribute('units','M'); 
    electrolyte_data_node.appendChild(o2_node);   

    s2_node = docNode.createElement('S2conc');
    s2_node.setTextContent(num2str(electNode.s2))
    s2_node.setAttribute('units','M'); 
    electrolyte_data_node.appendChild(s2_node);     

    temp_node = docNode.createElement('Temperature');
    temp_node.setTextContent(num2str(electNode.temp))
    temp_node.setAttribute('units','C'); 
    electrolyte_data_node.appendChild(temp_node);

    conductivity_node = docNode.createElement('Conductivity');
    conductivity_node.setTextContent(num2str(electNode.sigma))
    conductivity_node.setAttribute('units','S/m'); 
    electrolyte_data_node.appendChild(conductivity_node);    

    flow_node = docNode.createElement('Flow');
    flow_node.setTextContent(num2str(electNode.flow))
    flow_node.setAttribute('units','m/s'); 
    electrolyte_data_node.appendChild(flow_node);       
    %======================================================================
    %======================================================================     
    exp_mdl_data_node = docNode.createElement('Data');
    docNode.getDocumentElement.appendChild(exp_mdl_data_node);    

    for i = 1:metalNode.N

        pointNode = docNode.createElement('point');
        pointNode.setTextContent(num2str(i))
        pointNode.setAttribute('units','unitless'); 
        exp_mdl_data_node.appendChild(pointNode); 

        anodicNode = docNode.createElement('anodici');
        anodicNode.setTextContent(num2str(dataNode(i).ianode))
        anodicNode.setAttribute('units','A/m2'); 
        exp_mdl_data_node.appendChild(anodicNode);   

        orrNode = docNode.createElement('orri');
        orrNode.setTextContent(num2str(dataNode(i).iorr))
        orrNode.setAttribute('units','A/m2'); 
        exp_mdl_data_node.appendChild(orrNode);    
        
        herNode = docNode.createElement('heri');
        herNode.setTextContent(num2str(dataNode(i).iher))
        herNode.setAttribute('units','A/m2'); 
        exp_mdl_data_node.appendChild(herNode);  

        totNode = docNode.createElement('totali');
        totNode.setTextContent(num2str(dataNode(i).itot))
        totNode.setAttribute('units','A/m2'); 
        exp_mdl_data_node.appendChild(totNode);

        potNode = docNode.createElement('Vapp');
        potNode.setTextContent(num2str(dataNode(i).v))
        potNode.setAttribute('units','Vsce'); 
        exp_mdl_data_node.appendChild(potNode);         
    end
   


    xmlwrite(fn,docNode)


    
    % phone_number_node = docNode.createElement('PhoneNumber');
    % phone_number_text = docNode.createTextNode('(508) 647-7000');
    % phone_number_node.appendChild(phone_number_text);
    % pc_type_node.appendChild(phone_number_node);
    % 
    % address_node = docNode.createElement('Address');
    % address_node.setTextContent('3 Apple Hill Dr, Natick MA')
    % % set an attribute directly
    % address_node.setAttribute('type','work');  
    % pc_type_node.appendChild(address_node);    
    
end