classdef PipeLoopExperiments
    %PipeLoopConstants stores data about the pipe loop potentials and
    %corrosion currents
    %   A class to keep properties related to the pipe loop experiments.
    %   Loads data from a CSV file and stores the potentials and currents
    %   in different structures.

    properties
        F = 96485; %coul/mol
        R = 8.314; % J/mol K
        convLtom3 = 0.001;   
        convTctoTK = 273.15;
        convHrtoSec = 3600;
        convAgAgCltoSHE = 0.199;
        convModeltoExp = 0.097558;
        convmAtoA = 1.0e-3;
        time
        corrCurrent
        pipeLoopsCathodic aPipeLoop
        pipeLoopsAnodic aPipeLoop
    end

    methods
        function obj = PipeLoopExperiments(filename)
            %PipeLoopConstants Construct an instance of PipeLoopConstants
            %   Constructor that accepts a filename.  Needs to be a CSV
            %   file, which is read using the readtable function, and
            %   assumes there are 6 header lines in the file.  The time and
            %   currents are stored in class arrays while the potentials
            %   are stored in 'aPipeLoop' structures.
            %==============================
            % Load the data
            %==============================
            allData = readtable(filename,'NumHeaderLines',6); %
            %==============================
            obj.time = allData.Var2;    
            %==============================
            % Galvanic current between I625 Pipe 1 and CuNi Pipe 1
            %==============================
            corrCurrent = [allData.Var3,allData.Var20,allData.Var37,allData.Var64,allData.Var81, ...
                allData.Var98, allData.Var126,allData.Var143,allData.Var160];
            corrCurrent(corrCurrent==0) = nan;  
            obj.corrCurrent = corrCurrent;
            %==============================
            % Corrosion potentials along Pipe Loop 1
            %=============================
            % Pipe loop 1
            v1 = [allData.Var4,allData.Var5,allData.Var6,allData.Var7,allData.Var8, ...
                allData.Var9,allData.Var10,allData.Var11];
            v2 = [allData.Var12,allData.Var13,allData.Var14,allData.Var15,allData.Var16, ...
                allData.Var17,allData.Var18];

            %==============================
            % Corrosion potentials along Pipe Loop 2
            %=============================
            v3 = [allData.Var21,allData.Var22,allData.Var23,allData.Var24,allData.Var25, ...
                allData.Var26,allData.Var27]; %,allData.Var28
            v4 = [allData.Var29,allData.Var30,allData.Var31,allData.Var32,allData.Var33, ...
                allData.Var34,allData.Var35];   
            
            %==============================
            % Corrosion potentials along Pipe Loop 3
            %=============================
            v5 = [allData.Var38,allData.Var39,allData.Var40,allData.Var41,allData.Var42, ...
                allData.Var43,allData.Var44,allData.Var45];
            v6 = [allData.Var46,allData.Var47,allData.Var48,allData.Var49,allData.Var50, ...
                allData.Var51,allData.Var52,allData.Var53];   
            %==============================
            % Corrosion potentials along Pipe Loop 4
            %=============================
            v7 = [allData.Var65,allData.Var66,allData.Var67,allData.Var68,allData.Var69, ...
                allData.Var70,allData.Var71,allData.Var72];
            v8 = [allData.Var73,allData.Var74,allData.Var75,allData.Var76,allData.Var77, ...
                allData.Var78,allData.Var79,allData.Var80];    
            %==============================
            % Corrosion potentials along Pipe Loop 5
            %=============================
            v9 = [allData.Var82,allData.Var83,allData.Var84,allData.Var85,allData.Var86, ...
                allData.Var87,allData.Var88,allData.Var89];
            v10 = [allData.Var90,allData.Var91,allData.Var92,allData.Var93,allData.Var94, ...
                allData.Var95,allData.Var96,allData.Var116];   
            %==============================
            % Corrosion potentials along Pipe Loop 6
            %=============================
            v11 = [allData.Var99,allData.Var100,allData.Var101,allData.Var102,allData.Var103, ...
                allData.Var104,allData.Var105,allData.Var106];
            v12 = [allData.Var107,allData.Var108,allData.Var109,allData.Var110,allData.Var111, ...
                allData.Var112,allData.Var113,allData.Var114];     
            %==============================
            % Corrosion potentials along Pipe Loop 7
            %=============================
            v13 = [allData.Var127,allData.Var128,allData.Var129,allData.Var130,allData.Var131, ...
                allData.Var132,allData.Var133]; %
            v14 = [allData.Var135,allData.Var136,allData.Var137,allData.Var138,allData.Var139, ...
                allData.Var140,allData.Var141,allData.Var142];      
            %==============================
            % Corrosion potentials along Pipe Loop 8
            %=============================
            v15 = [allData.Var144,allData.Var145,allData.Var146,allData.Var147];
            v16 = [allData.Var152,allData.Var153,allData.Var154,allData.Var155,allData.Var156, ...
                allData.Var157,allData.Var158, allData.Var159];   
            %==============================
            % Corrosion potentials along Pipe Loop 9
            %=============================
            v17 = [allData.Var161,allData.Var162,allData.Var163,allData.Var164];
            v18 = [allData.Var169,allData.Var170,allData.Var171,allData.Var172,allData.Var173, ...
                allData.Var174,allData.Var175, allData.Var176];   

            VsCathodic = [v1,v3,v5,v7,v9,v11,v13,v15,v17];
            VsAnodic = [v2,v4,v6,v8,v10,v12,v14,v16,v18];

            cLocs = -[0,0.1,10,20,40,60,80,110];
            aLocs = [10,20,30,40,60,80,105];
            cL1 = numel(cLocs);
            aL1 = numel(aLocs);
            obj.pipeLoopsCathodic(1) = aPipeLoop(1,'I625',10,4,cLocs,12,-[0.1,10],VsCathodic(:,1:cL1));
            obj.pipeLoopsAnodic(1) = aPipeLoop(1,'CuNi',10,4,aLocs,12,10,VsAnodic(:,1:aL1));

            cLocs = -[0,10,20,40,60,80,110];
            aLocs = [10,20,30,40,60,80,105];
            cL2 = numel(cLocs) + cL1;
            aL2 = numel(aLocs) + aL1;            
            obj.pipeLoopsCathodic(2) = aPipeLoop(2,'I625',10,4,cLocs,3,[],VsCathodic(:,cL1+1:cL2));
            obj.pipeLoopsAnodic(2) = aPipeLoop(2,'CuNi',10,4,aLocs,3,[],VsAnodic(:,aL1+1:aL2));

            cLocs = -[0,6,10,20,40,60,80,114];
            aLocs = [6,10,20,30,40,60,80,114];        
            cL3 = numel(cLocs) + cL2;
            aL3 = numel(aLocs) + aL2;             
            obj.pipeLoopsCathodic(3) = aPipeLoop(3,'I625',10,2,cLocs,9,[],VsCathodic(:,cL2+1:cL3));
            obj.pipeLoopsAnodic(3) = aPipeLoop(3,'CuNi',10,2,aLocs,9,[],VsAnodic(:,aL2+1:aL3)); 

            cLocs = -[0,6,10,20,40,60,80,114];
            aLocs = [6,10,20,30,40,60,80,114];    
            cL4 = numel(cLocs) + cL3;
            aL4 = numel(aLocs) + aL3;             
            obj.pipeLoopsCathodic(4) = aPipeLoop(4,'I625',10,1.5,cLocs,3,[],VsCathodic(:,cL3+1:cL4));
            obj.pipeLoopsAnodic(4) = aPipeLoop(4,'CuNi',10,1.5,aLocs,3,[],VsAnodic(:,aL3+1:aL4)); 

            cLocs = -[0,6,10,20,40,60,80,114];
            aLocs = [6,10,20,30,40,60,80,114];   
            cL5 = numel(cLocs) + cL4;
            aL5 = numel(aLocs) + aL4;            
            obj.pipeLoopsCathodic(5) = aPipeLoop(5,'Ti',10,2.0,cLocs,3,[],VsCathodic(:,cL4+1:cL5));
            obj.pipeLoopsAnodic(5) = aPipeLoop(5,'CuNi',10,2.0,aLocs,3,[],VsAnodic(:,aL4+1:aL5));             

            cLocs = -[0,6,10,20,40,60,80,114];
            aLocs = [6,10,20,30,40,60,80,114];    
            cL6 = numel(cLocs) + cL5;
            aL6 = numel(aLocs) + aL5;            
            obj.pipeLoopsCathodic(6) = aPipeLoop(6,'I625',10,2.0,cLocs,3,-[10,40],VsCathodic(:,cL5+1:cL6));
            obj.pipeLoopsAnodic(6) = aPipeLoop(6,'CuNi',10,2.0,aLocs,3,[10,40],VsAnodic(:,aL5+1:aL6));     

            cLocs = -[0,20,25,30,35,38,52];
            aLocs = [6,10,20,30,40,60,80,114];    
            cL7 = numel(cLocs) + cL6;
            aL7 = numel(aLocs) + aL6;            
            obj.pipeLoopsCathodic(7) = aPipeLoop(7,'I625HX',5,2.0,cLocs,3,[],VsCathodic(:,cL6+1:cL7));
            obj.pipeLoopsAnodic(7) = aPipeLoop(7,'CuNi',10,2.0,aLocs,3,[],VsAnodic(:,aL6+1:aL7));            

            cLocs = [0,6,10,19];
            aLocs = -[6,10,20,30,40,60,80,114]; 
            cL8 = numel(cLocs) + cL7;
            aL8 = numel(aLocs) + aL7;             
            obj.pipeLoopsCathodic(8) = aPipeLoop(8,'I625Short',2,2.0,cLocs,3,-[6,19],VsCathodic(:,cL7+1:cL8));
            obj.pipeLoopsAnodic(8) = aPipeLoop(8,'CuNi',10,2.0,aLocs,3,[10,40],VsAnodic(:,aL7+1:aL8)); 

            cLocs = -[0,6,10,18.5];
            aLocs = [6,10,20,30,40,60,80,114]; 
            cL9 = numel(cLocs) + cL8;
            aL9 = numel(aLocs) + aL8;              
            obj.pipeLoopsCathodic(9) = aPipeLoop(9,'TiShort',10,2.0,cLocs,3,[],VsCathodic(:,cL8+1:cL9));
            obj.pipeLoopsAnodic(9) = aPipeLoop(9,'CuNi',10,2.0,aLocs,3,[],VsAnodic(:,aL8+1:aL9));                
        end
    end
end