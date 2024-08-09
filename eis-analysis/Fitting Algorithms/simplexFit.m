classdef simplexFit
%simplexFit - class for performing nonlinear regression using the simplex
%technique
%
% This class contains a constructor and methods for performing the
% regression analysis.
%
%
%==========================================================================
% Author:   Steve Policastro, Ph.D., Materials Science
% Center for Corrosion Science and Engineering, U.S. Naval Research
% Laboratory
% email address: steven.policastro@nrl.navy.mil  
% Website: 
% June 2022; Last revision: 13 August 2022
%==========================================================================

    properties        
        fitFun
        cLimits
        iters
    end

    methods
        function obj = simplexFit(funF,conLim)
        %simplex - Constructor for the simplex class
        %
        % Class constructor.
        %
        % Syntax:  
        %
        % Inputs: 
        % fun = handle to the function to be fit
        % cL = 2 colum matrix with the upper and lower limits on the b
        % parameters
        %
        % Outputs: None
        %
        % Other m-files required: None
        % Subfunctions: none
        % MAT-files required: none
        %
        % See also: 
        %
        %==========================================================================
        % Author:   Steve Policastro, Ph.D., Materials Science
        % Center for Corrosion Science and Engineering, U.S. Naval Research
        % Laboratory
        % email address: steven.policastro@nrl.navy.mil  
        % Website: 
        % June 2022; Last revision: 13 August 2022
        %==========================================================================            
            obj.fitFun = funF;
            if ~isempty(conLim)
                obj.cLimits = conLim;
            else 
                obj.cLimits = 0.0;
            end
           
        end
        function [bFit,MSE,fParams] = fitFn(obj,data,b,aParams)
        % fit - nonlinear regression function fitting uisng the simplex
        % algorithm
        %
        % Fits a function to a dataset using the simplex algorithm.
        %
        % Syntax:  [bfit,MSE,fParams] = fit(x,yR,yI,b,aParams)
        %
        % Inputs: 
        % data(1) = independent variable
        % data(2) = dependent variable
        % ** OR **
        % data(2) = real component of the dependent variable, and 
        % data(3) = imaginary component of the dependent variable    
        %     
        % b = vector of fit parameters      
        % aParams = vector of 1s and 0s of the same length as b, indicating whether
        % the parameter in b is active, or not, respsectively 
        %  
        %
        % Outputs: 
        % bfit = vector of the final fit parameters
        % MSE = mean square erro of the final fit iterations
        % fParams = vector of the model parameters
        %
        %
        % Other m-files required: objectiveFunction, class must be
        % instantiated
        % Subfunctions: none
        % MAT-files required: none
        %
        % See also: 
        %
        %==========================================================================
        % Author:   Steve Policastro, Ph.D., Materials Science
        % Center for Corrosion Science and Engineering, U.S. Naval Research
        % Laboratory
        % email address: steven.policastro@nrl.navy.mil  
        % Website: 
        % June 2022; Last revision: 8 July 2022
        %==========================================================================   
            x = data(:,1);
            yR = data(:,2);       
            
            if numel(data) == 2
                yI = zeros(size(x));
            else
                yI = data(:,3);
            end

            % Count the number of active parameters for the fit
            kcount = 0;
            for i = 1:numel(aParams)
                if aParams(i) == 1
                    kcount = kcount + 1;
                end
            end

            nActiveParameters = kcount; %length(b);    
            nVertices = nActiveParameters+1;   
            maxIters = 2000;
            epsTol = 1.0e-12; %1.0e-6; %
            alphaMag = 2.0;
            betaCon = 0.25;
            gammaMag = 2.0;
            maxMSE = 1.0e8;
        
            widthRegion = zeros(1,nActiveParameters);
            bFit.coefficients = zeros(size(b));
            for i = 1:length(b)
                if aParams(i) == 0
                    bFit.coefficients(i) = b(i);
                end
            end
            bFit.converged = 0;
            
            Jmat = zeros(nVertices,1);
            lowMeanVertex = zeros(nActiveParameters,1);
            trackJMax = zeros(maxIters,1);
            trackJMin = zeros(maxIters,1);

            for iteration = 1:maxIters
        
                if iteration == 1
                    aStartingSimplex = obj.buildSimplexWidth(b,aParams);
                    % if ~isempty(obj.cLimits)                        
                    %     aStartingSimplex = obj.buildSimplexConstraint(b,aParams,obj.cLimits,obj.fitFun,x,yR,yI,nActiveParameters,nVertices);
                    % else
                    %     aStartingSimplex = obj.buildSimplexWidth(b,aParams);
                    % end
                    % Calculate the objective functions for each vertex
                    jdx = 1;
                    for i = 1:nVertices
                        aVertex = aStartingSimplex(i,:);
                        [~,zr,zi,~,gR,~] = obj.fitFun(obj.includeAllValues(aVertex,b,aParams),x);
                        if gR == true
                            Jmat(jdx) =  objectiveFunction(yR, yI, zr, zi, nActiveParameters);    
                        else
                            Jmat(jdx) = maxMSE;
                        end
                        jdx = jdx + 1;
                    end                    
                else
                    aStartingSimplex = aSimplex;
                    Jmat = Jsort;
                end
               
                % Find the Jmin and Jmax points
                [Jsort,I] = sort(Jmat);
                trackJMin(iteration,1) = Jsort(1); %Jmin;
                aSimplex = aStartingSimplex(I,:);

                aSimplexLow = aSimplex(1:nVertices-1,:);
                lowMeanVertex = mean(aSimplexLow); %,2 <- This was from the old way of arranging the matrix of vertices.  Now, they are in rows, so we need the mean of each column...

                vReflect = obj.reflectVertex(aSimplex(nVertices,:),lowMeanVertex,b,alphaMag,aParams,obj.cLimits);
        
                % Calculate the objective functions for the reflected point
                aVertex = obj.includeAllValues(vReflect,b,aParams);
                [~,zr,zi,~,gR,fParams] = obj.fitFun(aVertex,x);
                if gR == true
                    JRef =  objectiveFunction(yR, yI, zr, zi, nActiveParameters);
                else
                    JRef = maxMSE;
                end

                % Now, determine what to do with the reflected vertex
                if JRef < Jsort(1)
                    % Try extending the reflected vertex
                    vExtend = obj.extendVertex(aSimplex(nVertices,:),lowMeanVertex,b,gammaMag,aParams,obj.cLimits);        
                    % Calculate the objective functions for the expanded point
                    aVertex = obj.includeAllValues(vExtend,b,aParams);
                    [~,zr,zi,~, gR, fParams] = obj.fitFun(aVertex,x);
                    if gR == true
                        JExt =  objectiveFunction(yR, yI, zr, zi, nActiveParameters);

                        if JExt <= JRef %JExt <= Jmax && 
                            % Accept the extended point                     
                            aSimplex(nVertices,:) = vExtend;
                            Jsort(nVertices) = JExt;                   
                        else
                            % Accept the reflected point                    
                            aSimplex(nVertices,:) = vReflect;
                            Jsort(nVertices) = JRef;                 
                        end
                    else
                        Jsort(nVertices) = maxMSE;
                    end                 
                elseif JRef < Jsort(nVertices-1) && JRef >= Jsort(1)
                    % Accept the reflected point           
                    aSimplex(nVertices,:) = vReflect;
                    Jsort(nVertices) = JRef;                   
                elseif JRef >= Jsort(end-1)
                    % Try contracting the vertex with the worst objective function
                    % value
                    if JRef < Jsort(nVertices)
                        vContract = obj.contractOutsideVertex(vReflect,lowMeanVertex,b,betaCon,aParams,obj.cLimits);
                        aVertex = obj.includeAllValues(vContract,b,aParams);
                        [~,zr,zi,~, gR, ~] = obj.fitFun(aVertex,x);
                        if gR == true
                            JCon =  objectiveFunction(yR, yI, zr, zi, nActiveParameters);
                            if JCon < JRef
                                aSimplex(nVertices,:) = vContract;
                                Jsort(nVertices) = JCon;                                
                            else

                                % Shrink all vertices                       
                                for i = 2:nVertices
                                    for j = 1:nActiveParameters
                                        dist = aSimplex(i,j) - aSimplex(1,j);
                                        tPoint = aSimplex(1,j) + betaCon*dist;
                                        aSimplex(i,j) = tPoint;
                                    end
                                end
                               
                                % disp('Shrink!')

                            end
                        else
                            Jsort(nVertices) = maxMSE;
                        end
                    elseif JRef >= Jsort(nVertices)
                        vContract = obj.contractInsideVertex(aSimplex(nVertices,:),lowMeanVertex,b,betaCon,aParams,obj.cLimits);
                        aVertex = obj.includeAllValues(vContract,b,aParams);
                        [~,zr,zi,~, gR, ~] = obj.fitFun(aVertex,x);
                        if gR == true
                            JCon =  objectiveFunction(yR, yI, zr, zi, nActiveParameters);
                            if JCon < Jsort(nVertices)
                                aSimplex(nVertices,:) = vContract;
                                Jsort(nVertices) = JCon;                                
                            else

                                % Shrink all vertices                       
                                for i = 2:nVertices
                                    for j = 1:nActiveParameters
                                        dist = aSimplex(i,j) - aSimplex(1,j);
                                        tPoint = aSimplex(1,j) + betaCon*dist;
                                        aSimplex(i,j) = tPoint;
                                    end
                                end
                                                               
                                % disp('Shrink!')
                            end
                        else
                            Jsort(nVertices) = maxMSE;
                        end                        

                    end

                end
                trackJMax(iteration,1) = max(Jsort); %(nVertices);
          
                if iteration > 1
                    JmaxCheck = abs(Jsort(nVertices) - trackJMax(iteration-1,1))/Jsort(nVertices);
                    JminCheck = abs(Jsort(1) - trackJMin(iteration-1,1))/Jsort(1);
        
                    if JmaxCheck < epsTol && JminCheck < epsTol %
                        lowMeanVertex = aSimplex(1,:); %mean(aSimplex); %,2
                        aVertex = obj.includeAllValues(lowMeanVertex,b,aParams);
                        bFit.coefficients = aVertex; 
                        bFit.converged = 1;   
                        MSE = Jsort(1);
                        obj.iters = iteration;
                        break;
                    end
                    if (Jsort(1) == 0)
                        lowMeanVertex = aSimplex(1,:); %mean(aSimplex); %,2
                        aVertex = obj.includeAllValues(lowMeanVertex,b,aParams);
                        bFit.coefficients = aVertex; 
                        bFit.converged = 1;   
                        MSE = Jsort(1);
                        obj.iters = iteration;
                        break;                        
                    end
                end        
            end
        
            if iteration == maxIters
                lowMeanVertex = aSimplex(1,:); %mean(aSimplex); %,2
                nPs = logical(aParams);
                aVertex = lowMeanVertex(nPs);
                % aVertex = b(~nPs);
                % aVertex = obj.includeAllValues(avgVertex,b,aParams);
                bFit.coefficients = aVertex; 
                bFit.converged = 1;    
                MSE = maxMSE;
                obj.iters = maxIters;
            end

            iters2 = trackJMax > 0;
            iterations = 1:1:maxIters;
            iterations = iterations(iters2);

            iters3 = trackJMin > 0;
            iterations2 = 1:1:maxIters;
            iterations2 = iterations2(iters3);

            % figure(850)
            % hold on
            % plot(iterations,trackJMax(trackJMax > 0),'-bo')
            % plot(iterations2,trackJMin(trackJMin > 0),'-r+')
            % ax = gca;
            % ax.YScale = 'log';
            % hold off
        end
    end
    
    methods (Access=private)      
        function a = includeAllValues(~,aVertex,b,parametersArray)    
        % includeAllValues - This function includes the stationary parameters into a vertex
            a = zeros(size(b));
            nPs = logical(parametersArray)';
            a(nPs) = aVertex;
            a(~nPs) = b(~nPs);
        end
        function a = buildSimplexWidth(~,b,aP)
        %buildSimplexWidth - Builds an initial simplext by randomly selecting
        %vertices within the confines of the db vector limits
            limits = 0.1;
            nParameters = numel(b);
            nVertices = nParameters + 1;
            
            a = zeros(nVertices,nParameters);
            bsimp = rand(nVertices,nParameters);
            a(1,:) = b;
            for i = 2:nVertices
                for j = 1:nParameters
                    lowerLimit = (1.0-limits) * b(j); %using the MATLAB approcah algorithm for fminsearch
                    upperLimit = (1.0+limits) * b(j); %using the MATLAB approcah algorithm for fminsearch     
                    if aP(j) == 1
                        if bsimp(i,j) >= 0.5
                             a(i,j) = lowerLimit;
                        else
                            a(i,j) = upperLimit;
                        end                 
                    end                    
                end
            end
        end
        function a = buildSimplexConstraint(obj,b,aP,cL,fun,x,yR,yI,nActiveParameters,nRandVertices)
        %buildSimplexConstraint - builds an initial simplex by randomly selecting
        %vertex points between the constrained limits of the active b parameters.
            nRandVertices = 1000; % %
            nVertices = numel(b) + 1;
            maxMSE = 1.0e8; 
            limits = 0.1;

            nParameters = length(b);
            a = zeros(nVertices,nActiveParameters);

            rPoints = rand(nRandVertices,nActiveParameters);  %,nVertices
            bTemp = zeros(nRandVertices,nActiveParameters);
            J = ones(1,nRandVertices).*maxMSE;            
        
            aP2 = logical(aP);             
            lowerLimit = (1.0-limits).*b; %using the MATLAB approach algorithm for fminsearch
            upperLimit = (1.0+limits).*b; %using the MATLAB approach algorithm for fminsearch

            % First candidate vertex is our initial guess
            bTemp(1,:) = b;
            % aCheck = ~abs(b) > 0;            
            % for j = 1:numel(aCheck)
            %     if aCheck(j) == 1
            %         b(j) = 2.5e-3;
            %         bTemp(j,1) = 2.5e-3; %using the MATLAB approach algorithm for fminsearch
            %     end
            % end
            % j = 1;
            for i = 2:nRandVertices
                for j = 1:nActiveParameters
                    bTemp(i,j) = rPoints(i,j) * (cL(j,2) - cL(j,1));
                % bTemp(:,i) = b;
                % if rPoints(j) >= 0.5
                %     bTemp(j,i) = upperLimit(j); %using the MATLAB approach algorithm for fminsearch
                % else
                %     bTemp(j,i) = lowerLimit(j); %using the MATLAB approach algorithm for fminsearch
                % end
                end
                % j = j + 1;
            end

            % Calculate the objective functions for each candidate vertex
            for i = 1:nRandVertices
                [~,zr,zi,~,gR,~] = fun(obj.includeAllValues(bTemp(i,:),b,aP),x);
                if gR == true
                    J(1,i) = objectiveFunction(yR, yI, zr, zi, nActiveParameters);    
                end
            end
            % Sort the objective function array...
            [~,I] = sort(J);
            a(1:nVertices,:) = bTemp(I(1:nVertices),:);
        end
        
        function a = reflectVertex(~,aVertex,avgVertex,~,alpha,aP,cL)
            %reflectVertex - This function reflects a vertex across the average a distance alpha
            % a = zeros(size(aVertex));
            aPs = logical(aP);
            lowerLimit = cL(aPs,1);
            upperLimit = cL(aPs,2);

            a = (alpha.*avgVertex) - aVertex;
            gC = a > upperLimit';
            a(gC) = upperLimit(gC);

            lC = a < 0; %lowerLimit;
            a(lC) = lowerLimit(lC);
        end    
        function a = contractOutsideVertex(~,aVertex,avgVertex,~,beta,aP,cL)
            %contractVertex - This function contracts a vertex toward the average a distance beta
            % a = zeros(size(aVertex));
            aPs = logical(aP);
            lowerLimit = cL(aPs,1);
            upperLimit = cL(aPs,2);

            % distNum = aVertex - avgVertex;
            distNum = avgVertex - aVertex;

            a = avgVertex + beta.*distNum;
            gC = a > upperLimit';
            a(gC) = upperLimit(gC);

            lC = a < lowerLimit';
            a(lC) = lowerLimit(lC); 
        end 
        function a = contractInsideVertex(~,aVertex,avgVertex,~,beta,aP,cL)
            %contractVertex - This function contracts a vertex toward the average a distance beta
            % a = zeros(size(aVertex));
            aPs = logical(aP);
            lowerLimit = cL(aPs,1);
            upperLimit = cL(aPs,2);

            distNum = aVertex - avgVertex;
            % distNum = avgVertex - aVertex;

            a = avgVertex + beta.*distNum;
            gC = a > upperLimit';
            a(gC) = upperLimit(gC);

            lC = a < lowerLimit';
            a(lC) = lowerLimit(lC); 
        end        
        function a = extendVertex(~,aVertex,avgVertex,~,gamma,aP,cL)
            %extendVertex - This function reflects a vertex across the average a distance alpha
            % a = zeros(size(aVertex));
            aPs = logical(aP);
            lowerLimit = cL(aPs,1);
            upperLimit = cL(aPs,2);

            distNum = avgVertex - aVertex;            
            a = avgVertex + gamma.*distNum;
            gC = a > upperLimit';
            a(gC) = upperLimit(gC);

            lC = a < lowerLimit';
            a(lC) = lowerLimit(lC); 
        end        
    end
end