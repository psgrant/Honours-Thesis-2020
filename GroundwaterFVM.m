% Nonlinear groundwater flow problem in MATLAB
clear all
clc
mex GroundWaterFlux_MEX.c
% mex generatePsiPAndkP_MEX.c

%%

load mesh.mat % Load Mesh
load turbo_colormap.mat % Load custom colormap
simtim = zeros(5,1);

elements = elements(:,2:4);
nodeClass = Nodedata(:,1);
nodeLocations = Nodedata(:,2:3);
numNodes = size(nodeClass,1);
numElements = size(elements,1);
[lower upper] = bandwidth(connectionMatrix(RCM,RCM));

% Assign which node is in which material
% 1 - Alluvim, 2 - Confining layer, 3 - Coal, 4 - Volcanics
simtime = zeros(8,1);
simCount = 1;
for sim = 1:1
    % initalisation
    switch sim
        
        case 1
            preCondSelector = 1;
        case 2
            preCondSelector = 2;
        case 3
            preCondSelector = 3;
        case 4
            preCondSelector = 4;
        case 5
            preCondSelector = 5;
        case 6
            preCondSelector = 6;
        case 7
            preCondSelector = 7;
        case 8
            preCondSelector = 8;
        otherwise
            
    end
    
    elementMaterials = zeros(size(nodeClass));
    
    for element = 1:numElements
        
        
        x = elecentx(element);
        z = elecenty(element);
        % Alluvium (0 <= x <= 50, 40 <= z <= 100)
        if (((x <= 50) && (40 <= z)) || ((50 <= x) && (x <= 300) && (50 <= z)))
            
            elementMaterials(element) = 1;
            
        elseif (((50 <= x) && (x <= 350)) && ((40 <= z) && (z <= 50)))
            % Confining
            elementMaterials(element) = 2;
            
            
        elseif ((z <= 40) || ((350 <= x) && ((40 <= z) && (z <= 50))))
            % Coal
            elementMaterials(element) = 3;
            
        elseif ((300 <= x) && ((50 <= z)))
            % Volcanics
            elementMaterials(element) = 4;
            
        else
            error('Element center not defined')
        end
        
    end
    
    
    % Plot elements in their respective materials
    % plotElements(elecentx,elecenty,elementMaterials,elements,Nodedata)
    
    % Hydraulic conductivities for each element
    
    KXX = [3.185 0.092 4.095 4.95];
    KZZ = [0.91 0.016 1.17 0.21];
    NPAR = [1.51 1.3964 2.239 2];
    ALPHA = [1.43 1.04 2.8 2.5];
    RES = [0.01 0.106 0.01 0.01];
    SAT = [0.33 0.486 0.1 0.05];
    
    
    Rb = 75;
    Hr = 85;
    deltaHxR = 40;
    Kr = 0.05;
    
    
    homogenousAquifer = 0;
    if homogenousAquifer == 1
        
        KXX = mean(KXX) .* ones(1,4);
        KZZ = mean(KZZ) .* ones(1,4);
        NPAR = mean(NPAR) .* ones(1,4);
        ALPHA = mean(ALPHA) .* ones(1,4);
        RES = mean(RES) .* ones(1,4);
        SAT = mean(SAT) .* ones(1,4);
        
    end
    
    
    
    
    Kxx = KXX(1) * (elementMaterials == 1) + KXX(2) * (elementMaterials == 2)...
        + KXX(3) * (elementMaterials == 3) + KXX(4) * (elementMaterials == 4);
    
    Kzz = KZZ(1) * (elementMaterials == 1) + KZZ(2) * (elementMaterials == 2)...
        + KZZ(3) * (elementMaterials == 3) + KZZ(4) * (elementMaterials == 4);
    
    % N and m paramerters
    npar = NPAR(1) * (elementMaterials == 1) + NPAR(2) * (elementMaterials == 2)...
        + NPAR(3) * (elementMaterials == 3) + NPAR(4) * (elementMaterials == 4);
    
    mpar = 1 - 1./npar;
    
    % Alpha
    alphaPar = ALPHA(1) * (elementMaterials == 1) + ALPHA(2) * (elementMaterials == 2)...
        + ALPHA(3) * (elementMaterials == 3) + ALPHA(4) * (elementMaterials == 4);
    
    % residual and saturated water contents
    psiRes = RES(1) * (elementMaterials == 1) + RES(2) * (elementMaterials == 2)...
        + RES(3) * (elementMaterials == 3) + RES(4) * (elementMaterials == 4);
    
    psiSat = SAT(1) * (elementMaterials == 1) + SAT(2) * (elementMaterials == 2)...
        + SAT(3) * (elementMaterials == 3) + SAT(4) * (elementMaterials == 4);
    
    % domain of x and z
    L1 = 500;
    L2 = 100;
    
    % Inital conditions
    hbot = -5;
    htop = -10;
    
    h = hbot + (htop - hbot).* nodeLocations(:,2) /L2;
    hOld = h;
    h0 = h;
    
    % Used for plotting
    
    [psiP, kP] = updatePsiAndk(h,alphaPar,npar,mpar,elements,numElements,numNodes,SCVAreas,CVAreas,psiRes,psiSat);
    [maxSat, ~] = updatePsiAndk(zeros(size(h)),alphaPar,npar,mpar,elements,numElements,numNodes,SCVAreas,CVAreas,psiRes,psiSat);
    
    % Update node Classes
    nodeClass(nodeClass == 3) = 1; % Interior
    nodeClass(linenodes{1}(:,1)) = 2; % Bottom
    nodeClass(linenodes{2}(:,1)) = 3; % Right z < 50
    nodeClass(linenodes{3}(:,1)) = 3; % Right z > 50
    nodeClass(linenodes{4}(:,1)) = 4; % Top x > 300
    nodeClass(linenodes{5}(:,1)) = 4; % Top x < 300
    nodeClass(cornernodes(5,1)) = 4; % Top x = 300
    nodeClass(cornernodes(7,1)) = 11; % Left z = 77.5
    nodeClass(linenodes{6}(:,1)) = 5; % Left 100 > z > 77.5
    nodeClass(linenodes{7}(:,1)) = 6; % Left 40 < z < 77.5
    nodeClass(linenodes{8}(:,1)) = 6; % Left z < 40
    nodeClass(cornernodes(1,1)) = 7; % Bottom Left
    nodeClass(cornernodes(2,1)) = 8; % Bottom Right
    nodeClass(cornernodes(6,1)) = 10; % Top left
    nodeClass((nodeLocations(:,1) == 50) & ((nodeLocations(:,2) >= 55) & (nodeLocations(:,2) <= 75))) = 12; % pumping
    
    nodeClass(((50 <= nodeLocations(:,1)) & (nodeLocations(:,1) <= 100) & (nodeLocations(:,2) > 85))) = 13;
    nodeClass(((100 < nodeLocations(:,1)) & (nodeLocations(:,1) <= 300) & (nodeLocations(:,2) > 95))) = 14;
    nodeClass(((nodeLocations(:,1) > 300)) & (nodeLocations(:,2) > 90)) = 15; % Evapotranspiration not top boundary
    
    
    nodeClass(cornernodes(4,1)) = 9; % Top right
    
   
    logicCond = (Nodedata(:,2) == 50) & ((Nodedata(:,3) >= 55) & (Nodedata(:,3) <= 75));
    pumpingNodes = nonzeros((1:numNodes)'.*logicCond);
    pumpingElements = [];
    for i = 1:size(pumpingNodes,1)
        pumpingElements = [pumpingElements; nonzeros((1:numElements)' .* sum(elements == pumpingNodes(i),2))];
    end
    pumpingElements = unique(pumpingElements);
    
    
    % Node Boundary Lengths
    nodeBLengths = zeros(size(nodeClass));
    nodeBLengths(linenodes{1}(:,1)) = BoundaryLengths{1}(1) * 2; % Bottom Boundary
    nodeBLengths(linenodes{2}(:,1)) = BoundaryLengths{2}(1) * 2; % Right Boundary
    nodeBLengths(linenodes{3}(:,1)) = BoundaryLengths{3}(1) * 2; % Right Boundary
    nodeBLengths(linenodes{4}(:,1)) = BoundaryLengths{4}(1) * 2; % Top Boundary
    nodeBLengths(linenodes{5}(:,1)) = BoundaryLengths{5}(1) * 2; % Top Boundary
    nodeBLengths(linenodes{6}(:,1)) = BoundaryLengths{6}(1) * 2; % Left Boundary
    nodeBLengths(linenodes{7}(:,1)) = BoundaryLengths{7}(1) * 2; % Left Boundary
    nodeBLengths(linenodes{8}(:,1)) = BoundaryLengths{8}(1) * 2; % left Boundary
    
    % Node Corner lengths
    
    nodeCLengths = zeros(size(cornernodes,1),2);
    nodeCornerOrder = [ 8 1 2 3 4 5 6 7];
    
    for i = 1:size(cornernodes,1)
        
        nodeCLengths(i,1) = BoundaryLengths{i}(1);
        nodeCLengths(i,2) = BoundaryLengths{nodeCornerOrder(i)}(1);
        
    end
    
    % The node between volcanics and alluvium
    nodeBLengths(cornernodes(5,1)) = sum(nodeCLengths(5,:));
    %%
    % Parameters
    upwinding = 0;
    theta = 1;
    dt = 1;
    numYears = 3;
    minDeltaT = 0.0001;
    maxDeltaT = 10;
    upDtAfterConvergences = 8;
    iterCounter = 1;
    curTime = 0;
    
    fluxLimiting = 0;
    updateKzz = 1;
    
    
    
    % chnage in dt array
    dtTimeArray = zeros(1000,1);
    % Current time array
    timeArray = zeros(1000,1);
    
    % Used to bring plot to front only for the first time it is called
    firstPlot = 1;
    newtonFail = 0;
    
    
    
    psiPSave = zeros(size(timeArray));
    rainArray = psiPSave;
    totalArea = sum(CVAreas);
    exitFlag = [0 0 0];
    rainConst = 1 / 1000; % mm/day
    
    % Use mex file 1 = MEX, 0 = MATLAB
    mex_file = 0;
    rainfall = rainConst *  (1 + cos(2*pi*curTime / 365)); % mm/day
    pumpingRatesPerNode = rainConst / length(pumpingNodes) * 0.2; % mm/day/node
    pumpingTurnedOn = 1;
    parameters = [0, fluxLimiting, dt, theta, curTime, rainfall, Rb, Hr, deltaHxR, Kr,pumpingRatesPerNode, pumpingTurnedOn];
    tic
    while curTime < 365 * numYears + minDeltaT
        plotFig = 0; % used if the figure needs to be updated
        fprintf('Percentage of subsimulation %g: %g ====== dt = %g\n',sim, round(100*curTime/1095,2),dt)
        rainfall = rainConst *  (1 + cos(2*pi*curTime / 365)); % mm/day
        % Update params
        parameters(3) = dt;
        parameters(5) = curTime;
        parameters(6) = rainfall;
        %     parameters(11) = rainfall / length(pumpingNodes) * 0.5;
        HOld = hOld + nodeLocations(:,2);
        H = h + nodeLocations(:,2);
        
%         if (curTime > 365 * pumpingTurnedOn && updateKzz == 1)
%             updateKzz = 0;
%             Kzz(pumpingElements) = Kzz(pumpingElements) * 62.5;
%             %         Kxx(pumpingElements) = Kxx(pumpingElements) * 62.5;
%         end
        
        [psiPOld, kPOld] = updatePsiAndk(hOld,alphaPar,npar,mpar,elements,numElements,numNodes,SCVAreas,CVAreas,psiRes,psiSat);
        
        shapeFuncsOld = shapefunctions(elements,nodeLocations,HOld);
        
        
        % Assign function to an anonymous function
        if mex_file == 0
            aquiferFlux = @(h) groundwaterFlux(h,hOld,BoundaryElements,nodeBLengths,nodeCLengths,...
                CVAreas,elements,Kxx,Kzz,nodeLocations,nodeClass,alphaPar,npar,mpar,psiRes,psiSat,...
                SCVAreas,SCVFaceLengths,SCVNx,SCVNy,numNodes,numElements,linenodes,...
                shapeFuncsOld,psiPOld,kPOld,theta,dt,curTime,curTime-dt,rainfall,pumpingRatesPerNode);
        else
            aquiferFlux = @(h) GroundWaterFlux_MEX(h,hOld,nodeLocations,nodeClass...
                ,SCVFaceLengths,SCVNx,SCVNy,CVAreas,elements,Kxx,Kzz,alphaPar...
                ,npar,mpar,psiRes,psiSat,SCVAreas,nodeBLengths,nodeCLengths...
                ,shapeFuncsOld, psiPOld, kPOld,parameters);
            
        end
        % Call newton's solver
        [x_out, exitFlag] = Full_Newton_Solver_LS(aquiferFlux,h,exitFlag,lower,upper,preCondSelector);
        
        
        
        if exitFlag(1) == 1 % Newton Failed to Converge
            %         fprintf(2,'Newton Method has failed at dt = %g,trying dt = %g, Current time: %g \n',dt,dt/2,(curTime))
            exitFlag(2:3) = 0;
            dt = dt / 2; % Half the timestep and try again
            
        elseif (exitFlag(1) == 0) && (exitFlag(2) > upDtAfterConvergences)
            
            %         fprintf('Newton Method effective at dt = %g,trying dt = %g, Current time: %g \n',dt,min(maxDeltaT, 1.2 * dt), (curTime))
            dt = min(maxDeltaT, 1.2 * dt);%Speed up
            % Increse dt and update newtonexitflags
            curTime = curTime + dt;
            timeArray(iterCounter) = curTime;
            plotFig = 1;
            exitFlag(2:3) = 0;
            
            h = x_out;
            hOld = h;
            
            [psiP, ~] = updatePsiAndk(h,alphaPar,npar,mpar,elements,numElements,...
                numNodes,SCVAreas,CVAreas,psiRes,psiSat);
            
            psiPSave(iterCounter) = sum(psiP .* CVAreas) / totalArea;
            rainArray(iterCounter) = rainfall;
            dtTimeArray(iterCounter) = dt;
            iterCounter = iterCounter + 1;
            
        elseif (exitFlag(1) == 0) && (exitFlag(2) <= upDtAfterConvergences)
            % Newton still converges, but hasn't within the last few iters
            %         fprintf('Regular dt: %g, Current time: %g\n',(dt),curTime);
            
            curTime = curTime + dt;
            timeArray(iterCounter) = curTime;
            plotFig = 1;
            
            h = x_out;
            hOld = h;
            
            [psiP, ~] = updatePsiAndk(h,alphaPar,npar,mpar,elements,numElements,...
                numNodes,SCVAreas,CVAreas,psiRes,psiSat);
            psiPSave(iterCounter) = sum(psiP .* CVAreas) / totalArea;
            rainArray(iterCounter) = rainfall;
            dtTimeArray(iterCounter) = dt;
            iterCounter = iterCounter + 1;
            
        end
        
        if dt < minDeltaT
            newtonFail = 1;
            break
            
        end
        
        
        % preallocate more memory in a chunk to increase speed
        if iterCounter == size(psiPSave,1)
            psiPSave(end:end+1000,1) = 0;
            rainArray(end:end+1000,1) = 0;
            dtTimeArray(end:end+1000,1) = 0;
            timeArray(end:(end+1000),1) = 0;
        end
        
%         if mod(iterCounter,3) == 0
%             
%             % Plotting
%             if plotFig == 1
%                 subplot(3,3,[1 2])
%                 trimesh(elements,nodeLocations(:,1),nodeLocations(:,2),psiP,'facecolor',[0 0 0]);
%                 shading interp
%                 colorbar
%                 colormap(turbo)
%                 %                             title(['\fontsize{14}Volumetric \psi,', num2str(round(curTime/365,2)), ' Years'])
%                 xlabel('\fontsize{14}x [m]')
%                 ylabel('\fontsize{14}z [m]')
%                 xlim([0 500])
%                 ylim([0 100])
%                 pbaspect([1 0.35 1])
%                 caxis([0 0.5])
%                 view(2);
%                 
%                 subplot(3,3,[7 8])
%                 trimesh(elements,nodeLocations(:,1),nodeLocations(:,2),h,'facecolor',[0 0 0]);
%                 shading interp
%                 colorbar
%                 colormap(turbo)
%                 xlim([0 500])
%                 ylim([0 100])
%                 %                             title(['\fontsize{14} Pressure head h, Current time = ', num2str(round(curTime/365,2)), ' Years'])
%                 xlabel('\fontsize{14}x [m]')
%                 ylabel('\fontsize{14}z [m]')
%                 pbaspect([1 0.35 1])
%                 caxis([-20 0])
%                 view(2);
%                 
%                 subplot(3,3,[4 5])
%                 
%                 trimesh(elements,nodeLocations(:,1),nodeLocations(:,2),psiP ./ maxSat,'facecolor',[0 0 0]);
%                 shading interp
%                 colorbar
%                 colormap(turbo)
%                 xlim([0 500])
%                 ylim([0 100])
%                 %                             title(['\fontsize{14} Relative Saturation Current time = ', num2str(round(curTime/365,2)), ' Years'])
%                 xlabel('\fontsize{14}x [m]')
%                 ylabel('\fontsize{14}z [m]')
%                 pbaspect([1 0.35 1])
%                 caxis([0 1])
%                 view(2);
%                 
%                 
%                 subplot(3,3,3)
%                 
%                 plot(timeArray(1:iterCounter -1),dtTimeArray(1:iterCounter -1),'b','linewidth',1.5,'markersize',10)
%                 xlabel('\fontsize{14}Time [Days]')
%                 ylabel('\fontsize{14}dt [days]');
%                 %                             title(['\fontsize{14}dt = ', num2str(round(dt,4))])
%                 xlim([0 curTime]);
%                 ylim([0, ceil(max(dtTimeArray))])
%                 
%                 
%                 subplot(3,3,9)
%                 
%                 plot(timeArray(1:iterCounter-1),psiPSave(1:iterCounter-1),'b','linewidth',1.5)
%                 xlabel('\fontsize{14}Time [Days]')
%                 ylabel('\fontsize{14}Volumetic water content');
%                 %                             title('\fontsize{14}Volumetric water content')
%                 ylim([floor(min(psiPSave)*100)/100, ceil(max(psiPSave)*100)/100])
%                 xlim([0 curTime]);
%                 
%                 subplot(3,3,6)
%                 
%                 plot(timeArray(1:iterCounter-1),rainArray((1:iterCounter-1)) * 1000,'b','linewidth',1.5)
%                 xlabel('\fontsize{14}Time [Days]')
%                 ylabel('\fontsize{14}Rainfall [mm/d]');
%                 %                             title('\fontsize{14}Rainfall')
%                 ylim([0 max(rainArray)*1250])
%                 xlim([0 curTime]);
%                 
% %                 
%             end
%             
%             drawnow
%             if firstPlot == 1
%                 shg
%                 firstPlot = 2;
%             end
%             
%         end
%         
%         if simCount == 2
%             simCount = 1;
%         end
        if newtonFail == 1
            simtime(preCondSelector,simCount) = 0;
        else
            simtime(preCondSelector,simCount) = toc;
        end
%         simCount = simCount + 1;
        
    end
end
    %%
    simtime
    mean(simtime,2)/60
    %%
    
    toc
    
    
    %% End Visualisation
    
    
    timeArray = nonzeros(timeArray);
    dtTimeArray = nonzeros(dtTimeArray);
    psiPSave = nonzeros(psiPSave);
    
    
    subplot(3,3,[1 2])
    trimesh(elements,nodeLocations(:,1),nodeLocations(:,2),psiP,'facecolor',[0 0 0]);
    shading interp
    colorbar
    colormap(turbo)
%     title(['\fontsize{14}Volumetric \psi, Current time = ', num2str(round(curTime/365,2)), ' Years'])
    xlabel('\fontsize{14}x [m]')
    ylabel('\fontsize{14}z [m]')
    xlim([0 500])
    ylim([0 100])
    pbaspect([1 0.2 1])
    caxis([0 0.5])
    view(2);
    
    subplot(3,3,[7 8])
    trimesh(elements,nodeLocations(:,1),nodeLocations(:,2),h,'facecolor',[0 0 0]);
    shading interp
    colorbar
    colormap(turbo)
    xlim([0 500])
    ylim([0 100])
%     title(['\fontsize{14} Pressure head h, Current time = ', num2str(round(curTime/365,2)), ' Years'])
    xlabel('\fontsize{14}x [m]')
    ylabel('\fontsize{14}z [m]')
    pbaspect([1 0.2 1])
    caxis([-20 0])
    view(2);
    
    subplot(3,3,[4 5])
    
    trimesh(elements,nodeLocations(:,1),nodeLocations(:,2),psiP ./ maxSat,'facecolor',[0 0 0]);
    shading interp
    colorbar
    colormap(turbo)
    xlim([0 500])
    ylim([0 100])
%     title(['\fontsize{14} Relative Saturation Current time = ', num2str(round(curTime/365,2)), ' Years'])
    xlabel('\fontsize{14}x [m]')
    ylabel('\fontsize{14}z [m]')
    pbaspect([1 0.2 1])
    caxis([0 1])
    view(2);
    
    
    subplot(3,3,3)
    
    plot(timeArray(1:iterCounter -1),dtTimeArray(1:iterCounter -1),'b','linewidth',1.5,'markersize',10)
    xlabel('\fontsize{14}Time [Days]')
    ylabel('\fontsize{14}dt [days]');
%     title('\fontsize{14}Change in dt')
    xlim([0 curTime]);
    ylim([0, ceil(max(dtTimeArray))])
    
    
    subplot(3,3,9)
    
    plot(timeArray(1:iterCounter-1),psiPSave(1:iterCounter-1),'b','linewidth',1.5)
    xlabel('\fontsize{14}Time [Days]')
    ylabel('\fontsize{14}Water content []');
%     title('\fontsize{14}Volumetric water content over time')
    ylim([floor(min(psiPSave)*100)/100, ceil(max(psiPSave)*100)/100])
    xlim([0 curTime]);
    
    subplot(3,3,6)
    
    plot(timeArray(1:iterCounter-1),(1+cos(2*pi/365*timeArray(1:iterCounter-1))),'b','linewidth',1.5)
    xlabel('\fontsize{14}Time [Days]')
    ylabel('\fontsize{14}Rainfall [mm/d]');
%     title('\fontsize{14}Rainfall')
    ylim([0 max(rainArray)*1250])
    xlim([0 curTime]);
    
    shg
    
    
    
    
    
    
    
    
    
