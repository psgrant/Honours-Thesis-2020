% initaliation
clear all
clc
load mesh.mat
load turbo_colormap.mat
mex Diff_Advec.c


Dx = 0.1;
Dy = 0.1;
timescale = 0:0.01:1;
xc = 0;
yc = 0;
vx = 1;
vy = 1;
deltat = timescale(2) - timescale(1);
theta = 0.5;
x = Nodedata(:,2);
y = Nodedata(:,3);
[X,Y] = meshgrid(x,y);
T = elements(:,2:4);

% Inital condition
phi = exp(-((Nodedata(:,2)-xc).^2 / (Dx)) - ((Nodedata(:,3)-yc).^2 / (Dy)));
phiold = phi;
phiSave = zeros(size(phi,1),length(timescale));
phiSaveExact = phiSave;
z = connectionMatrix(RCM,RCM);
[lower upper] = bandwidth(z);
MatrixBandwidth = lower + upper + 1;
elements = elements(:,2:4);
FluxConst = 0;
exitflag = [0 0 0];
mex_file = 0;
breakout = 0;
simtime = zeros(8,5);
simCount = 1;
%%
for sim = 1:40
    sim
    switch sim
        % Select the differnt preconditoners
        case 1
            preCondSelector = 1;
        case 6
            preCondSelector = 2;
        case 11
            preCondSelector = 3;
        case 16
            preCondSelector = 4;
        case 21
            preCondSelector = 5;
        case 26
            preCondSelector = 6;
        case 31
            preCondSelector = 7;
        case 36
            preCondSelector = 8;
        otherwise
            
    end
    
    % Inital condition
    phi = exp(-((Nodedata(:,2)-xc).^2 / (Dx)) - ((Nodedata(:,3)-yc).^2 / (Dy)));
    phiold = phi;
    phiSave = zeros(size(phi,1),length(timescale));
    phiSaveExact = phiSave;
    
    
    % For each timestep
    tic
    for t = 1:length(timescale)
        % return the current time
        time = timescale(t);
        % save the volumetric average for each timestep
        phiSave(:,t) = phi.*CVAreas;
        x = Nodedata(:,2);
        y = Nodedata(:,3);
        phiExact = 1/(4*time+1) * exp(-(x-vx*time-xc).^2./(Dx*(4*time+1))-(y-vy*time-yc).^2./(Dy*(4*time+1)));
        phiSaveExact(:,t) = phiExact .* CVAreas;
        
%         % Uncomment for live plot
        %         subplot(121)
        %         TriPlot1 = trimesh(T,x,y,phi,'facecolor',[0 0 0]);
        %         xlabel('x');
        %         ylabel('y')
        %         shading interp;
        %         colormap(turbo)
        %         axis square
        %         title(['\fontsize{14}Numerical t = ', num2str(time)])
        %         caxis([0 1]);
        %         zlim([0 1])
        %         colorbar
        %         view(2)
        %
        %         subplot(122)
        %         TriPlot1 = trimesh(T,x,y,phiExact,'facecolor',[0 0 0]);
        %         xlabel('x');
        %         ylabel('y')
        %         shading interp;
        %         colormap(turbo)
        %         axis square
        %         title(['\fontsize{14}Exact t = ', num2str(time)])
        %         caxis([0 1]);
        %         zlim([0 1])
        %         colorbar
        %         view(2)
        %         drawnow
        %
        
        parameters = [Dx,Dy,deltat,theta,xc,yc,time,vx,vy];
        %     fprintf(2,' Current Time %f/%f \n',timescale(t),timescale(end))
        
        if mex_file == 1
            % Setup up flux functions in terms of phi
            
            DiffusionFluxes = @(phi) Diff_Advec(phi,phiold,Nodedata,SCVFaceLengths,...
                SCVNx, SCVNy, CVAreas,elements,parameters);
        else
            [ShapeFuncsOld] = shapefunctions(elements,Nodedata(:,2:3),phiold);
            DiffusionFluxes = @(X) Diffusion_Advection_Flux_Dirichlet(X,phiold,Nodedata...
                ,SCVFaceLengths,SCVNx,SCVNy,CVAreas,elements,BoundaryElements,BoundaryLengths...
                ,Dx,Dy,deltat,theta,xc,yc,time,FluxConst,ShapeFuncsOld,vx,vy);
        end
        % Apply newtons iterations
        
        
        [x_out, exitflag] = Full_Newton_Solver_LS(DiffusionFluxes,phi,[0 0 0],lower,upper,preCondSelector);
        
        if exitflag(1) == 1
            breakout = 1;
            break;
            
        end
        
        % Update phi and phi old
        phi = x_out;
        phiold = phi;

        
        
        
    end
    
    
    if simCount == 6
        simCount = 1;
    end
    
    simtime(preCondSelector,simCount) = toc;
    
    
    if breakout == 1
        simtime(preCondSelector,simCount) = NaN;
        breakout = 0;
    end
    simCount = simCount + 1;
end
%% print results
simtime
mean(simtime,2)
%%
clf
subplot(2,2,1)

TriPlot1 = trimesh(T,x,y,phi,'facecolor',[0 0 0]);
xlabel('x');
ylabel('y')
shading interp;
colormap(turbo)
axis square
title(['\fontsize{14}Numerical t = ', num2str(time)])
caxis([0 1]);
zlim([0 1])
colorbar
view(2)


subplot(2,2,2)
TriPlot1 = trimesh(T,x,y,phiExact,'facecolor',[0 0 0]);
xlabel('x');
ylabel('y')
shading interp;
colormap(turbo)
axis square
title(['\fontsize{14}Exact t = ', num2str(time)])
caxis([0 1]);
zlim([0 1])
colorbar
view(2)


subplot(2,2,[3 4])
phiMean = sum(phiSave,1) / sum(CVAreas);
phiExactMean = sum(phiSaveExact,1) / sum(CVAreas);

plot(timescale,phiMean,'r','linewidth',2);
hold on
phiExpected = 1 - FluxConst * timescale * 4;
plot(timescale,phiExactMean,'k--','linewidth',2);
hold off
xlabel('\fontsize{18}Time'); ylabel('\fontsize{18}\phi');

legend('\fontsize{12}Numerical', '\fontsize{12}Exact')
% Calcualte error
diffplot = max(phiMean) - min(phiMean);
ylim([0, 0.5])
error = phiMean - phiExactMean;
error = [max(abs(error)); norm(error)]
shg





