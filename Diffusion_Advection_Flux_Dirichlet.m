function [F] =  Diffusion_Advection_Flux_Dirichlet(phi,phiold,Nodedata,SCVFaceLengths,SCVNx,SCVNy,CVAreas,elements,BoundaryElements,BoundaryLengths,Dx,Dy,deltat,theta,xc,yc,time,FluxConst,ShapeFuncsOld,vx,vy)


Fluxout = zeros(size(Nodedata,1),1);

% generate shape functions for phi and oldphi
% Todo: more effecient storage
[ShapeFuncs] = shapefunctions(elements,Nodedata(:,2:3),phi);
NodeOrder = [2 3 1];

% For each element
for e = 1:size(elements,1)
    
    
    % Isolate nodes in element
    Element = elements(e,:);
    
    % For each node in the element
    for n = 1:3
        
        % Isolate node of interest
        Node =  Element(n);
        OtherNode = Element(NodeOrder(n));
        
        
        % Deciding on which phi to use for upwinding
       
        
        phiNode = 0.5 * (phi(Node) + phi(OtherNode));
        phioldNode = 0.5 * (phiold(Node) + phiold(OtherNode));
        % Implicit part of the flux (theta = 1)
        % Add the flux for each node in each element
        Implicit = theta * (SCVNx(e,n) * SCVFaceLengths(e,n)...
            * (Dx * ShapeFuncs(e,1) - vx * (phiNode))...
            + (SCVNy(e,n) * SCVFaceLengths(e,n)...
            * (Dy * ShapeFuncs(e,2) - vy *  (phiNode))));
        
        % Implicit part of the flux (theta = 0)
        % Add the flux for each node in each element
        Explicit = (1-theta) * (SCVNx(e,n) * SCVFaceLengths(e,n)...
            * (Dx * ShapeFuncsOld(e,1) - vx * (phioldNode))...
            + (SCVNy(e,n) * SCVFaceLengths(e,n)...
            * (Dy * ShapeFuncsOld(e,2) - vy * (phioldNode))));
        
        
        % Add them to the flux term
        Fluxout(Node) = Fluxout(Node) + (Implicit + Explicit);
        Fluxout(OtherNode) = Fluxout(OtherNode) - (Implicit + Explicit);
        
        
        
        
        
        
        
    end
    
    
end

% Calculate constant for that node
FVMConst = deltat ./ CVAreas;

F = phi - phiold - FVMConst .* Fluxout;
% Loop through all the boundary nodes and apply a no flux boundary
% condition

% Loop through each node and apply a Dirichlet boundary condition

% For each node
for Node = 1:size(Nodedata,1)
    
    %     if the Node is a corner or boundary node
    if Nodedata(Node,1) == 1 || Nodedata(Node,1) == 2
        
        % Retrive x, and y positions of the node
        x = Nodedata(Node,2);
        y = Nodedata(Node,3);
        
        % Calculate time dependent Dirichlet BC
        F(Node) = phi(Node) - 1/(4*time+1) * exp(-(x-vx*time-xc)^2/(Dx*(4*time+1))-(y-vy*time-yc)^2/(Dy*(4*time+1)));
        
    end
    
end







end


