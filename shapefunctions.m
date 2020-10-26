function [ShapeFuncs] = shapefunctions(elements,Nodedata,phi)
% This function retuns the shape function constans B and C found in the
% equation u = A + Bx + Cy. When taking the grad of this function you
% retrive grad(u) = B xhat + C yhat. ShapeFuncs is an array where the
% columns ar B and C and the rows corespond with each element. 
% THE GRAD OF SHAPE FUNCTIONS ARE CONSTANT THROUGHOUT THE WHOLE ELEMENT.

% Initalise the shape fuction array
numele = size(elements,1);
ShapeFuncs = zeros(numele,2);
% Loop through all elements and use cramers rule to find B and C
% (we do not care about A).
for e = 1:numele
    
    % Which nodes are in an element
    Nodes = elements(e,:);
    NodesLoc = Nodedata(Nodes,:);
    Nodex = NodesLoc(:,1);
    Nodey = NodesLoc(:,2);
    phiNodes = phi(Nodes);
    % Calculate the determinate of the divisor

    A = Nodex(1) * (Nodey(2) - Nodey(3)) + Nodex(2) * (Nodey(3) - Nodey(1)) ...
        + Nodex(3) * (Nodey(1) - Nodey(2));
    
    % using cramers rule to find B and C
    B = phiNodes(1) * (Nodey(2) - Nodey(3)) + phiNodes(2) * (Nodey(3) - Nodey(1)) ...
        + phiNodes(3) * (Nodey(1) - Nodey(2));
    B = B / A;
    
    C = phiNodes(1) * (Nodex(3) - Nodex(2)) + phiNodes(2) * (Nodex(1) - Nodex(3)) ...
        + phiNodes(3) * (Nodex(2) - Nodex(1));
    C = C / A;
    
    ShapeFuncs(e,:) = [B C];
    
end