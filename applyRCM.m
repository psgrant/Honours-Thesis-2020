function [Nodedata,RCM,matrixSPY] = applyRCM(numnodes,numelements,Nodedata,elements,plotfig)
%% Applies the RCM algorithm to a supplied mesh

matrixSPY = sparse(zeros(numnodes));
etemp = elements(:,2:4);
for i = 1:numnodes
    elementnumbers = nonzeros(sum((etemp == i),2).*(1:numelements)');
    
    for j = elementnumbers
        nodesconnecting = nonzeros(etemp(j,:).*(etemp(j,:) ~= i));
        matrixSPY(i,nodesconnecting) = 1;
    end
end

RCM = symrcm(matrixSPY);
% Plot interaction matrix before and after RCM ordering
if plotfig == 1
    
    figure
    subplot(121);
    spy(matrixSPY)
    bandwidth(matrixSPY)
    title('Before Ordering')
    subplot(122)
    NodeConnections = matrixSPY(RCM,RCM);
    spy(NodeConnections)
    bandwidth(matrixSPY(RCM,RCM))
    title('After Ordering')
end

% Apply RCM ordering to the nodedata array
Nodedata = Nodedata(RCM,:);
Nodedata(:,1) = (1:numnodes)';
