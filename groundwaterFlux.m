function [F] = groundwaterFlux(h,hOld,BoundaryElements,nodeBLengths,nodeCLengths,...
    CVAreas,elements,Kxx,Kzz,nodeLocations,nodeClass,alphaPar,npar,mpar,psiRes,psiSat,...
    SCVAreas,SCVFaceLengths,SCVNx,SCVNy,numNodes,numElements,linenodes,...
    shapeFuncsOld,psiPOld,kPOld,theta,deltat,curTime,oldTime,rainfall,pumpingRatesPerNode)


% Create the F function for the Newton's solver


flux = zeros(numNodes,1);
F = flux;
H = nodeLocations(:,2) + h;
shapeFuncs = shapefunctions(elements,nodeLocations,H);
nodeOrder = [2 3 1];
% loop over every element
[psiP, kP] = updatePsiAndk(h,alphaPar,npar,mpar,elements,numElements,...
    numNodes,SCVAreas,CVAreas,psiRes,psiSat);

for e = 1:numElements
    
    element = elements(e,:);
    
    % Loop over the 3 nodes in the element
    for n = 1:3
        
        % Isolate node and adjacanct node numbers
        node = element(n);
        otherNode = element(nodeOrder(n));
         % Arithmetic Averaging
            
            kNode = 0.5 * (kP(node) + kP(otherNode));
            kNodeOld = 0.5 * (kPOld(node) + kPOld(otherNode));
       
        
        
        implicit = theta * (SCVFaceLengths(e,n)...
            * ((SCVNx(e,n) * kNode * Kxx(e) * shapeFuncs(e,1))...
            + (SCVNy(e,n) * kNode * Kzz(e) * shapeFuncs(e,2))));
        
        explicit = (1 - theta) * (SCVFaceLengths(e,n)...
            * ((SCVNx(e,n) * kNodeOld * Kxx(e) * shapeFuncsOld(e,1))...
            + (SCVNy(e,n) * kNodeOld * Kzz(e) * shapeFuncsOld(e,2))));
        
        flux(node) = flux(node) + implicit + explicit;
        flux(otherNode) = flux(otherNode) - (implicit + explicit);
        
        
        
    end
end

FVMConst = deltat ./ CVAreas;
F = psiP - psiPOld - FVMConst .* flux;


Rain  = rainfall;
Rb = 75;
Hr = 85;
deltaHxR = 40;
Kr = 0.05;
% Rainfall bc
for node = 1:numNodes
    
    if H(node) < 100
        switch nodeClass(node) % Rainfall
            case 4 % Top boundary
                
                F(node) = F(node) - FVMConst(node) * nodeBLengths(node) * (Rain);
                
            case 5 % River boundary
                riverFlux = 0;
                if H(node) <= 100
                    riverFlux = Kr * (Hr - H(node)) / deltaHxR;
                    if H(node) < Rb
                        riverFlux = Kr * (Hr - Rb) / deltaHxR;
                    end
                end
                
                F(node) = F(node) - FVMConst(node) * nodeBLengths(node) * riverFlux;
                
            case 9 % top right corner
                
                F(node) = F(node) - FVMConst(node) * nodeCLengths(4,1)  * (Rain);
                
            case 10 % top left corner
                
                riverFlux = 0;
                if H(node) <= 100
                    riverFlux = Kr * (Hr - H(node)) / deltaHxR;
                    if H(node) < Rb
                        riverFlux = Kr * (Hr - Rb) / deltaHxR;
                    end
                end
                
                F(node) = F(node) - FVMConst(node) *...
                    (nodeCLengths(6,2)  * (Rain) + nodeCLengths(6,1) * riverFlux);
                
            case 11 % Bottom Boundary
                
                riverFlux = 0;
                if H(node) <= 100
                    riverFlux = Kr * (Hr - H(node)) / deltaHxR;
                    if H(node) < Rb
                        riverFlux = Kr * (Hr - Rb) / deltaHxR;
                    end
                end
                F(node) = F(node) - FVMConst(node) * nodeCLengths(7,2) * riverFlux;
                
            case 12 % Town Bore hole
                if curTime > 365 * 1
                    F(node) = F(node) + deltat * pumpingRatesPerNode;
                end
                
            case 13 % repairng vegetation zone
                if psiP(node) > 0.165
                   F(node) = F(node) + deltat *  0.05 * Rain *(( nodeLocations(node,2) - 85)^2)/225;
                end
                
                if nodeLocations(node,2) == 100
                    F(node) = F(node) - FVMConst(node) * nodeBLengths(node) * (Rain);
                end
                
                case 14 % Crop zone
                if psiP(node) > 0.165
                   F(node) = F(node) + deltat *  0.025 * Rain *(( nodeLocations(node,2) - 95)^2)/25;
                end
                
                if nodeLocations(node,2) == 100
                    F(node) = F(node) - FVMConst(node) * nodeBLengths(node) * (Rain);
                end
                
                case 15 % Pine plantation
                if psiP(node) > 0.025
                   F(node) = F(node) + deltat *  0.0035 * Rain *(( nodeLocations(node,2) - 90)^2)/100;
                end
                
                if nodeLocations(node,2) == 100
                    F(node) = F(node) - FVMConst(node) * nodeBLengths(node) * (Rain);
                end
            otherwise
                
        end
    end
end








