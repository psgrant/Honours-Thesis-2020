function [psiP, kP] = updatePsiAndk(h,alphaPar,npar,mpar,elements,numElements,numNodes,SCVAreas,CVAreas,psiRes,psiSat)

psiP = zeros(numNodes,1);
kP = psiP;

% Calculating hydrological parameters S, psiP, kP (psi and k at node P)
for element = 1:numElements
    
    for n = 1:3
        node = elements(element,n);
        
        if h(node) < 0
            
            S = (1 + (-alphaPar(element) * h(node)) ^ npar(element)) ^ (-mpar(element));
            
            psiP(node) = psiP(node) + (psiRes(element) + S *(psiSat(element) - psiRes(element))) * SCVAreas(element,n);
            
            kP(node) = kP(node) + ((sqrt(S) * (1 - (1 - S^(1/mpar(element)))^mpar(element))^2))* SCVAreas(element,n);
            
        else
            kP(node) = kP(node) + SCVAreas(element,n);
            psiP(node) = psiP(node) + psiSat(element) * SCVAreas(element,n);
            
        end
        
        
        
        
    end
    
end
psiP = psiP ./ CVAreas;
kP = kP ./ CVAreas;