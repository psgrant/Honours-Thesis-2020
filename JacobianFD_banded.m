function J = JacobianFD_banded(Ffunc, x, Fx0,lower,upper)
% JACOBIANFD Computes the Finite Difference Jacobian of Ffunc evaluated at x
FUllbandwidth = lower + upper + 1;

N = size(x,1);
J = (zeros(N,N));


for j = 1:FUllbandwidth
    s = zeros(N,1);
    for i = j:FUllbandwidth:N
        s(i) = s(i) + 1;
    end
    
    if norm(x,2) == 0
        h = sqrt(eps)/norm(s);
    else
        h = sqrt(eps)*norm(x,2)/norm(s);
    end
    
    J_band = (Ffunc(x + h*s) - Fx0)/h;
    
    for k = j:FUllbandwidth:N
        
        row_start_pos = max(k-lower,1);
        row_end_pos = min(k+upper,N);
        
        J(row_start_pos:row_end_pos,k) = J_band(row_start_pos:row_end_pos);
        
        
    end
    
end
[ii,jj,rr] = find(J);
J = sparse(ii,jj,rr);

end

