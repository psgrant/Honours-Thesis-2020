function [x,m] = pGMRESv02(A,b,x0,L,U,p,tol,MaxKrylov,diagnostics)
%Calculates the General Minimal Residual Method (GMRES) for a given matrix
%A and solution vector b.
%Input: Matrix A, vector b
%output: solution vector x and number of loops m

%Initialise
N = size(A,1);
H = zeros(MaxKrylov+1,MaxKrylov);
V = zeros(N,MaxKrylov+1);

r = b - A*x0;
beta = norm(r,2);
V(:,1) = r/beta;
m=0;
rnorm=inf;

% [L,U,p] = lu(M);

while rnorm > beta*tol && m <= MaxKrylov
   
    m=m+1;
    %Arnoldi (Modified Gram-Schmidt)
    V(:,m+1) = A*(U\(L\(p*V(:,m))));
    
    %V1(:,m+1) = A*(M\V1(:,m));
    for j = 1:m
        H(j,m) = V(:,j)'*V(:,m+1);
        V(:,m+1) = V(:,m+1) - H(j,m)*V(:,j);
    end
    H(m+1,m) = norm(V(:,m+1),2);
    
    % Check for breakdown
    if abs(H(m+1,m)) < 1e-14
%         fprintf('Invariant Krylov Subspace detected at m=%g\n',m);
        y = H(1:m,1:m) \ ([beta; zeros(m-1,1)]); % Invariant space detected
        break;
    else
        V(:,m+1) = V(:,m+1)/H(m+1,m);
    end
    % Solve small m dimensional least squares problem for y
    rhs= [beta; zeros(m,1)];
    y = H(1:m+1,1:m) \ rhs;
    % determine residual norm
    rnorm = norm(rhs-H(1:m+1,1:m)*y);
       
    if diagnostics, fprintf('m=%g ||r_m||=%g tol=%g\n',m,rnorm,beta*tol); end
    
end

%% Compute approximate solution
x = x0 + (U\(L\(p*((V(:,1:m)*y)))));

end


