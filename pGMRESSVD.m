function [x,m] = pGMRESSVD(A,b,x0,M,tol,MaxKrylov,diagnostics)
%Calculate the General Minimal Residual Method (GMRES) for a given
%A and solution vector b.
%Input: Matrix A, vector b
%output: Solution vector x and number of loops m

%Initialise
N = size(A,1);
H = zeros(MaxKrylov+1, MaxKrylov);
V = zeros(N, MaxKrylov+1);
% y = zeros(1,size(H,2));

r = b - A*x0;
beta = norm(r,2);
V(:,1) = r/beta;
m = 0;
rnorm = inf;
% x0 = zeros(length(b));

% L = eye(size(M)) + tril(M,-1);
% U = triu(M);

while rnorm > beta*tol && m <= MaxKrylov
    
    m = m+1;
    %Arnoldi (Modified Gram-Schmidt)
    V(:,m+1) = A*(M\V(:,m));
    for j = 1:m
        H(j,m) = V(:,j)'*V(:,m+1);
        V(:,m+1) = V(:,m+1) - H(j,m)*V(:,j);
    end
    H(m+1,m) = norm(V(:,m+1),2);
    
    %Check for breakdown
    if abs(H(m+1,m)) < 1e-14
        %         fprintf('Invariant Krylov Subspace detected at m = %g \n',m);
        y = H(1:m,1:m) \ ([beta; zeros(m-1,1)]); %Invariant space detected
        break;
    else
        V(:,m+1) = V(:,m+1)/H(m+1,m);
    end
    
    %Solve small m dimensional least squares problem for y

    
    rhs= [beta; zeros(m-1,1)];
    [U,S,V1] = svd(H(1:m+1,1:m));
    y = V1(1:m,1:m) * diag((diag(1./S(1:m,1:m)))) * (U(1:m,1:m)') * rhs;
    rnorm = beta*abs(U(1,m+1));
    
%     
%             rhs= [beta; zeros(m,1)];
%         y = H(1:m+1,1:m) \ rhs;
% %         determine residual norm
%         rnorm = norm(rhs-H(1:m+1,1:m)*y);
    
    %    if diagnostics, fprintf('m=%g ||r_m||=%g tol=%g\n',m,rnorm,beta*tol); end
    
    
end

%Compute approximate solution
x = x0 + M\(V(:,1:m)*y);

end




