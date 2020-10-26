function  [x_out, exitflag,ShapeFuncs] = Full_Newton_Solver_LS(F,x,exitflag,lower,upper,preCondSelector)
% F = function handle, x = initial guess
% fprintf('\n Inexact Newton Method:\n');
exitflag(1) = 0;
exitflag(2) = exitflag(2) + 1;
x_return = x;
Resid = F(x);
MaxIters = 15; % Maximum newton steps.
MaxKrylov = 20;
errvecIN    = zeros(MaxIters,1);
errvecIN(1) = norm(Resid,2);
k = 1; err = Inf;
tola = 1e-6; tolr = 1e-4;
errold=norm(F(x),2);
errvecS(1) = errold;
m = 4;
skipMyGMRES = 0;
matlabGMRES = 0;
tol=tola + tolr*errvecIN(1);%Scale Relative tolerance by initial error.
rho = 1;
lambda_min = 1e-5;
alpha = 1e-4;
setup.udiag = 0;setup.thresh = 1;setup.type = 'nofill';setup.milu = 'off';setup.droptol = 0.1;
counter = 0;
% Newton Iteration
% Iterate our initial guess y until it satisfies F(y) = 0

while err > tol && k < MaxIters
    
    if mod(k,m) == 0 || rho > 0.5% Update the Jacobian every m iterations
        %                 fprintf(2,'Jacobian updated!\n');
        %Banded Jacobian
        J = JacobianFD_banded(F,x,Resid,lower,upper); % This function generates the banded jacobian.
        
        % Select preconditoner
        switch preCondSelector
            
            case 1
                M = J;
            case 2
                M = diag(diag(J));
            case 3
                M = diag(diag(J)) + tril(J,-1);
            case 4
                M = 0.25 * diag(diag(J)) + tril(J,-1);
            case 5
                M = 1.5 * diag(diag(J)) + tril(J,-1);
            case 6
                setup.type = 'ilutp';
                setup.milu = 'row';
                setup.droptol = 1e-15;
                [L,U] = ilu(J,setup);
            case 7
                skipMyGMRES = 1;
            case 8
                matlabGMRES = 1;
                skipMyGMRES = 1;
            otherwise
                M = J;
        end
    end
    
    if skipMyGMRES == 0
        % GMRES
        if preCondSelector == 6
            % nEeded for L U P
            [dx,~] = CW_pGMRES_SVD(J,-Resid,x,L,U,tol,MaxKrylov,false);
        else
            % USed for M
            [dx,m] = pGMRESSVD(J,-Resid,x,M,tol,MaxKrylov);
        end
        
    else
        %Backslash
        if matlabGMRES == 0
            dx = J\(-Resid);
        else % MATLAB GMRES
            [dx,~] = gmres(J,-Resid);
        end
    end
    
    % Linesearching
    lambda = 1;
    x_d = x + lambda*dx;
    Fxd = F(x_d);
    while (Fxd' * Fxd) >= (1-2*alpha*lambda)*(Resid' * Resid) && lambda > lambda_min
        sigma = 0.5;
        lambda = sigma*lambda;
        x_d = x + lambda*dx;
        Fxd = F(x_d);
    end
    x = x_d;
    % Error calculation
    Resid = F(x);
    errold = err;
    err = norm(Resid,2);
    
    rho=err/errold;
%         fprintf('%2.0f: ||F(y)|| = %1.8e | Threshold: %1.5e: Rho: %1.5e \n',counter,err,tol,rho);
    counter = counter + 1;
    k = k + 1;
    errvecIN(k) = err; % Store error at each Newton Iteration for plotting purposes.
    
    
    if err>tol && k == MaxIters
        exitflag(1) = 1;
        exitflag(2) = 0;
        x_out = x_return;
        return
        %error('Newton Method has failed as it exceeded the maximum number of iterations and tolerance was not met. \n MaxIters = %2.0f',k)
    end
end

x_out = x;

end
