function [Decision,iterationsSPG]= SpectralProjectedGradient(A,p,x0,tol,maxit, M, L,sigmin,sigmax)  
%% Spectral Projected Gradient 
%
% Syntax: 
%       [Decision,iterationsSPG]= SpectralProjectedGradient(A,p,x0,tol,maxit, M, L,sigmin,sigmax)
%
% Input: 
%         A: set of n points in R^m
%         p: a query point in R^m 
%         x0: initial point (Ax0 = p0)
%         tol: stop criterion to find a p_epsilon-solution
%         maxit: Maximum number of iterations
%         M:  for the nonmonotone line search in SPG
%         L: norm(A)^2 constante Lipschits for stopping criteria SPG 
%         sigmin: lowest acceptable spectral parameter
%         sigmax: highest acceptable spectral parameter
%
% Output: 
%         Decision: 1 when p \in conv(A)  
%         Decision: -1 when p \notin conv(A)
%         Decision: 0 when the Algorithm reaches maxit
%         Decision: -2 when the Algorithm fail on linear search  
%         iterationsSPG: number of iterations
%
%% Initialization 
    Ap = A'*p;
    Decision = 0;
    x = x0;
    %IterationsSPG = 1 (First iteration SPG)
    z = A*x;
    diff = z - p;
    normdiff = norm(diff,2);
    grad = A'*diff;                          %Gradient current iteration
    fn = 0.5*normdiff^2;
    fnsave = fn*ones(M,1);                   %Save objective M function values
    xold = x;                                %x_k
    gradold = grad;                           %Grad x_k
    iterationsSPG = 1;
    normp = norm(p,2);

%% loop principal
while  (norm (z-p)> tol)  %while pk it's not a p_epsilon-solution
    Ok = -1;
    if iterationsSPG>1
        sdif = x-xold;                 %x_{k+1}-x_k
        ydif = grad-gradold;           % Gradient(x_{k+1})- Gradient(x_k)
        sig = min (sigmax, max(sigmin, (sdif'*ydif)/(sdif'*sdif)));   % spectral parameter
    else
        sig = 1;                       %Firt Iteration (Classical Gradient Project) 
    end
    proj = simplex_proj(x-(1/sig)*grad);    %Unit simplex projection
    zz = A*proj;
    %Find descent direction
    d = proj - x;                      %Direction

    %Stop Criterion equivalent to TA        
    if norm(d,2)<= norm(zz-p,2)*(tol)/(3*L*sqrt(2))
        Decision = 1;
        fprintf( 'SPG:  The point p is INSIDE of conv(A)\n' );
        return
    end

    %Criteria for stopping out (checking if z is a witness)       
    normz=norm(zz,2);
    if any(Ap - A'*zz >= (normp^2 - normz^2)/2)
        Ok=1;
    end     
    if Ok~= 1
      Decision=-1;
      fprintf('SPG: The point p is OUTSIDE of conv(A)\n');
      return
    end

    % Nonmonotone line-search
    j= mod (iterationsSPG,M);
    fnsave(j+1,1) = fn;    
    fmax= max(fnsave);  

    xt = x + d;
    z = A*xt;  
    diff = z-p;                                
    normdiff = norm(diff,2);                     %Norm of objective function                                                 %Gradiente nova iteração
    fn = 0.5*normdiff^2;                         %Objective function value
    gtd = grad'*d;                               

    %Linear search to find tk and xnew                                          
    t=1;
    [ xt, diff, fn, flaglin] = nonmonotoneArmijocriterion( fn, fmax, t, 1e-4, gtd, x,d,p, xt,diff, A);
    if flaglin < 0
        Decision = -2;
        fprintf('SPG: linear search failed!\n')
        return
    end

    xold= x;
    gradold= grad;
    x=xt;
    grad = A'*diff;                          %Gradiente iteração atual
    iterationsSPG = iterationsSPG+1;                                     %Contador iteração

    if (iterationsSPG > maxit)
        fprintf('SPG: número máximo de iterações atingido!\n')
        return
    end  
end
    Decision = 1;
    fprintf('SPG: The point p is INSIDE of conv(A)\n')
    return;