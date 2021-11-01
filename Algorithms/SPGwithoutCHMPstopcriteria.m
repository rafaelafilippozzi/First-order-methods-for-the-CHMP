function [normdiff, Decision]= SPGwithoutCHMPstopcriteria(A,p,x0,tol,maxit, M, sigmin,sigmax)  
%% Spectral Projected Gradient without CHMP stop criteria
%
% Syntax: 
%       [normdiff]= SPGprojecao_semcriteriotestemunha(A,p,x0,tol,maxit, M, sigmin,sigmax)
%
% Input: 
%         A: set of n points in R^m
%         p: a query point in R^m 
%         x0: initial point (Ax0 = p0)
%         tol: stop criterion to find a p_epsilon-solution
%         maxit: Maximum number of iterations
%         M:  for the nonmonotone line search in SPG
%         sigmin: lowest acceptable spectral parameter
%         sigmax: highest acceptable spectral parameter
%
% Output:
%         normdiff: Distance between p and the projection of p on the conv(A)
%         Decision: 1 when p \in conv(A)  
%         Decision: -1 when p \notin conv(A)
%         Decision: 0 when the Algorithm reaches maxit
%         Decision: -2 when the Algorithm fail on linear search  
%
%% Initialization 
    x = x0;
    Decision = 0;
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

%% loop principal
while  (norm (z-p)> tol)  %while pk it's not a p_epsilon-solution
    if iterationsSPG>1
        sdif = x-xold;                 %x_{k+1}-x_k
        ydif = grad-gradold;           % Gradient(x_{k+1})- Gradient(x_k)
        sig = min (sigmax, max(sigmin, (sdif'*ydif)/(sdif'*sdif)));   % spectral parameter
    else
        sig = 1;                       %Firt Iteration (Classical Gradient Project) 
    end
    proj = simplex_proj(x-(1/sig)*grad);    %Unit simplex projection
    %Find descent direction
    d = proj - x;                      %Direction      
    if(norm(d,2) < tol)
      fprintf('SPG: The point p is OUTSIDE of conv(A)\n');
      Decision = -1;
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
        fprintf('SPG: linear search failed!\n')
        Decision = -2;
        return
    end

    xold= x;
    gradold= grad;
    x=xt;
    grad = A'*diff;                          %Gradiente iteração atual
    iterationsSPG = iterationsSPG+1;                                     %Contador iteração

    if (iterationsSPG > maxit)
        fprintf('SPG: número máximo de iterações atingido!\n')
        Decision = 0;
        return
    end  
end
    fprintf('SPG: The point p is INSIDE of conv(A)\n')
    return;

