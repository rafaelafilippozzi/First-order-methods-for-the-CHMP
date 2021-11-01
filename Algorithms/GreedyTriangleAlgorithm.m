function [v,beta,Decision,iterationsGT] = GreedyTriangleAlgorithm(A,p,p0,tol,maxit)
%% Greedy Triangle Algorithm 
%
% Syntax: 
%       [v,beta,Decision,iterationsGT] = GreedyTriangleAlgorithm(A,p,p0,tol,maxit)
%
% Input: 
%         A: set of n points in R^m
%         p: a query point in R^m 
%         p0: initial point
%         tol: stop criterion to find a p_epsilon-solution
%         maxit: Maximum number of iterations
%
% Output: 
%         v: vector normal to the separating hyperplane (v=0 if p belongs to conv(A))
%         beta: vector normal to separator hyperplane (beta= v'x)
%         Decision: 1 when p \in conv(A)  
%         Decision: -1 when p \notin conv(A)
%         Decision: 0 when the Algorithm reaches maxit 
%         iterationsGT: number of iterations
%
%% Initialization
    [m, ~] = size(A);
    v = zeros(m,1);
    beta = 0;
    iterationsGT = 0;
    Decision = 0;
    pk = p0;
    pkp = pk - p;
    optm = norm(pkp,2);
    normp2 = norm(p,2)^2;

%% loop principal    
while (optm > tol)               %while pk it's not a p_epsilon-solution
    prod = A'*pkp;
    [~, idxj] = min(prod);       %Choose vj \in \argmin {vi^T(pk - p): v_i \in A/{pk}}
    j = idxj(1);
    vj = A(:,j);
    normpk2 = norm(pk,2)^2; 
    diff = vj - pk;
    dist = norm(diff, 2);  
    if 2*vj'*pkp <= (normpk2 - normp2)           % Simple pivot characterization 
        alpha = min(1,(pkp'*(-diff))/(dist^2));  %Linear search
        pk = (1-alpha)*pk + alpha*vj;            %new iteration
        pkp = pk - p;
        optm = norm(pkp,2);
        iterationsGT=iterationsGT+1;
        if (iterationsGT > maxit)
            fprintf('GT: maximum number of iterations reached!\n')
            return
        end
    else 
        v = p-pk;                               
        normv = norm(v)^2/2;
        beta = v'*pk + normv;                    
        Decision = -1;
        fprintf('GT: The point p is OUTSIDE of conv(A)\n')
        return
    end
end
        Decision = 1;
        fprintf('GT: The point p is INSIDE of conv(A)\n')
        return;
