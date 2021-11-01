function [v,beta,Decision,iterationsta] = TriangleAlgorithm(A,p,p0,tol,maxit)
%% Triangle Algorithm with random pivots 
%
% Syntax:
%        [v,beta,Decision,iterationsta] = TriangleAlgorithmwithrandompivots(A,p,p0,tol,maxit) 
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
%         iterationsta: number of iterations
%
%% Initialization
    rng(123)
    [m, ~] = size(A);
    v = zeros(m,1);
    beta = 0;
    iterationsta = 0;
    Decision = 0;
    pk = p0;
    pkp = pk - p;
    optm = norm(pkp,2);
    normp2 = norm(p,2)^2;

%% loop principal    
while (optm > tol)                               %while pk it's not a p_epsilon-solution
    normpk2 = norm(pk,2)^2; 
    prod = A'*pkp;
    idxj = find(2*prod - normpk2 + normp2 < 0);  %find simple pivots
    if ~isempty(idxj)                            %if exist pivot simple
        k = length(idxj);
        j = idxj(randi(k,1));                    %randomize the choice of pivot
        vj = A(:,j);
        diff = vj - pk;
        dist = norm(diff, 2);
        alpha = min(1,(pkp'*(-diff))/(dist^2));  %Linear search
        pk = (1-alpha)*pk + alpha*vj;            %new iteration
        pkp = pk - p;
        optm = norm(pkp,2);
        iterationsta=iterationsta+1;
        if (iterationsta > maxit)
            fprintf('TA: maximum number of iterations reached!\n')
            return
        end
    else
        v = p-pk;                               
        normv = norm(v)^2/2;
        beta = v'*pk + normv;                    
        Decision = -1;
        fprintf('TA: The point p is OUTSIDE of conv(A)\n')
        return
    end
end
        Decision = 1;
        fprintf('TA: The point p is INSIDE of conv(A)\n')
        return;
