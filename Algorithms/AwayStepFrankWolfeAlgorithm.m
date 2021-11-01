function [v,beta,Decision,iterationsASFW, countas] = AwayStepFrankWolfeAlgorithm(A,p,p0,i, tol,maxit)
%% Away Step Frank Wolfe Algorithm 
%
% Syntax: 
%       [v,beta,Decision,iterationsASFW, countas] = AwayStepFrankWolfeAlgorithm(A,p,p0,i, tol,maxit)
%
% Input: 
%         A: set of n points in R^m
%         p: a query point in R^m 
%         p0: initial point
%         i: column index of p0
%         tol: stop criterion to find a p_epsilon-solution
%         maxit: Maximum number of iterations
%
% Output: 
%         v: vector normal to the separating hyperplane (v=0 if p belongs to conv(A))
%         beta: vector normal to separator hyperplane (beta= v'x)
%         Decision: 1 when p \in conv(A)  
%         Decision: -1 when p \notin conv(A)
%         Decision: 0 when the Algorithm reaches maxit 
%         iterationsASFW: number of iterations
%         countas: step away counter
%
%% Initialization
    [m, n] = size(A);
    a = sparse(n,1); % vector with convex combination coefficients
    a(i) = 1;
    v = zeros(m,1);
    beta = 0;
    iterationsASFW = 0;
    Decision = 0;
    countas = 0;
    pk = p0;
    pkp = pk - p;
    optm = norm(pkp,2);
    normp2 = norm(p,2)^2;

%% loop principal    
while (optm > tol)               %while pk it's not a p_epsilon-solution
    U = find(a>0);
    prod = A'*pkp;
 
    % Frank-Wolfe point 
    [~, idxj] = min(prod);
    j = idxj(1);
    vj = A(:,j);  
    % Away point
    [~, idxk] = max(prod(U));
    k = U(idxk(1));
    wk = A(:,k);
    
    normpk2 = norm(pk,2)^2;     
    diff = pk - vj;
    dist = norm(diff, 2);
    
    if 2*vj'*pkp <= (normpk2 - normp2)          % Simple pivot characterization 
        dfw = -diff;                            % FW direction
        gtdfw = pkp'*dfw;                       % gradient in FW direction
        da = pk - wk;                           % Away direction
        gtda = pkp'*da;                         % gradient in Away direction
        away=0;
       
        if (gtdfw <= gtda)                      %Choose FW Step 
            amax = 1;
            alpha = min(amax, -gtdfw/(dist^2)); %Line search
            d = dfw;                            %Direction
        else                                    %Choose Away Step 
            away=1;
            countas = countas +1;
            amax = a(k)/(1 - a(k));             %to ensure the point remains on the convA
            alpha = min(amax, -gtda/(dist^2));  %Line search
            d = da;                             %Direction
        end
        pk = pk + alpha*d;                      %New iterate       
        % update convex combination coefficients
        if ~away
            % FW step
            if (alpha == amax)
                a = 0*a;
                a(j) = 1;
            else
                a = (1-alpha)*a;
                a(j) = a(j) + alpha;
            end
        else
            % away step
            a = (1 + alpha)*a;
            if (alpha == amax)
                a(k) = 0;
            else
                a(k) = a(k) - alpha;
            end
        end        
        pkp = pk - p;
        optm = norm(pkp,2);
        iterationsASFW=iterationsASFW+1;
            if (iterationsASFW > maxit)
            fprintf('ASFW: maximum number of iterations reached!\n')
            return
            end
    else 
        v = p-pk;                               
        normv = norm(v)^2/2;
        beta = v'*pk + normv;                    
        Decision = -1;
        fprintf('ASFW: The point p is OUTSIDE of conv(A)\n')
        return
    end
end
        Decision = 1;
        fprintf('ASFW: The point p is INSIDE of conv(A)\n')
        return;
