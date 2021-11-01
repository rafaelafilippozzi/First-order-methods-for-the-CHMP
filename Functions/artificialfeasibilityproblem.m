function [averagetime, averageiteration, parameters, countfail] =  artificialfeasibilityproblem (m,n,viability, tolfix)
%% Artificial  
%         Generates artificial instances of linear programming feasibility
%         problem. The feasible set of a linear programming problem can be
%         given by Omega= {x \in \Rn: Ax=b, x>=0 and e^Tx <=N}
%         
% Syntax:
%         [averagetime, averageiteration, parameters, countfail] = artificialchmp(nvar,pvar)
%
% Input: 
%         m == the number of rows of matrix A
%         n == the number of columns of matrix A
%         viability == 1 problem is viable
%         viability == -1 problem is inviable
%         tolfix == tolerance epsilon
%
% Output:
%         averagetime: average time of r examples of all algorithms
%         averageiteration: average iterations of r examples of all algorithms
%         parameters:
%         countfail: Counter if there is failure in the algorithms (maxit,decision or linear search) 


rng(11);                         %Initialize seed for experiments to be reproducible
r = 10;                          %Number of times it will run with the pair (N,n)
M = 60;                          %M for the nonmonotone line search in SPG
maxit = 10^6;                    %Maximum number of iterations
N = 1200;                        %Boundary e^Tx<=N 
sigmin = 1e-10;                  %lowest acceptable spectral parameter
sigmax = 1e+10;                  %highest acceptable spectral parameter
Rmatrix = [];                    %Matrix with max of distance \tilde{A(:,i)} for p
delta0matrix = [];               %Matrix with min of distance \tilde{A(:,i)} for p
Deltaprojmatrix = [];            %Distance between proj of p in convA and p
counterawaystep = zeros(r,1);    % counter away step

%Vector for average times
    timeTA = zeros(r,1);                  %For Triangle Algorithm (TA)
    timeSPG = zeros(r,1);                 %For Spectral Projected Gradient (SPG)
    timeGT = zeros(r,1);                  %For Greedy Triangle (GT)
    timeASFW = zeros(r,1);                %For Away Step Frank Wolfe (ASFW) 
    timelinprog = zeros(r,1);             %For linprog

%Vector for average iterations  
    iterationsTA = zeros(r,1);             %For TA
    iterationsSPG = zeros(r,1);            %For SPG
    iterationsGT = zeros(r,1);             %For GT 
    iterationsASFW = zeros(r,1);           %For ASFW 
    iterationslinprog = zeros(r,1);        %For linprog

    %Fail every r rounds
        failTA = 0;                                  %For TA
        failSPG = 0;                                 %For SPG
        failGT = 0;                                  %For GT
        failASFW = 0;                                %For ASFW        
        faillinprog = 0;                             %For linprog

    % Fail decision
        faildecisionTA = 0;                          %For TA
        faildecisionSPG = 0;                         %For SPG
        faildecisionGT = 0;                          %For GT
        faildecisionASFW = 0;                        %For ASFW

for k = 1:r
    A = randn(m,n);            % returns an m-by-n matrix containing pseudorandom values drawn    from the standard normal distribution
    % normalizing A
    d = zeros(size(A,2),1);
    for i=1:size(A,2)
        d(i) = 1/norm(A(:,i),2);
        A(:,i)= A(:,i).*d(i);
    end
    A = A +1;                   %centered at ones(n,1)
    x = rand(n,1);              %random vector x with each entry from a uniform distribution in (0,1)
    b = A*x;    
    if viability == 1                
        viabilitylinprog = 1;   %for linprog fail
    else     
        b(1,1) = -b(1,1);       %multiply the first coordinate of b by?1
        viabilitylinprog = -2;  %for linprog fail 
    end
    
        %% The convex hull membership problem equivalence of Feasibility problem  
    Atilde = [A zeros(m,1) -b; ones(n,1)' 1 -N; zeros(1,n) 0  1];
    p = [zeros(m,1); 0; 1/(N+1)];
    
    %Initiation
    Raux = [];
    for l=1:(n+2)
        x=norm(Atilde(:,l)-p);
        Raux=[Raux x];
    end
    R=max(Raux);
    Rmatrix = [Rmatrix; R];
    [delta0,I]= min(Raux);
    delta0matrix = [delta0matrix ; delta0];
    i=I(1);                 %indice column p0
    p0 = Atilde(:,i);       %initial point for TA, GT and ASFW 
    x0 = zeros(n+2,1);
    x0(i)=1;                %initial point for SPG 
    L = norm(Atilde,2)^2;   %constante Lipschits for stopping criteria SPG 
    tol = tolfix*R;         %tolerance epsilon*R
    
       
    %Calculating projection in case p is out to estimate complexity
    if viability == -1
        [Deltaproj,~] = SPGwithoutCHMPstopcriteria(Atilde,p,x0,1e-7,maxit, M, sigmin,sigmax);
        Deltaprojmatrix = [Deltaprojmatrix Deltaproj]; 
    end
%% Algorithms 
    %Triangle Algorithm with Random Pivots
    tic, [~,~,Decision,iterationstav] = TriangleAlgorithm(Atilde,p,p0,tol,maxit);
    timeTAv = toc;
    if Decision == viability    
            timeTA(k,1) = timeTAv;
            iterationsTA(k,1) = iterationstav ;
    else
        failTA = failTA + 1;  %If TA reached the maxit or wrong decision 
        if Decision == -1*viability %If TA wrong decision
            faildecisionTA= faildecisionTA + 1;
        end
    end      

    % Greedy Triangle Algorithm 
    tic, [~,~,Decision,iterationsGTv] = GreedyTriangleAlgorithm(Atilde,p,p0,tol,maxit);
    timeGTv = toc;
    if Decision == viability    
            timeGT(k,1) = timeGTv;
            iterationsGT(k,1) = iterationsGTv ;
    else
        failGT = failGT + 1;  %If GT reached the maxit or wrong decision 
        if Decision == -1*viability %If GT wrong decision
            faildecisionGT= faildecisionGT + 1;
        end
    end      

    %Away Step Frank Wolfe Algorithm          
    tic,[~,~,Decision,iterationsASFWv, countas] = AwayStepFrankWolfeAlgorithm(Atilde,p,p0,i, tol,maxit);
    timeASFWv = toc;
    if Decision == viability    
       timeASFW(k,1) = timeASFWv;
       iterationsASFW(k,1) = iterationsASFWv ;
       counterawaystep (k,1) = countas;
    else
        failASFW = failASFW + 1;  %If ASFW reached the maxit or wrong decision 
        if Decision == -1*viability %If TA wrong decision
            faildecisionASFW= faildecisionASFW + 1;
        end
    end      

    %Spectral Projected Gradient  (SPG)
    tic,[Decision,iterationsSPGv] = SpectralProjectedGradient(Atilde,p,x0,tol,maxit, M, L,sigmin,sigmax);
    timeSPGv = toc;
    if Decision == viability    
       timeSPG(k,1) = timeSPGv;
       iterationsSPG(k,1) = iterationsSPGv ;
    else
        failSPG = failSPG + 1;  %If SPG reached the maxit or linear search failed or wrong decision 
        if Decision == -1*viability %If SPG wrong decision
            faildecisionSPG= faildecisionSPG + 1;
        end
    end       


    % Linprog
    options = optimoptions('linprog','TolCon',1e-3,'Preprocess','none');
    tic,  [~,~,flaglinprogaux, output] =  linprog([], ones(1,size(A,2)),N, A,b,zeros(size(A,2),1),[],options);
    timelinprogv = toc;

    if flaglinprogaux == viabilitylinprog    
       timelinprog(k,1) = timelinprogv;
       iterationsSPG(k,1) = output.iterations; 
    else
        faillinprog = faillinprog + 1;  %If SPG reached the maxit or linear search failed or wrong decision 
    end       
end            
    %TA
    averagetimeTA = sum(timeTA)/(r-failTA);
    averageiterationTA = sum(iterationsTA)/(r-failTA);
    countfaildecTA =  faildecisionTA;
    countfailTA = failTA;       

    %GT
    averagetimeGT = sum(timeGT)/(r-failGT);
    averageiterationGT = sum(iterationsGT)/(r-failGT);
    countfaildecGT =  faildecisionGT;
    countfailGT = failGT;

    %ASFW
    averagetimeASFW = sum(timeASFW)/(r-failASFW);
    averageiterationASFW = sum(iterationsASFW)/(r-failASFW);
    averagecounterawaystep = sum(counterawaystep)/(r-failASFW);
    countfaildecASFW =  faildecisionASFW;
    countfailASFW = failASFW;

    %SPG
    averagetimeSPG = sum(timeSPG)/(r-failSPG);
    averageiterationSPG = sum(iterationsSPG)/(r-failSPG);
    countfaildecSPG =  faildecisionSPG;
    countfailSPG = failSPG;

    %linprog
    averagetimelinprog = sum(timelinprog)/(r-faillinprog);
    averageiterationslinprog = sum(iterationslinprog)/(r-faillinprog);
    countfaillinprog = faillinprog;

    %% Structs
   %times:
   averagetime.TA = averagetimeTA;  
   averagetime.SPG = averagetimeSPG;
   averagetime.GT =  averagetimeGT; 
   averagetime.ASFW = averagetimeASFW;
   averagetime.linprog = averagetimelinprog;

   %iterations:
   averageiteration.TA = averageiterationTA;  
   averageiteration.SPG = averageiterationSPG;
   averageiteration.GT =  averageiterationGT; 
   averageiteration.ASFW = averageiterationASFW;
   averageiteration.linprog = averageiterationslinprog;

   %parameters
   parameters.deltaproj = Deltaprojmatrix;
   parameters.R = Rmatrix;
   parameters.delta0 = delta0matrix;
   parameters.counterAS = averagecounterawaystep;  

   %Countfail 
   countfail.allTA = countfailTA;
   countfail.allSPG = countfailSPG;
   countfail.allGT = countfailGT;
   countfail.allASFW = countfailASFW;
   countfail.alllinprog = countfaillinprog;
   countfail.decTA = countfaildecTA;
   countfail.decSPG = countfaildecSPG;
   countfail.decGT = countfaildecGT;
   countfail.decASFW = countfaildecASFW;

   