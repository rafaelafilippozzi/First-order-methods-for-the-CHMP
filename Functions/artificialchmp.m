function [averagetime, averageiteration, parameters, countfail] = artificialchmp(nvar,pvar)
%% Artificial CHMP  
%         Generates artificial instances of CHMP that consists in deciding 
%         whether p lies of conv(A) where each element of A={v1,...,vn} is 
%         randomly generated according to a uniform distribution in the unit
%         ball of Rm and the query point p in Rm is generated as chosen in Input
%
% Syntax:
%         [averagetime, averageiteration, parameters, countfail] = artificialchmp(nvar,pvar)
%
% Input: 
%         nvar: Amount of points. 
%         nvar == 1 (500,1000,1500...,10000), nvar == 2 (10000,20000,30000...,100000), 
%         pvar: Cases of generating p. 
%         pvar == 1 p in the relative interior of conv(A);  
%         pvar == 2 p in conv(A) with visibility factor close to zero;
%         pvar == 3 p lies outside of conv(A) and away from the boundary; 
%         pvar == 4 p lies outside of conv(A) and with visibility factor close to zero.
%
% Output:
%         averagetime: average time of r examples of all algorithms
%         averageiteration: average iterations of r examples of all algorithms
%         parameters: Deltaprojmatrix, Rmatrix, delta0matrix, averagecounterawaystep
%         countfail: Counter if there is failure in the algorithms (maxit,decision or linear search) 

%rng(123);                      %Initialize seed for experiments to be reproducible
tolfix = 1e-4;                 %Chosen epsilon tolerance 
m = 100;                       %Dimension 
r = 10;                        %Number of times it will run with the same pair (m,n)
M = 15;                        %M for the nonmonotone line search in SPG
sigmin = 1e-8;                 %lowest acceptable spectral parameter
sigmax = 1e+8;                 %highest acceptable spectral parameter
Rmatrix = [];                  %Matrix with max of distance A(:,i) for p
delta0matrix = [];             %Matrix with min of distance A(:,i) for p
Deltaprojmatrix = [];          %Distance between proj of p in convA and p 
nend = 10;                     %number of multiples for chosen nvar
averagecounterawaystep = zeros(nend,1);          %Average counter away step

%Matrix for minimum, maximum, total time averages 
    averagetimeTA = zeros(nend,3);                  %For Triangle Algorithm (TA)
    averagetimeSPG = zeros(nend,3);                 %For Spectral Projected Gradient (SPG)
    averagetimeGT = zeros(nend,3);                  %For Greedy Triangle (GT)
    averagetimeASFW = zeros(nend,3);                %For Away Step Frank Wolfe (ASFW) 

%Vector for average iterations  
    averageiterationTA = zeros(nend,1);             %For TA
    averageiterationSPG = zeros(nend,1);            %For SPG
    averageiterationGT = zeros(nend,1);             %For GT 
    averageiterationASFW = zeros(nend,1);           %For ASFW 

%Counter if there is failure of decision
    countfaildecTA = zeros(nend,1);                 %For TA      
    countfaildecSPG = zeros(nend,1);                %For SPG      
    countfaildecGT = zeros(nend,1);                 %For GT      
    countfaildecASFW = zeros(nend,1);               %For ASFW    

%Quantidade de falhas
    countfailTA = zeros(nend,1);                    %For TA 
    countfailSPG = zeros(nend,1);                   %For SPG 
    countfailGT = zeros(nend,1);                    %For GT 
    countfailASFW = zeros(nend,1);                  %For ASFW 

for j=1:nend
    if nvar ==1 
        n = 500*j;
    else
        n = 10000*j;
    end

    maxit = min(max(1000*n, 10000),1e+6);            %Maximum number of iterations 
    counterawatstepforrrounds = zeros(r,1);          %For counter away step        
    fprintf('\n ======= (n = %d) ====== \n',n);


    %Time every r rounds
        timeTA = zeros(r,1);                         %For TA
        timeSPG = zeros(r,1);                        %For SPG
        timeGT = zeros(r,1);                         %For GT
        timeASFW = zeros(r,1);                       %For ASFW

    %Iterations every r rounds
        iterationsTA = zeros(r,1);                   %For TA
        iterationsSPG = zeros(r,1);                  %For SPG
        iterationsGT = zeros(r,1);                   %For GT
        iterationsASFW = zeros(r,1);                 %For ASFW

    %Fail every r rounds
        failTA = 0;                                  %For TA
        failSPG = 0;                                 %For SPG
        failGT = 0;                                  %For GT
        failASFW = 0;                                %For ASFW        

    % Fail decision
        faildecisionTA = 0;                          %For TA
        faildecisionSPG = 0;                         %For SPG
        faildecisionGT = 0;                          %For GT
        faildecisionASFW = 0;                        %For ASFW
    for k = 1:r
        rng(10*k)                           %Initialize seed for experiments to be reproducible
        [ A ] = generateArandom( m,n );     %each element of A={v1,...,vn} randomly generated according to a uniform distribution in the unit ball of Rm
        %Generate a query point (m,1)   
        if pvar == 1
           p = max(2*ones(m,1)+ randn(m,1),0);
           p = p/(n*norm(p));              %p with a lot of chance of being inside the relative interior of the convA
           answer  = 1;
        else
           %choose two points of A for which e'vi are the highest
           aux = sum(A,1);    
           [~, idx] = sort(aux,'descend');
           qidx = idx(1);
           lidx = idx(2);               
            if pvar == 3
               p = 1.5*(A(:,qidx)+ A(:,lidx))/2;
               answer  = -1;  
            else
               p = 0.5*A(:,qidx) + 0.5*A(:,lidx); 
               %for qidx or lidx not be the starting point:
               dist = norm(A(:,qidx) - A(:,lidx),2);
               d0 = 0.9*(dist/2);
               normp = norm(p,2);
               alpha =  normp - d0;
               vs = alpha*(p/normp); %p is closer to vs than to qidx or lidx
               A = [A vs]; 
               n=n+1;
               if pvar == 2 
                  answer  = 1;
               else if pvar == 4
                       p = 1.01*(A(:,qidx)+ A(:,lidx))/2;
                       answer  = -1;
                    end
               end
            end
        end

        %Initiation
        Raux = [];
        for l=1:n
            x=norm(A(:,l)-p);
            Raux=[Raux x];
        end
        R=max(Raux);
        Rmatrix = [Rmatrix; R];
        [delta0,I]= min(Raux);
        delta0matrix = [delta0matrix ; delta0];
        i=I(1);           %indice column p0
        p0 = A(:,i);      %initial point for TA, GT and ASFW 
        x0 = zeros(n,1);
        x0(i)=1;          %initial point for SPG 
        L = norm(A,2)^2;  %constante Lipschits for stopping criteria SPG 
        tol = tolfix*R;   %tolerance epsilon*R

        %Calculating projection in case p is out to estimate complexity
        if answer == -1
        [Deltaproj,~] = SPGwithoutCHMPstopcriteria(A,p,x0,tol,maxit, M, sigmin,sigmax);
        Deltaprojmatrix = [Deltaprojmatrix Deltaproj]; 
        end

%% Algorithms 
        %Triangle Algorithm with Random Pivots
        tic, [~,~,Decision,iterationstav] = TriangleAlgorithm(A,p,p0,tol,maxit);
        timeTAv = toc;
        if Decision == answer    
                timeTA(k,1) = timeTAv;
                iterationsTA(k,1) = iterationstav ;
        else
            failTA = failTA + 1;  %If TA reached the maxit or wrong decision 
            if Decision == -1*answer %If TA wrong decision
                faildecisionTA= faildecisionTA + 1;
            end
        end      

        % Greedy Triangle Algorithm 
        tic, [~,~,Decision,iterationsGTv] = GreedyTriangleAlgorithm(A,p,p0,tol,maxit);
        timeGTv = toc;
        if Decision == answer    
                timeGT(k,1) = timeGTv;
                iterationsGT(k,1) = iterationsGTv ;
        else
            failGT = failGT + 1;  %If GT reached the maxit or wrong decision 
            if Decision == -1*answer %If GT wrong decision
                faildecisionGT= faildecisionGT + 1;
            end
        end      

        %Away Step Frank Wolfe Algorithm          
        tic,[~,~,Decision,iterationsASFWv, countas] = AwayStepFrankWolfeAlgorithm(A,p,p0,i, tol,maxit);
        timeASFWv = toc;
        if Decision == answer    
           timeASFW(k,1) = timeASFWv;
           iterationsASFW(k,1) = iterationsASFWv ;
           counterawatstepforrrounds (k,1) = countas;
        else
            failASFW = failASFW + 1;  %If ASFW reached the maxit or wrong decision 
            if Decision == -1*answer %If TA wrong decision
                faildecisionASFW= faildecisionASFW + 1;
            end
        end      

        %Spectral Projected Gradient  (SPG)
        tic,[Decision,iterationsSPGv] = SpectralProjectedGradient(A,p,x0,tol,maxit, M, L,sigmin,sigmax);
        timeSPGv = toc;
        if Decision == answer    
           timeSPG(k,1) = timeSPGv;
           iterationsSPG(k,1) = iterationsSPGv ;
        else
            failSPG = failSPG + 1;  %If SPG reached the maxit or linear search failed or wrong decision 
            if Decision == -1*answer %If SPG wrong decision
                faildecisionSPG= faildecisionSPG + 1;
            end
        end                 
    end        
             fprintf('\n ..................... \n \n');               
       
        %TA
        averagetimeTA(j,1) = min(timeTA);
        averagetimeTA(j,2) = max(timeTA);
        averagetimeTA(j,3) = sum(timeTA)/(r-faildecisionTA);
        averageiterationTA(j,1) = sum(iterationsTA)/(r-faildecisionTA);
        countfaildecTA(j,1) =  faildecisionTA;
        countfailTA(j,1) = failTA;       
        
        %GT
        averagetimeGT(j,1) = min(timeGT);
        averagetimeGT(j,2) = max(timeGT);
        averagetimeGT(j,3) = sum(timeGT)/(r-faildecisionGT);
        averageiterationGT(j,1) = sum(iterationsGT)/(r-faildecisionGT);
        countfaildecGT(j,1) =  faildecisionGT;
        countfailGT(j,1) = failGT;
        
        %ASFW
        averagetimeASFW(j,1) = min(timeASFW);
        averagetimeASFW(j,2) = max(timeASFW);
        averagetimeASFW(j,3) = sum(timeASFW)/(r-faildecisionASFW);
        averageiterationASFW(j,1) = sum(iterationsASFW)/(r-faildecisionASFW);
        averagecounterawaystep(j,1) = sum(counterawatstepforrrounds)/(r-faildecisionASFW);
        countfaildecASFW(j,1) =  faildecisionASFW;
        countfailASFW(j,1) = failASFW;
        
        %SPG
        averagetimeSPG(j,1) = min(timeSPG);
        averagetimeSPG(j,2) = max(timeSPG);
        averagetimeSPG(j,3) = sum(timeSPG)/(r-faildecisionSPG);
        averageiterationSPG(j,1) = sum(iterationsSPG)/(r-faildecisionSPG);
        countfaildecSPG(j,1) =  faildecisionSPG;
        countfailSPG(j,1) = failSPG;
        
        %% Structs
       %times:
       averagetime.TA = averagetimeTA;  
       averagetime.SPG = averagetimeSPG;
       averagetime.GT =  averagetimeGT; 
       averagetime.ASFW = averagetimeASFW;

       %iterations:
       averageiteration.TA = averageiterationTA;  
       averageiteration.SPG = averageiterationSPG;
       averageiteration.GT =  averageiterationGT; 
       averageiteration.ASFW = averageiterationASFW;

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
       countfail.decTA = countfaildecTA;
       countfail.decSPG = countfaildecSPG;
       countfail.decGT = countfaildecGT;
       countfail.decASFW = countfaildecASFW;

        fprintf(' ===================== \n');
end
