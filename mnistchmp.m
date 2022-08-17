function [averagetime, averageiteration, norms, parameters, countfail] = mnistchmp(datatrain,datatest)
%% MNIST CHMP  
%         
% Syntax:
%         [averagetime, averageiteration, parameters, countfail] = artificialchmp(nvar,pvar)
%
% Input: 
%         datatrain: Matrix with datatrain points. 
%         pdatatest: Point datatest to test
%
% Output:
%         averagetime: average time of r examples of all algorithms
%         averageiteration: average iterations of r examples of all algorithms
%         parameters: Deltaprojmatrix, Rmatrix, delta0matrix, averagecounterawaystep
%         countfail: Counter if there is failure in the algorithms (maxit,decision or linear search) 

%rng(123);                      %Initialize seed for experiments to be reproducible
tolfix = 1e-8;                 %Chosen epsilon tolerance 
m = 784;                       %Dimension 
n = 60000;
r = 1;                        %Number of times it will run with the same pair (m,n)
M = 3;                        %M for the nonmonotone line search in SPG
sigmin = 1e-8;                 %lowest acceptable spectral parameter
sigmax = 1e+8;                 %highest acceptable spectral parameter
Rmatrix = [];                  %Matrix with max of distance A(:,i) for p
delta0matrix = [];             %Matrix with min of distance A(:,i) for p
Deltaprojmatrix = [];          %Distance between proj of p in convA and p 
nend = 10000;                     %number of multiples for chosen nvar
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

%Vector for average iterations  
    averageiterationTA = zeros(nend,1);             %For TA
    averageiterationSPG = zeros(nend,1);            %For SPG
    averageiterationGT = zeros(nend,1);             %For GT 
    averageiterationASFW = zeros(nend,1);           %For ASFW 
    
%Vector for average iterations  
    averagenpkpTA = zeros(nend,1);             %For TA
    averagenpkpSPG = zeros(nend,1);            %For SPG
    averagenpkpGT = zeros(nend,1);             %For GT 
    averagenpkpASFW = zeros(nend,1);           %For ASFW 
    
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

[ A ] = datatrain;
for j=1:nend  
    p = datatest(:,j);
    maxit = min(max(1000*n, 10000),1e+6);            %Maximum number of iterations 
    counterawatstepforrrounds = zeros(r,1);          %For counter away step        
    fprintf('\n ======= (j = %d) ====== \n',j);


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
    
    %Normas pk-p    
        npkpTA = zeros(r,1);                   %For TA
        npkpSPG = zeros(r,1);                  %For SPG
        npkpGT = zeros(r,1);                   %For GT
        npkpASFW = zeros(r,1);                 %For ASFW
    
    
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
        
        answer = -1;

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
        
   for k = 1:r
        rng(10*k)                           %Initialize seed for experiments to be reproducible
        %% Algorithms 
        %Triangle Algorithm with Random Pivots
        tic, [vta,~,Decision,iterationstav] = TriangleAlgorithm(A,p,p0,tol,maxit);
        timeTAv = toc;
        if Decision == answer    
                timeTA(k,1) = timeTAv;
                iterationsTA(k,1) = iterationstav ;
                npkpTA(k,1) = norm(vta,'fro');
        else
            failTA = failTA + 1;  %If TA reached the maxit or wrong decision 
            if Decision == -1*answer %If TA wrong decision
                faildecisionTA= faildecisionTA + 1;
            end
        end      

        % Greedy Triangle Algorithm 
        tic, [vgt,~,Decision,iterationsGTv] = GreedyTriangleAlgorithm(A,p,p0,tol,maxit);
        timeGTv = toc;
        if Decision == answer    
                timeGT(k,1) = timeGTv;
                iterationsGT(k,1) = iterationsGTv;
                npkpGT(k,1) = norm(vgt,'fro');
        else
            failGT = failGT + 1;  %If GT reached the maxit or wrong decision 
            if Decision == -1*answer %If GT wrong decision
                faildecisionGT= faildecisionGT + 1;
            end
        end      

        %Away Step Frank Wolfe Algorithm          
        tic,[vasfw,~,Decision,iterationsASFWv, countas] = AwayStepFrankWolfeAlgorithm(A,p,p0,i, tol,maxit);
        timeASFWv = toc;
        if Decision == answer    
           timeASFW(k,1) = timeASFWv;
           iterationsASFW(k,1) = iterationsASFWv ;
           counterawatstepforrrounds (k,1) = countas;
           npkpASFW(k,1) = norm(vasfw,'fro');
        else
            failASFW = failASFW + 1;  %If ASFW reached the maxit or wrong decision 
            if Decision == -1*answer %If TA wrong decision
                faildecisionASFW= faildecisionASFW + 1;
            end
        end      

        %Spectral Projected Gradient  (SPG)
        tic,[Decision,iterationsSPGv,vspg] = SpectralProjectedGradient(A,p,x0,tol,maxit, M, L,sigmin,sigmax);
        timeSPGv = toc;
        if Decision == answer    
           timeSPG(k,1) = timeSPGv;
           iterationsSPG(k,1) = iterationsSPGv ;
           npkpSPG(k,1) = norm(vspg,'fro');
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
        averagenpkpTA(j,1) = sum(npkpTA)/(r-faildecisionTA);
        
        %GT
        averagetimeGT(j,1) = min(timeGT);
        averagetimeGT(j,2) = max(timeGT);
        averagetimeGT(j,3) = sum(timeGT)/(r-faildecisionGT);
        averageiterationGT(j,1) = sum(iterationsGT)/(r-faildecisionGT);
        countfaildecGT(j,1) =  faildecisionGT;
        countfailGT(j,1) = failGT;
        averagenpkpGT(j,1) = sum(npkpGT)/(r-faildecisionGT);
        
        %ASFW
        averagetimeASFW(j,1) = min(timeASFW);
        averagetimeASFW(j,2) = max(timeASFW);
        averagetimeASFW(j,3) = sum(timeASFW)/(r-faildecisionASFW);
        averageiterationASFW(j,1) = sum(iterationsASFW)/(r-faildecisionASFW);
        averagecounterawaystep(j,1) = sum(counterawatstepforrrounds)/(r-faildecisionASFW);
        countfaildecASFW(j,1) =  faildecisionASFW;
        countfailASFW(j,1) = failASFW;
        averagenpkpASFW(j,1) = sum(npkpASFW)/(r-faildecisionASFW);
        
        %SPG
        averagetimeSPG(j,1) = min(timeSPG);
        averagetimeSPG(j,2) = max(timeSPG);
        averagetimeSPG(j,3) = sum(timeSPG)/(r-faildecisionSPG);
        averageiterationSPG(j,1) = sum(iterationsSPG)/(r-faildecisionSPG);
        countfaildecSPG(j,1) =  faildecisionSPG;
        countfailSPG(j,1) = failSPG;
        averagenpkpSPG(j,1) = sum(npkpSPG)/(r-faildecisionSPG);
        
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
       
       %norms: 
       norms.TA = averagenpkpTA;
       norms.SPG = averagenpkpSPG;
       norms.GT = averagenpkpGT;
       norms.ASFW = averagenpkpASFW;
       
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
