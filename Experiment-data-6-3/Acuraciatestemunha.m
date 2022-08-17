function [averagetime, averageiteration, norms, parameters, countfail, labels,acuracia] = Acuraciatestemunha(datatrain,datatest)
%Acuracia com testemunha
tolfix = 1e-8;                 %Chosen epsilon tolerance 
m = 784;                       %Dimension 
r = 10; 
M = 3;                         %M for the nonmonotone line search in SPG
sigmin = 1e-8;                 %lowest acceptable spectral parameter
sigmax = 1e+8;                 %highest acceptable spectral parameter
nend = 10000;                  %number of multiples for chosen nvar
            

    Rmatrix = zeros(nend,r);                  %Matrix with max of distance A(:,i) for p
    delta0matrix = zeros(nend,r);             %Matrix with min of distance A(:,i) for p
    Deltaprojmatrix = zeros(nend,r);          %Distance between proj of p in convA and p 
    averagecounterawaystep = zeros(nend,r);          %Average counter away step

        %Vector for time
            averagetimeTA = zeros(nend,r);                  %For Triangle Algorithm (TA)
            averagetimeSPG = zeros(nend,r);                 %For Spectral Projected Gradient (SPG)
            averagetimeGT = zeros(nend,r);                  %For Greedy Triangle (GT)
            averagetimeASFW = zeros(nend,r);                %For Away Step Frank Wolfe (ASFW) 
            averagetimeProj = zeros(nend,r);                %For SPG sem criterio de parada 

            
        %Vector for iterations  
            averageiterationTA = zeros(nend,r);             %For TA
            averageiterationSPG = zeros(nend,r);            %For SPG
            averageiterationGT = zeros(nend,r);             %For GT 
            averageiterationASFW = zeros(nend,r);           %For ASFW 

        %Vector norm pkp  
            averagenpkpTA = zeros(nend,r);             %For TA
            averagenpkpSPG = zeros(nend,r);            %For SPG
            averagenpkpGT = zeros(nend,r);             %For GT 
            averagenpkpASFW = zeros(nend,r);           %For ASFW 

        %Counter if there is failure of decision
            countfaildecTA = zeros(nend,r);                 %For TA      
            countfaildecSPG = zeros(nend,r);                %For SPG      
            countfaildecGT = zeros(nend,r);                 %For GT      
            countfaildecASFW = zeros(nend,r);               %For ASFW    

        %Quantidade de falhas
            countfailTA = zeros(nend,r);                    %For TA 
            countfailSPG = zeros(nend,r);                   %For SPG 
            countfailGT = zeros(nend,r);                    %For GT 
            countfailASFW = zeros(nend,r);                  %For ASFW 
    
for k = 0:9
    I = find(datatrain(1,:)==k);
    n = size(I,2);
    maxit = min(max(1000*n, 10000),1e+6);            %Maximum number of iterations
    [ A ] = datatrain(2:end,I);
    L = norm(A,2)^2;               %constante Lipschits for stopping criteria SPG   
    r = k+1;
    for j=1:nend  
        p = datatest(2:end,j);
         
        counterawatstepforrrounds = 0;          %For counter away step        
        fprintf('\n ======= (j = %d) ====== \n',j);

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
            Rmatrix(j,r) = R;
            [delta0,I]= min(Raux);
            delta0matrix(j,r) = delta0;
            i=I(1);           %indice column p0
            p0 = A(:,i);      %initial point for TA, GT and ASFW 
            x0 = zeros(n,1);
            x0(i)=1;          %initial point for SPG 
            tol = tolfix*R;   %tolerance epsilon*R

            %Calculating projection in case p is out to estimate complexity
            if answer == -1
            tic,[Deltaproj,~] = SPGwithoutCHMPstopcriteria(A,p,x0,tol,maxit, M, sigmin,sigmax);
            averagetimeProj(j,r) = toc;
            Deltaprojmatrix(j,r) = Deltaproj; 
            end
            
      %% Algorithms 
            %Triangle Algorithm with Random Pivots
            tic, [vta,~,Decision,iterationstav] = TriangleAlgorithm(A,p,p0,tol,maxit);
            timeTAv = toc;
            if Decision == answer    
                    averagetimeTA(j,r) = timeTAv;
                    averageiterationTA(j,r) = iterationstav ;
                    averagenpkpTA(j,r) = norm(vta,'fro');
            else
                countfailTA(j,r) =  1;  %If TA reached the maxit or wrong decision 
                if Decision == -1*answer %If TA wrong decision
                    countfaildecTA(j,r)=0;
                end
            end      

            % Greedy Triangle Algorithm 
            tic, [vgt,~,Decision,iterationsGTv] = GreedyTriangleAlgorithm(A,p,p0,tol,maxit);
            timeGTv = toc;
            if Decision == answer    
                    averagetimeGT(j,r)  = timeGTv;
                    averageiterationGT(j,r) = iterationsGTv ;
                    averagenpkpGT(j,r) = norm(vgt,'fro');
            else
                countfailGT(j,r) =  1; %If GT reached the maxit or wrong decision 
                if Decision == -1*answer %If GT wrong decision
                    countfaildecGT(j,r)=0;
                end
            end      

            %Away Step Frank Wolfe Algorithm          
            tic,[vasfw,~,Decision,iterationsASFWv, countas] = AwayStepFrankWolfeAlgorithm(A,p,p0,i, tol,maxit);
            timeASFWv = toc;
            if Decision == answer    
               averagetimeASFW(j,r)  = timeASFWv;
               averageiterationASFW(j,r) = iterationsASFWv ;
               averagecounterawaystep(j,r) = countas;
               averagenpkpASFW(j,r) = norm(vasfw,'fro');
            else
                 countfailASFW(j,r) =  1;  %If ASFW reached the maxit or wrong decision 
                if Decision == -1*answer %If ASFW wrong decision
                    countfaildecASFW(j,r)=0;
                end
            end      

            %Spectral Projected Gradient  (SPG)
            tic,[Decision,iterationsSPGv,vspg] = SpectralProjectedGradient(A,p,x0,tol,maxit, M, L,sigmin,sigmax);
            timeSPGv = toc;
            if Decision == answer    
               averagetimeSPG(j,r)  = timeSPGv;
               averageiterationSPG(j,r) = iterationsSPGv ;
               averagenpkpSPG(j,r) = norm(vspg,'fro');
            else
                 countfailSPG(j,r) =  1;  %If SPG reached the maxit or linear search failed or wrong decision 
                if Decision == -1*answer %If SPG wrong decision
                    countfaildecSPG(j,r)=0;
                end
            end                   
                  fprintf('\n ..................... \n \n');               
   end
end
        %Aux Labels
        auxlabelTA = zeros(1,nend); 
        for i=1:nend
           [~,b] = min(averagenpkpTA(i,:));
           auxlabelTA(1,i) = b-1;
           [~,b] = min(averagenpkpSPG(i,:));
           auxlabelSPG(1,i) = b-1;
           [~,b] = min(averagenpkpASFW(i,:));
           auxlabelASFW(1,i) = b-1;
           [~,b] = min(averagenpkpGT(i,:));
           auxlabelGT(1,i) = b-1;
           [~,b] = min(Deltaprojmatrix(i,:));
           auxlabelProj(1,i) = b-1;
        end
        %Aux Acuracia
            [a]=find(labels.TA ~= labelsdatatest);
            AcuraciaTA= 100- size(a,2)/100;
            [a]=find(labels.ASFW ~= labelsdatatest);
            AcuraciaASFW= 100- size(a,2)/100;
            [a]=find(labels.GT ~= labelsdatatest);
            AcuraciaGT= 100- size(a,2)/100;
            [a]=find(labels.SPG ~= labelsdatatest);
            AcuraciaSPG= 100- size(a,2)/100;
            [a]=find(labels.Proj ~= labelsdatatest);
            AcuraciaProj= 100- size(a,2)/100;
            
            %% Structs
           %times:
           averagetime.TA = averagetimeTA;  
           averagetime.SPG = averagetimeSPG;
           averagetime.GT =  averagetimeGT; 
           averagetime.ASFW = averagetimeASFW;
           averagetime.PROJ = averagetimeProj;

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
           norms.Proj = Deltaprojmatrix;

           %parameters
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

           %labels
           labels.TA = auxlabelTA;
           labels.GT = auxlabelGT;
           labels.SPG = auxlabelSPG;
           labels.ASFW = auxlabelASFW;
           labels.Proj = auxlabelProj;           
           
           %Acuracia
           acuracia.TA = AcuraciaTA;
           acuracia.GT = AcuraciaGT;
           acuracia.SPG = AcuraciaSPG;
           acuracia.ASFW = AcuraciaASFW;
           acuracia.Proj = AcuraciaProj;
           
            fprintf(' ===================== \n');
end