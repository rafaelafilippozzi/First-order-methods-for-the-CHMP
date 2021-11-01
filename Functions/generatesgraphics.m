 function [n] = generatesgraphics(pvar,n,averagetime,countfail)

if n==1
 n=500:500:5000;
else
 n=10000:10000:100000;
end
legendaux = [];
%Tempo
if countfail.allTA == zeros(10,1) 
   semilogy( n, averagetime.TA (:,3), 'k-.', 'LineWidth', 3)
   legendaux = [legendaux ;'TA  '];
   hold on
end
if countfail.allGT == zeros(10,1)
   semilogy( n, averagetime.GT (:,3),'k--', 'LineWidth',3)
   legendaux = [legendaux ;'GT  '];
   hold on
end
if countfail.allSPG == zeros(10,1)
   semilogy( n, averagetime.SPG(:,3), 'k-', 'LineWidth', 3)
   legendaux = [legendaux ; 'SPG '];
   hold on
end
if countfail.allASFW == zeros(10,1)
   semilogy( n, averagetime.ASFW(:,3), 'k:', 'LineWidth',3)
   legendaux = [legendaux ; 'ASFW'];
   hold on
end 
xlabel ('n')
ylabel ('time(s)')
[m, ~] = size (legendaux);
if m == 4
legend(legendaux(1,:),legendaux(2,:),legendaux(3,:),legendaux(4,:))
else if m == 3
        legend(legendaux(1,:),legendaux(2,:),legendaux(3,:))
     else if m == 2
          legend(legendaux(1,:),legendaux(2,:))
          else 
           legend(legendaux(1,:))
          end
     end
end


