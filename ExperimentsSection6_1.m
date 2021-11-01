%Experiments_section_6_1
Titles = ['Case a'; 'Case b'; 'Case c'; 'Case d'];
Size = ['pequeno'; 'grande '];
for pvar = 1:4
    for nvar = 1:2 
        [averagetime, averageiteration, parameters, countfail] = artificialchmp(nvar,pvar);
        [n] = generatesgraphics(pvar,nvar,averagetime,countfail);
        title(Titles(pvar,:));
        hold off
        save ([Titles(pvar,:),Size(nvar,:)])
        saveas(gcf,[Titles(pvar,:),Size(nvar,:)])    
    end 
end
