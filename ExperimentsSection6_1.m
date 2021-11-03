%Experiments_section_6_1
addpath ('Functions')
addpath ('Algorithms')
Titles = ['Case a'; 'Case b'; 'Case c'; 'Case d'];
Size = ['pequeno'; 'grande '];
mkdir visitordata6-1;
for pvar = 1:4
    for nvar = 1:2                
        [averagetime, averageiteration, parameters, countfail] = artificialchmp(nvar,pvar);
        [n] = generatesgraphics(pvar,nvar,averagetime,countfail);
        title(Titles(pvar,:));
        hold off
        namesave = fullfile('visitordata6-1',[Titles(pvar,:),Size(nvar,:)]);
        save (namesave);
        saveas(gcf,namesave)    
    end 
end
