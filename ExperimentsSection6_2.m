%Experiments_section_6_2
addpath ('Functions')
addpath ('Algorithms')
m = [50; 50; 100;200];
n = [200; 2000; 500 ;2000];
tolfix = [1e-6; 1e-7];
viability = [1;-1];
Titlesm = ['m_50__'; 'm_50__'; 'm_100_'; 'm_200_'];
Titlesn = ['n_200__'; 'n_2000_'; 'n_500__'; 'n_2000_'];
Titlesviability = ['feasible__'; 'infeasible'];
Titletol = ['_1e-6';'_1e-7'];
mkdir visitordata6-2;
for i= 1:4  
    for j = 1:2 
        for k = 1:2
            [averagetime, averageiteration, parameters, countfail] =  artificialfeasibilityproblem (m(i),n(i),viability(j), tolfix(k));  
            namesave = fullfile('visitordata6-2',[Titlesm(i,:),Titlesn(i,:),Titlesviability(j,:),Titletol(k,:)]);
            save (namesave);
        end       
    end
end
