function SDEPMOEA(Global)
% <algorithm> <S>
% SDE^+ - MOEA

%% Generate the sampling points and random population
Population = Global.Initialization();

%% Optimization
while Global.NotTermination(Population)
    [~,SDEP_fitness,~] = SDE_plus_indicator(Population.objs, 1);          
    MatingPool = TournamentSelection(2,Global.N,-SDEP_fitness);  
    Offspring  = GA(Population(MatingPool),{1,20,1,20}); 
    Population = SDEMOEA_EnvironmentalSelection([Population,Offspring],Global.N);
end

end

function [SDE_fitness,SDEP_fitness, Distance] = SDE_plus_indicator(PopObj, f)
    N = size(PopObj,1);
    PopObj = (PopObj-repmat(min(PopObj,[],1),N,1))./repmat(max(PopObj,[],1)-min(PopObj,[],1),N,1);
    Distance = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = 1 : N
            if(i == j)
                continue;
            end
            Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    SDEP_fitness = zeros(1, N);
    if f == 1
        [a,Remain] = sortrows(sort(Distance, 2),'descend');
        SDEP_fitness(Remain) = N:-1:1;
    end
    SDE_fitness = min(Distance,[],2); 
end

function Population = SDEMOEA_EnvironmentalSelection(Population,N)

    [SDE_fitness, ~, SDE_Distance] = SDE_plus_indicator(Population.objs, 0);
    if(sum(SDE_fitness ~= 0) <= N)
        %  One-time evaluation
        [~,Remain] = sortrows(sort(SDE_Distance, 2),'descend');
        Remain = Remain(1:N);
    else
        % Dynamic evaluation
        Choose = SDE_fitness ~= 0;
        while sum(Choose) ~= N
            Remain = find(Choose);
            mat = min(SDE_Distance(Choose,Choose),2);
            [~,x]  = min(min(mat,[],2));
            Choose(Remain(x)) = false;
        end
    end
    Population = Population(Remain);
end
