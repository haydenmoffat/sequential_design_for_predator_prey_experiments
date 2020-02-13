% Solves the ode for a range of parameters and initial conditions for
% a given model and time
% theta is the current particle values
% init is the initial prey population
% time is the time that the predator has access to the prey
% M is the model of interest

function n = ode_solver(theta, init, time, M)
tdom= [0 time];
p = [theta(:,1,:), theta(:,2,:)];
n=zeros(size(theta,1),length(init));
for f = 1:size(theta,1)
    for g=1:length(init)
        if ismember(M, [1 5])
            [~,u] = ode45(@(tt,u) ode_M2(tt,u,p(f,:)), tdom, init(g));
        elseif ismember(M, [2 6])
            [~,u] = ode45(@(tt,u) ode_M3(tt,u,p(f,:)), tdom, init(g));
        elseif ismember(M, [3 7])
            u = init(g)*exp(-p(f,1)*time);
        elseif ismember(M, [4 8])
            u = init(g) / (init(g)*p(f,1)*time + 1);
        end
        n(f,g) = u(length(u));
        if n(f,g) < 0
            n(f,g) = 0;
        end
    end
    
end
end
