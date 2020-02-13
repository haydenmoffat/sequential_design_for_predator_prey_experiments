% Samples N particles from the prior distributions of each model
% specified. Here we consider a log normal prior on each of the parameters where 
% the mean and standard deviation of the natural logarithm of the variables 
% are -1.4 and 1.35, respectively.

theta = zeros(N, 3, K);
for p = 1:K
    if ismember(Models(p),[1 2])
        theta(:,:,p) = normrnd(-1.4, 1.35, N, 3);
    
    elseif ismember(Models(p),[3 4])
        theta(:,:,p) = [normrnd(-1.4, 1.35, N, 1), log(zeros(N,1)),normrnd(-1.4, 1.35, N, 1)];
    
    elseif ismember(Models(p),[5 6])
        theta(:,:,p) = [normrnd(-1.4, 1.35, N, 2), log(zeros(N,1))];
   
    elseif ismember(Models(p), [7 8])
        theta(:,:,p) = [normrnd(-1.4, 1.35, N, 1), log(zeros(N,2))];
    end 
end
