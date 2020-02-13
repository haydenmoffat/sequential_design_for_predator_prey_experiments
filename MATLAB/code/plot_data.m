% Plotting example datasets. The sample means (o-o-o)
% for each value N are also shown.

% Select dataset to use
D = 1;
% Datasets:
% 1- Example dataset from Papanikolau et al. (2016)
% 2- The orginal data from Hassell et al. (1977)

% Load the dataset
if D == 1
    raw_data = load('papanikolau_data.mat');
    data = raw_data.data;
    
elseif D == 2
    raw_data = load('hassell_data.mat');
    data = [raw_data.N.' raw_data.n.'];
end

% Plot the data with jitter on the y axis
scatter(data(:,1), data(:,2) + normrnd(0, 0.1, size(data, 1), 1), 60, '.');
hold on;
grid on;
means = grpstats(data, data(:,1));
plot(means(:,1),means(:,2), '-o');
xlabel('Number of prey available (N)');
ylabel('Number predated (n)');
legend('Individual samples', 'Means for each value of N', 'Location', 'northwest')

