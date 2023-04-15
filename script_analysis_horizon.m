clc
clear
close all

% Markov chain
alpha = 0.99;
P = [[alpha, 1 - alpha, 0.]; 
    [(1 - alpha) / 2, alpha, (1 - alpha) / 2]; 
    [0., 1 - alpha, alpha]];
S = size(P, 1);
mu = [1., 0., 0.];
f = [1, -1, 2]';
gamma = 0.99;

[mumu, PP, ff] = build_generators(mu, P, f);

% true value
pi_gamma = (1 - gamma) * mu / (eye(S) - gamma * P);
med = pi_gamma * f;
eigenvalues = eig(P);
eigenvalues = abs(eigenvalues(1:end-1));
beta = max(eigenvalues);


% estimation
batches = 1:100:1000;

% seeding
n_seeds = 20;
seeds = [1:n_seeds];

% 10
horizon_1 = 10;
batches_1 = 1:50:1000;
curves_1 = zeros(1, length(batches_1));
confs_1 = zeros(1, length(batches_1));
errs = zeros(n_seeds, length(batches_1));
for s = seeds
    rng(s);
    estimates = [];
    for n_trajectories = batches_1
       estimates(end + 1) = finite_horizon_non_corrected(mumu, PP, gamma, ff, n_trajectories, horizon_1);
    end
    errors = abs(estimates - med);
    errs(s, :) = errors;
end
curves_1(1, :) = mean(errs);
confs_1(1, :) = 2 * std(errs) / sqrt(length(seeds));

% 30
horizon_2 = 30;
batches_2 = 1:15:330;
curves_2 = zeros(1, length(batches_2));
confs_2 = zeros(1, length(batches_2));
errs = zeros(n_seeds, length(batches_2));
for s = seeds
    rng(s);
    estimates = [];
    for n_trajectories = batches_2
       estimates(end + 1) = finite_horizon_non_corrected(mumu, PP, gamma, ff, n_trajectories, horizon_2);
    end
    errors = abs(estimates - med);
    errs(s, :) = errors;
end
curves_2(1, :) = mean(errs);
confs_2(1, :) = 2 * std(errs) / sqrt(length(seeds));

% 50
horizon_3 = 50;
batches_3 = 1:10:200;
curves_3 = zeros(1, length(batches_3));
confs_3 = zeros(1, length(batches_3));
errs = zeros(n_seeds, length(batches_3));
for s = seeds
    rng(s);
    estimates = [];
    for n_trajectories = batches_3
       estimates(end + 1) = finite_horizon_non_corrected(mumu, PP, gamma, ff, n_trajectories, horizon_3);
    end
    errors = abs(estimates - med);
    errs(s, :) = errors;
end
curves_3(1, :) = mean(errs);
confs_3(1, :) = 2 * std(errs) / sqrt(length(seeds));

fig = figure();
hold on
errorbar(batches_1 * horizon_1, curves_1, confs_1);
errorbar(batches_2 * horizon_2, curves_2, confs_2);
errorbar(batches_3 * horizon_3, curves_3, confs_3);
ylabel('error');
xlabel('N');
legend('h=10', 'h=30', 'h=50');