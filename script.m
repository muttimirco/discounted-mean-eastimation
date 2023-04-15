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
horizon = 10;
batches = [1:50:1001]; % batches = [1:50:1000];

% seeding
n_seeds = 20;
seeds = [1:n_seeds];
curves = zeros(4, length(batches));
confs = zeros(4, length(batches));
errs = zeros(n_seeds, length(batches));

% all sample
for s = seeds
    rng(s);
    estimates = [];
    for n_trajectories = batches
       estimates(end + 1) = all_sample_per_trajectory(mumu, PP, gamma, ff, n_trajectories, horizon * n_trajectories);
    end
    errors = abs(estimates - med);
    errs(s, :) = errors;
end
curves(1, :) = mean(errs);
confs(1, :) = 2 * std(errs) / sqrt(length(seeds));

% one sample
for s = seeds
    rng(s);
    estimates = [];
    for n_trajectories = batches
       estimates(end + 1) = one_sample_per_trajectory(mumu, PP, gamma, ff, n_trajectories, horizon * n_trajectories);
    end
    errors = abs(estimates - med);
    errs(s, :) = errors;
end
curves(2, :) = mean(errs);
confs(2, :) = 2 * std(errs) / sqrt(length(seeds));

% corrected
for s = seeds
    rng(s);
    estimates = [];
    for n_trajectories = batches
       estimates(end + 1) = finite_horizon_corrected(mumu, PP, gamma, ff, n_trajectories, horizon);
    end
    errors = abs(estimates - med);
    errs(s, :) = errors;
end
curves(3, :) = mean(errs);
confs(3, :) = 2 * std(errs) / sqrt(length(seeds));

% non-corrected
for s = seeds
    rng(s);
    estimates = [];
    for n_trajectories = batches
       estimates(end + 1) = finite_horizon_non_corrected(mumu, PP, gamma, ff, n_trajectories, horizon);
    end
    errors = abs(estimates - med);
    errs(s, :) = errors;
end
curves(4, :) = mean(errs);
confs(4, :) = 2 * std(errs) / sqrt(length(seeds));

% figure
fig = figure();
hold on
errorbar(batches * horizon, curves(1, :), confs(1, :));
errorbar(batches * horizon, curves(2, :), confs(2, :));
errorbar(batches * horizon, curves(3, :), confs(3, :));
errorbar(batches * horizon, curves(4, :), confs(4, :));
ylabel('error');
xlabel('N');
legend('all-sample', 'one-sample', 'corrected', 'non-corrected');