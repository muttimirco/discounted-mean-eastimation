function [mm, ss, N, M, Hs] = one_sample_per_trajectory(mu, P, gamma, f, n_trajectories, n_samples)
    total_sum = 0;
    N = 0;
    Hs = zeros(n_trajectories);
    for ii = 1:n_trajectories
        horizon = geornd(1 - gamma);
        Hs(ii) = horizon;
        N = N + horizon;
        if N > n_samples
            mm = total_sum / ii;
            return
        end
        x = mu();
        for jj = 0:horizon-1
            x = P(x);
        end
        y = f(x);
        total_sum = total_sum + y;
    end
    mm = total_sum / n_trajectories;
    ss = 0;
    M = n_trajectories;
end