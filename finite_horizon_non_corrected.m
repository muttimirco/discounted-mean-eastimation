function [mm, ss, N, M, Hs] = finite_horizon_non_corrected(mu, P, gamma, f, n_trajectories, horizon)
    total_sum = 0;
    for ii = 1:n_trajectories
        cum_disc_sum = 0;
        x = mu();
        y = f(x);
        for jj = 0:horizon-1
            cum_disc_sum = cum_disc_sum + gamma^jj * y;
            x = P(x);
            y = f(x);
        end
        total_sum = total_sum + cum_disc_sum;
    end
    mm = (1 - gamma) * total_sum / n_trajectories;
    ss = 0;
    N = n_trajectories * horizon;
    M = n_trajectories;
    Hs = horizon * ones(n_trajectories);
end