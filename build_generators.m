function [mumu, PP, ff] = build_generators(mu, P, f)
    PP = @(x) mnrnd(1, x * P);
    mumu = @() mnrnd(1, mu);
    ff = @(x) x * f;
end