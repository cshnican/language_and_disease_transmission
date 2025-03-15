function dydt = trudgill_network_ODE(t, y, q, param)
    
    % input:
    % t: time range
    % y: 1 x 4n vector (n = number of populations)
    % q: # of people who commute between populations
    % param: parameters for disease transmission, etc.

    % number of populations
    num_populations = length(y) / 4;

    % reshape vectors:
    y = reshape(y, 4, num_populations); % each column represents a population

    S = y(1, :); % susceptable
    E = y(2, :); % exposed
    I = y(3, :); % infectious
    R = y(4, :); % recovered

    % total population in each group
    N = S + E + I + R;

    % parameters
    beta = param(1: num_populations); % transmission rate
    q0 = param(num_populations + 1); % contact reduction rate
    eps = param(num_populations + 2); % infectivity reduction for exposed group
    alpha = param(num_populations + 3); % recovery rate
    kappa = param(num_populations + 4); % progression rate
    f = param(num_populations + 5); % survivability


    % commuting times (unit: day)
    ts = 0.25; % start commuting (6:00)
    te = 0.75; % end commuting (18:00)
    T = te - ts; % commuting durarion
    delta_t = 1; % unit time

    % time-dependent term:
    qt = sin(2*pi/T*mod(t - ts, delta_t)) * (mod(t, delta_t) >= ts) * (mod(t, delta_t) <= te);

    % contact matrix
    Q =(pi * q) / T * qt;

    % initialize derivatives
    dS = zeros(1, num_populations);
    dE = zeros(1, num_populations);
    dI = zeros(1, num_populations);
    dR = zeros(1, num_populations);

    % compute derivatives
    for i = 1:num_populations

        % susceptible
        dS(i) = sum(Q(:, i)' .* S ./ N) - sum(Q(i, :) .* S(i) ./ N(i)) - beta(i) * S(i) / N(i) * (eps * E(i) + (1 - q0) * I(i));

        % Exposed
        dE(i) = sum(Q(:, i)' .* E ./ N) - sum(Q(i, :) .* E(i) ./ N(i)) + beta(i) * S(i) / N(i) * (eps * E(i) + (1 - q0) * I(i)) - kappa * E(i);

        % infectious
        dI(i) = sum(Q(:, i)' .* I ./ N) - sum(Q(i, :) .* I(i) ./ N(i)) + kappa * E(i) - alpha * I(i);

        % recovered
        dR(i) = sum(Q(:, i)' .* R ./ N) - sum(Q(i, :) .* R(i) ./ N(i)) + f * alpha * I(i);
    end

    % combine
    dydt = [dS; dE; dI; dR];
    dydt = dydt(:); % flatten