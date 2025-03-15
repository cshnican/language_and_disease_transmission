function dydt = four_population_spatiotemporal_ODE(t, y, q, param)
    S1 = y(1); % susceptible individuals in population 1
    E1 = y(2); % exposed individuals in population 1
    I1 = y(3); % infectious individuals in population 1
    R1 = y(4); % removed individuals in population 1

    S2 = y(5); % susceptible individuals in population 2
    E2 = y(6); % exposed individuals in population 2
    I2 = y(7); % infectious individuals in population 2
    R2 = y(8); % removed individuals in population 2

    S3 = y(9); % susceptible individuals in population 3
    E3 = y(10); % exposed individuals in population 3
    I3 = y(11); % infectious individuals in population 3
    R3 = y(12); % removed individuals in population 3

    S4 = y(13); % susceptible individuals in population 4
    E4 = y(14); % exposed individuals in population 4
    I4 = y(15); % infectious individuals in population 4
    R4 = y(16); % removed individuals in population 4

    % total number in each population
    N1 = S1 + E1 + I1 + R1;
    N2 = S2 + E2 + I2 + R2;
    N3 = S3 + E3 + I3 + R3;
    N4 = S4 + E4 + I4 + R4;


    % parameters
    beta1 = param(1); % transmission rate within population 1
    beta2 = param(2); % transmission rate within population 2
    beta3 = param(3); % transmission rate within population 3
    beta4 = param(4); % transmission rate within population 4

    q0 = param(5); % contact reduction rate by isolation
    eps = param(6); % infectivity reduction rate for exposed group
    alpha = param(7); % recovery rate
    kappa = param(8); % progression rate to infections
    f = param(9); % survivability

    ts = 0.25; % stating time of commuting in each day (6:00)
    te = 0.75; % ending time of commuting in each day (18:00)
    T = te - ts; % total commuting duration per day(days)
    delta_t = 1; % unit time (1 day)


    % time-dependent term in the contact matrix
    qt = sin(2*pi/T*mod(t - ts, delta_t)) * (mod(t, delta_t) >= ts) * (mod(t, delta_t) <= te);
    
    % contact matrix
    Q = (pi* q) / T * qt;

    % population 1
    dydt(1) = Q(2,1) * S2/N2 + Q(3,1) * S3/N3 + Q(4,1) * S4/N4 - (Q(1,2) + Q(1,3) + Q(1,4)) * S1/N1 - beta1 * S1/N1 * (eps * E1 + (1-q0) * I1);

    dydt(2) = Q(2,1) * E2/N2 + Q(3,1) * E3/N3 + Q(4,1) * E4/N4 - (Q(1,2) + Q(1,3) + Q(1,4)) * E1/N1 + beta1 * S1/N1 * (eps * E1 + (1-q0) * I1) - kappa * E1;

    dydt(3) = Q(2,1) * I2/N2 + Q(3,1) * I3/N3 + Q(4,1) * I4/N4 - (Q(1,2) + Q(1,3) + Q(1,4)) * I1/N1 + kappa * E1 - alpha * I1;

    dydt(4) = Q(2,1) * R2/N2 + Q(3,1) * R3/N3 + Q(4,1) * R4/N4 - (Q(1,2) + Q(1,3) + Q(1,4)) * R1/N1 + f * alpha * I1;

    % population 2
    dydt(5) = Q(1,2) * S1/N1 + Q(3,2) * S3/N3 + Q(4,2) * S4/N4 - (Q(2,1) + Q(2,3) + Q(2,4)) * S2/N2 - beta2 * S2/N2 * (eps * E2 + (1-q0) * I2);

    dydt(6) = Q(1,2) * E1/N1 + Q(3,2) * E3/N3 + Q(4,2) * E4/N4 - (Q(2,1) + Q(2,3) + Q(2,4)) * E2/N2 + beta2 * S2/N2 * (eps * E2 + (1-q0) * I2) - kappa * E2;

    dydt(7) = Q(1,2) * I1/N1 + Q(3,2) * I3/N3 + Q(4,2) * I4/N4 - (Q(2,1) + Q(2,3) + Q(2,4)) * I2/N2 + kappa * E2 - alpha * I2;

    dydt(8) = Q(1,2) * R1/N1 + Q(3,2) * R3/N3 + Q(4,2) * R4/N4 - (Q(2,1) + Q(2,3) + Q(2,4)) * R2/N2 + f * alpha * I2;

    % population 3
    dydt(9) = Q(1,3) * S1/N1 + Q(2,3) * S2/N2 + Q(4,3) * S4/N4 - (Q(3,1) + Q(3,2) + Q(3,4)) * S3/N3 - beta3 * S3/N3 * (eps * E3 + (1-q0) * I3);

    dydt(10) = Q(1,3) * E1/N1 + Q(2,3) * E2/N2 + Q(4,3) * E4/N4 - (Q(3,1) + Q(3,2) + Q(3,4)) * E3/N3 + beta3 * S3/N3 * (eps * E3 + (1-q0) * I3) - kappa * E3;

    dydt(11) = Q(1,3) * I1/N1 + Q(2,3) * I2/N2 + Q(4,3) * I4/N4 - (Q(3,1) + Q(3,2) + Q(3,4)) * I3/N3 + kappa * E3 - alpha * I3;

    dydt(12) = Q(1,3) * R1/N1 + Q(2,3) * R2/N2 + Q(4,3) * R4/N4 - (Q(3,1) + Q(3,2) + Q(3,4)) * R3/N3 + f * alpha * I3;

    % population 4
    dydt(13) = Q(1,4) * S1/N1 + Q(2,4) * S2/N2 + Q(3,4) * S3/N3 - (Q(4,1) + Q(4,2) + Q(4,3)) * S4/N4 - beta4 * S4/N4 * (eps * E4 + (1-q0) * I4);

    dydt(14) = Q(1,4) * E1/N1 + Q(2,4) * E2/N2 + Q(3,4) * E3/N3 - (Q(4,1) + Q(4,2) + Q(4,3)) * E4/N4 + beta4 * S4/N4 * (eps * E4 + (1-q0) * I4) - kappa * E4;

    dydt(15) = Q(1,4) * I1/N1 + Q(2,4) * I2/N2 + Q(3,4) * I3/N3 - (Q(4,1) + Q(4,2) + Q(4,3)) * I4/N4 + kappa * E4 - alpha * I4;

    dydt(16) = Q(1,4) * R1/N1 + Q(2,4) * R2/N2 + Q(3,4) * R3/N3 - (Q(4,1) + Q(4,2) + Q(4,3)) * R4/N4 + f * alpha * I4;

    dydt = dydt';


    