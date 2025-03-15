function sol = four_population_spatiotemporal_ODE_solver(y0, q, tspan)

    % parameters
    beta1 = 0.56; % transmission rate within population 1
    beta2 = 0.56; % transmission rate within population 2
    beta3 = 0.56; % transmission rate within population 3
    beta4 = 0.56; % transmission rate within population 4

    q0 = 0; % contact reduction rate by isolation
    eps = 0; % infectivity reduction rate for exposed group
    alpha = 1/6; % recovery rate
    kappa = 1/1.9; % progression rate to infections
    f = 1; % survivability

    param = [beta1, beta2, beta3, beta4, q0, eps, alpha, kappa, f];


    % commuting matrix
    % q(1,:) = [0 400 200 100];
    % q(2,:) = [100 0 300 400];
    % q(3,:) = [100 600 0 200];
    % q(4,:) = [500 200 400 0];

    %tspan = [0 1];

    %y0 = [10000 0 0 0 1000 0 0 0 1000 0 0 0 3000 0 0 0];

    options = odeset('RelTol',1e-9);

    sol = ode45(@(t,y)four_population_spatiotemporal_ODE(t,y,q, param),tspan,y0,options);

    t = sol.x;
    y = sol.y;

    % extract variables
    S1 = y(1,:);
    E1 = y(2,:);
    I1 = y(3,:);
    R1 = y(4,:);

    S2 = y(5,:);
    E2 = y(6,:);
    I2 = y(7,:);
    R2 = y(8,:);

    S3 = y(9,:);
    E3 = y(10,:);
    I3 = y(11,:);
    R3 = y(12,:);

    S4 = y(13,:);
    E4 = y(14,:);
    I4 = y(15,:);
    R4 = y(16,:);

    % total populations
    N1 = S1 + E1 + I1 + R1;
    N2 = S2 + E2 + I2 + R2;
    N3 = S3 + E3 + I3 + R3;
    N4 = S4 + E4 + I4 + R4;