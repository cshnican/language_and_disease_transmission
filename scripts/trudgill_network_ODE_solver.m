function sol = trudgill_network_ODE_solver(y0, q, tspan)
    n_populations = length(y0)/4;
   
    % parameters
    betas = 0.56 * ones(1, n_populations);
    
    q0 = 0; % contact reduction rate by isolation
    eps = 0; % infectivity reduction rate for exposed group
    alpha = 1/6; % recovery rate
    kappa = 1/1.9; % progression rate to infections
    f = 1; % survivability

    param = [betas, q0, eps, alpha, kappa, f];

    options = odeset('RelTol',1e-6);

    sol = ode45(@(t,y)trudgill_network_ODE(t,y,q, param),tspan,y0,options);