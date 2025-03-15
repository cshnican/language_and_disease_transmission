clear all
close all
clc


% terminology:
% population - groups of people speaking the same language
% subgroup - a group in a population that live closely together

rng(69, "twister");
[T1, T2, T3, T4, T5, T6, half_points] = run_scenarios_paper();
writetable(T1, '../output_table/scenario1.txt');
writetable(T2, '../output_table/scenario2.txt');
writetable(T3, '../output_table/scenario3.txt');
writetable(T4, '../output_table/scenario4.txt');
writetable(T5, '../output_table/scenario5.txt');
writetable(T6, '../output_table/scenario6.txt');
writetable(half_points, '../output_table/half_points_scenario1-6.txt');

function [T1, T2, T3, T4, T5, T6, half_points] = run_scenarios_paper()
    low_rates = 0.1:0.01:0.15;
    high_rates = 1:0.01:1.5; 

    % Scenario 1
    populations = [1000, 1000, 1000, 1000];
    ngroups = 4;
    
    sol1 = make_trudgill_network(populations, ngroups, high_rates, low_rates);
    half_points(1, :) = get_half_point(sol1, populations, ngroups);

    T1 = array2table([sol1.x',sol1.y']);
    T1.Properties.VariableNames = {'t', 'S1', 'E1', 'I1', 'R1',...
        'S2', 'E2', 'I2', 'R2',...
        'S3', 'E3', 'I3', 'R3',...
        'S4', 'E4', 'I4', 'R4',...
        'S5', 'E5', 'I5', 'R5',...
        'S6', 'E6', 'I6', 'R6',...
        'S7', 'E7', 'I7', 'R7',...
        'S8', 'E8', 'I8', 'R8',...
        'S9', 'E9', 'I9', 'R9',...
        'S10', 'E10', 'I10', 'R10',...
        'S11', 'E11', 'I11', 'R11',...
        'S12', 'E12', 'I12', 'R12',...
        'S13', 'E13', 'I13', 'R13',...
        'S14', 'E14', 'I14', 'R14',...
        'S15', 'E15', 'I15', 'R15',...
        'S16', 'E16', 'I16', 'R16'};


    % Scenario 2
    populations = [1000, 1000, 1000, 1000];
    ngroups = 4;
    
    sol2 = make_trudgill_network(populations, ngroups, high_rates, high_rates);
    half_points(2, :) = get_half_point(sol2, populations, ngroups);
    
    T2 = array2table([sol2.x',sol2.y']);
    T2.Properties.VariableNames = {'t', 'S1', 'E1', 'I1', 'R1',...
        'S2', 'E2', 'I2', 'R2',...
        'S3', 'E3', 'I3', 'R3',...
        'S4', 'E4', 'I4', 'R4',...
        'S5', 'E5', 'I5', 'R5',...
        'S6', 'E6', 'I6', 'R6',...
        'S7', 'E7', 'I7', 'R7',...
        'S8', 'E8', 'I8', 'R8',...
        'S9', 'E9', 'I9', 'R9',...
        'S10', 'E10', 'I10', 'R10',...
        'S11', 'E11', 'I11', 'R11',...
        'S12', 'E12', 'I12', 'R12',...
        'S13', 'E13', 'I13', 'R13',...
        'S14', 'E14', 'I14', 'R14',...
        'S15', 'E15', 'I15', 'R15',...
        'S16', 'E16', 'I16', 'R16'};

    % Scenario 3 / Trudgill Type 3
    populations = [1000, 1000, 1000, 1000];
    ngroups = 4;
    
    sol3 = make_trudgill_network(populations, ngroups, low_rates, low_rates);
    half_points(3, :) = get_half_point(sol3, populations, ngroups);
    
    T3 = array2table([sol3.x',sol3.y']);
    T3.Properties.VariableNames = {'t', 'S1', 'E1', 'I1', 'R1',...
        'S2', 'E2', 'I2', 'R2',...
        'S3', 'E3', 'I3', 'R3',...
        'S4', 'E4', 'I4', 'R4',...
        'S5', 'E5', 'I5', 'R5',...
        'S6', 'E6', 'I6', 'R6',...
        'S7', 'E7', 'I7', 'R7',...
        'S8', 'E8', 'I8', 'R8',...
        'S9', 'E9', 'I9', 'R9',...
        'S10', 'E10', 'I10', 'R10',...
        'S11', 'E11', 'I11', 'R11',...
        'S12', 'E12', 'I12', 'R12',...
        'S13', 'E13', 'I13', 'R13',...
        'S14', 'E14', 'I14', 'R14',...
        'S15', 'E15', 'I15', 'R15',...
        'S16', 'E16', 'I16', 'R16'};

    % Scenario 4 / Trudgill Type 4
    populations = [1000, 1000, 1000, 1000];
    ngroups = 4;
    
    sol4 = make_trudgill_network(populations, ngroups, low_rates, high_rates);
    half_points(4, :) = get_half_point(sol4, populations, ngroups);
    
    T4 = array2table([sol4.x',sol4.y']);
    T4.Properties.VariableNames = {'t', 'S1', 'E1', 'I1', 'R1',...
        'S2', 'E2', 'I2', 'R2',...
        'S3', 'E3', 'I3', 'R3',...
        'S4', 'E4', 'I4', 'R4',...
        'S5', 'E5', 'I5', 'R5',...
        'S6', 'E6', 'I6', 'R6',...
        'S7', 'E7', 'I7', 'R7',...
        'S8', 'E8', 'I8', 'R8',...
        'S9', 'E9', 'I9', 'R9',...
        'S10', 'E10', 'I10', 'R10',...
        'S11', 'E11', 'I11', 'R11',...
        'S12', 'E12', 'I12', 'R12',...
        'S13', 'E13', 'I13', 'R13',...
        'S14', 'E14', 'I14', 'R14',...
        'S15', 'E15', 'I15', 'R15',...
        'S16', 'E16', 'I16', 'R16'};

    % Scenario 5 / Trudgill Type 5
    populations = [4000, 4000, 4000, 4000];
    ngroups = 4;
    
    sol5 = make_trudgill_network(populations, ngroups, low_rates, low_rates);
    half_points(5, :) = get_half_point(sol5, populations, ngroups);
    
    T5 = array2table([sol5.x',sol5.y']);
    T5.Properties.VariableNames = {'t', 'S1', 'E1', 'I1', 'R1',...
        'S2', 'E2', 'I2', 'R2',...
        'S3', 'E3', 'I3', 'R3',...
        'S4', 'E4', 'I4', 'R4',...
        'S5', 'E5', 'I5', 'R5',...
        'S6', 'E6', 'I6', 'R6',...
        'S7', 'E7', 'I7', 'R7',...
        'S8', 'E8', 'I8', 'R8',...
        'S9', 'E9', 'I9', 'R9',...
        'S10', 'E10', 'I10', 'R10',...
        'S11', 'E11', 'I11', 'R11',...
        'S12', 'E12', 'I12', 'R12',...
        'S13', 'E13', 'I13', 'R13',...
        'S14', 'E14', 'I14', 'R14',...
        'S15', 'E15', 'I15', 'R15',...
        'S16', 'E16', 'I16', 'R16'};

    % Scenario 6 / Trudgill Type 6
    populations = [4000, 4000, 4000, 4000];
    ngroups = 4;
    
    sol6 = make_trudgill_network(populations, ngroups, low_rates, high_rates);
    half_points(6, :) = get_half_point(sol6, populations, ngroups);
    
    T6 = array2table([sol6.x',sol6.y']);
    T6.Properties.VariableNames = {'t', 'S1', 'E1', 'I1', 'R1',...
        'S2', 'E2', 'I2', 'R2',...
        'S3', 'E3', 'I3', 'R3',...
        'S4', 'E4', 'I4', 'R4',...
        'S5', 'E5', 'I5', 'R5',...
        'S6', 'E6', 'I6', 'R6',...
        'S7', 'E7', 'I7', 'R7',...
        'S8', 'E8', 'I8', 'R8',...
        'S9', 'E9', 'I9', 'R9',...
        'S10', 'E10', 'I10', 'R10',...
        'S11', 'E11', 'I11', 'R11',...
        'S12', 'E12', 'I12', 'R12',...
        'S13', 'E13', 'I13', 'R13',...
        'S14', 'E14', 'I14', 'R14',...
        'S15', 'E15', 'I15', 'R15',...
        'S16', 'E16', 'I16', 'R16'};


    half_points = array2table(half_points');
    half_points.Properties.VariableNames = {'Scenario1', 'Scenario2', 'Scenario3', 'Scenario4', 'Scenario5', 'Scenario6'};
end





function T = run_XS()

low_rates = 0.1:0.01:0.15;
high_rates = 1:0.01:1.5; 

% 2 populations containing 4 subpopulations each ~ X languages used as
% national language

populations = [4000, 4000];
ngroups = 4;

sol = make_trudgill_network_2pop(populations, ngroups, high_rates, high_rates);
half_points(1, :) = get_half_point(sol, populations, ngroups);

plot_things(sol, populations, ngroups)
sgtitle('X languages used as national language')

% 2 populations containing 4 subpopulations each ~ X languages used as
% lingua franca

populations = [4000, 4000];
ngroups = 4;

sol = make_trudgill_network_2pop(populations, ngroups, low_rates, high_rates);
half_points(2, :) = get_half_point(sol, populations, ngroups);

plot_things(sol, populations, ngroups)
sgtitle('X languages used as lingua franca')

% 2 populations containing 4 subpopulations each ~ S languages used in
% populated areas

populations = [4000, 4000];
ngroups = 4;

sol = make_trudgill_network_2pop(populations, ngroups, high_rates, low_rates);
half_points(3, :) = get_half_point(sol, populations, ngroups);

plot_things(sol, populations, ngroups)
sgtitle('S languages in populated areas')

% 2 populations containing 4 subpopulations each ~ S languages used in
% isolated groups

populations = [4000, 4000];
ngroups = 4;

sol = make_trudgill_network_2pop(populations, ngroups, low_rates, low_rates);
half_points(4, :) = get_half_point(sol, populations, ngroups);

plot_things(sol, populations, ngroups)
sgtitle('S languages in isolated groups')


T = array2table(half_points');
T.Properties.VariableNames = {'Scenario1', 'Scenario2', 'Scenario3', 'Scenario4'};

end





function T = run_simple()

% 4 big populations, with no subgroups

    populations = [4000, 4000, 4000, 4000];
    ngroups = 1;

    sol = scenario1();
    half_points_simple(1, :) = get_half_point(sol, populations, ngroups);
    
    
    plot_things(sol, populations(1), ngroups)
    sgtitle('1 big population, 4000 people')
    
    

% 4 groups of non-isolated esoteric populations with size 1000 each

    populations = [1000, 1000, 1000, 1000];
    ngroups = 1;

    sol = scenario2();
    half_points_simple(2, :) = get_half_point(sol, populations, ngroups);

    plot_things(sol, populations, ngroups)
    sgtitle('4 big populations, 1000 people each, high contact')


% 4 groups of isolated esoteric populations with size 1000 each

    populations = [1000, 1000, 1000, 1000];
    ngroups = 1;

    sol = scenario3();
    half_points_simple(3, :) = get_half_point(sol, populations, ngroups);

    plot_things(sol, populations, ngroups)
    sgtitle('4 big populations, 1000 people each, low contact')

    T = array2table(half_points_simple');
    T.Properties.VariableNames = {'Type1_simple', 'Type2_simple', 'Type3_simple'};
end

function T = run_trudgill()

low_rates = 0.1:0.01:0.15;
high_rates = 1:0.01:1.5; 

% tridgill_type_1
populations = [1000, 1000, 1000, 1000];
ngroups = 4;

sol = make_trudgill_network(populations, ngroups, high_rates, low_rates);
half_points(1, :) = get_half_point(sol, populations, ngroups);

plot_things(sol, populations, ngroups)
sgtitle('Scenario 2')

% tridgill_type_2
populations = [1000, 1000, 1000, 1000];
ngroups = 4;

sol = make_trudgill_network(populations, ngroups, high_rates, high_rates);
half_points(2, :) = get_half_point(sol, populations, ngroups);

plot_things(sol, populations, ngroups)
sgtitle('Scenario 1')

% tridgill_type_3
populations = [1000, 1000, 1000, 1000];
ngroups = 4;

sol = make_trudgill_network(populations, ngroups, low_rates, low_rates);
half_points(3, :) = get_half_point(sol, populations, ngroups);

plot_things(sol, populations, ngroups)
sgtitle('Scenario 3')

% tridgill_type_4
populations = [1000, 1000, 1000, 1000];
ngroups = 4;

sol = make_trudgill_network(populations, ngroups, low_rates, high_rates);
half_points(4, :) = get_half_point(sol, populations, ngroups);

plot_things(sol, populations, ngroups)
sgtitle('Scenario 4')

% tridgill_type_5
populations = [4000, 4000, 4000, 4000];
ngroups = 4;

sol = make_trudgill_network(populations, ngroups, low_rates, low_rates);
half_points(5, :) = get_half_point(sol, populations, ngroups);

plot_things(sol, populations, ngroups)
sgtitle('Scenario 5')

% trudgill_type_6
populations = [4000, 4000, 4000, 4000];
ngroups = 4;

sol = make_trudgill_network(populations, ngroups, low_rates, high_rates);
half_points(6, :) = get_half_point(sol, populations, ngroups);

plot_things(sol, populations, ngroups)
sgtitle('Scenario 6')

T = array2table(half_points');
T.Properties.VariableNames = {'Scenario2', 'Scenario1', 'Scenario3', 'Scenario4', 'Scenario5', 'Scenario6'};

end

function sol = scenario1()
    % Scenario 1: 4000 people
    q(1,:) = [0 0 0 0];
    q(2,:) = [0 0 0 0];
    q(3,:) = [0 0 0 0];
    q(4,:) = [0 0 0 0];
    
    y0 = [3999 1 0 0 ...
          3999 1 0 0 ...
          3999 1 0 0 ...
          3999 1 0 0];
    
    tspan = [0 100];
    
    sol = trudgill_network_ODE_solver(y0, q, tspan);
    
end


function sol =  scenario2
    % Scenario 2: 1000 each; non-isolated esoteric groups

    high_rates = 1:0.01:1.5; 
    
    q(1,:) = [0 datasample(high_rates, 1) datasample(high_rates, 1) datasample(high_rates, 1)];
    q(2,:) = [datasample(high_rates, 1) 0 datasample(high_rates, 1) datasample(high_rates, 1)];
    q(3,:) = [datasample(high_rates, 1) datasample(high_rates, 1) 0 datasample(high_rates, 1)];
    q(4,:) = [datasample(high_rates, 1) datasample(high_rates, 1) datasample(high_rates, 1) 0];
    
    
    y0 = [999 1 0 0 ...
          1000 0 0 0 ...
          1000 0 0 0 ...
          1000 0 0 0];
    
    tspan = [0 100];
    
    sol = trudgill_network_ODE_solver(y0, q, tspan);
    
end    


function sol = scenario3


    % Scenario 3: 1000 each; isolated esoteric groups
    
    low_rates = 0.01:0.001:0.02;
    
    q(1,:) = [0 datasample(low_rates, 1) datasample(low_rates, 1) datasample(low_rates, 1)];
    q(2,:) = [datasample(low_rates, 1) 0 datasample(low_rates, 1) datasample(low_rates, 1)];
    q(3,:) = [datasample(low_rates, 1) datasample(low_rates, 1) 0 datasample(low_rates, 1)];
    q(4,:) = [datasample(low_rates, 1) datasample(low_rates, 1) datasample(low_rates, 1) 0];
    
    
    y0 = [999 1 0 0 ...
          1000 0 0 0 ...
          1000 0 0 0 ...
          1000 0 0 0];
    
    tspan = [0 100];
    
    sol = trudgill_network_ODE_solver(y0, q, tspan);
    
end

function sol = make_trudgill_network_2pop(populations, ngroups, within_group_rates, between_group_rates)

    qii = [0, datasample(within_group_rates, 1), datasample(within_group_rates, 1), datasample(within_group_rates, 1);
        datasample(within_group_rates, 1), 0, datasample(within_group_rates, 1), datasample(within_group_rates, 1);
        datasample(within_group_rates, 1), datasample(within_group_rates, 1), 0, datasample(within_group_rates, 1);
        datasample(within_group_rates, 1), datasample(within_group_rates, 1), datasample(within_group_rates, 1), 0];

    qij = [datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1);
        datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1);
        datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1);
        datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1)];

    q = [qii qij;
        qij qii];
    
    y0 = [populations(1)/ngroups - 1, 1, 0, 0, populations(1)/ngroups, 0, 0, 0, populations(1)/ngroups, 0, 0, 0, populations(1)/ngroups, 0, 0, 0, ...
        populations(2)/ngroups, 0, 0, 0, populations(2)/ngroups, 0, 0, 0, populations(2)/ngroups, 0, 0, 0, populations(2)/ngroups, 0, 0, 0];
    
    tspan = [0 100];
    
    sol = trudgill_network_ODE_solver(y0, q, tspan);

end


function sol = make_trudgill_network(populations, ngroups, within_group_rates, between_group_rates)

    qii = [0, datasample(within_group_rates, 1), datasample(within_group_rates, 1), datasample(within_group_rates, 1);
        datasample(within_group_rates, 1), 0, datasample(within_group_rates, 1), datasample(within_group_rates, 1);
        datasample(within_group_rates, 1), datasample(within_group_rates, 1), 0, datasample(within_group_rates, 1);
        datasample(within_group_rates, 1), datasample(within_group_rates, 1), datasample(within_group_rates, 1), 0];

    qij = [datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1);
        datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1);
        datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1);
        datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1) datasample(between_group_rates, 1)];

    q = [qii qij qij qij;
        qij qii qij qij;
        qij qij qii qij;
        qij qij qij qii];
    
    y0 = [populations(1)/ngroups - 1, 1, 0, 0, populations(1)/ngroups, 0, 0, 0, populations(1)/ngroups, 0, 0, 0, populations(1)/ngroups, 0, 0, 0, ...
        populations(2)/ngroups, 0, 0, 0, populations(2)/ngroups, 0, 0, 0, populations(2)/ngroups, 0, 0, 0, populations(2)/ngroups, 0, 0, 0, ...
        populations(3)/ngroups, 0, 0, 0, populations(3)/ngroups, 0, 0, 0, populations(3)/ngroups, 0, 0, 0, populations(3)/ngroups, 0, 0, 0, ...
        populations(4)/ngroups, 0, 0, 0, populations(4)/ngroups, 0, 0, 0, populations(4)/ngroups, 0, 0, 0, populations(4)/ngroups, 0, 0, 0];
    
    tspan = [0 100];
    
    sol = trudgill_network_ODE_solver(y0, q, tspan);

end

function plot_things(sol, populations, ngroups)
        figure
        for k = 1:length(populations) * ngroups
            subplot(length(populations),ngroups,k)
            plot(sol.x, sol.y((4*k-3):(4*k), :), 'LineWidth', 2)
            title("Society #" + num2str(ceil(k/ngroups)) + ', Community #' + num2str(k))
        end
        lgd = legend('Susceptable', 'Exposed', 'Infected', 'Recovered');
        lgd.Position(1) = 0.45;
        lgd.Position(2) = 0;
end

function t_mid = get_half_point(sol, populations, ngroups)
    % returns the timepoint where the # of susceptible people drops below 50%
    % of the populationn
    S = sol.y((1:(length(populations) * ngroups))*4 - 3, :); 

    for i = 1:(length(populations) * ngroups)
        tgt_population = populations(ceil(i/4))/ngroups/2;
        temp = S(i, :);

        if ~isempty(find((temp - tgt_population) < 0, 1))
            idx = find((temp - tgt_population) < 0, 1);
            t_mid(i) = sol.x(idx);

        else
            t_mid(i) = NaN;

        end

    end
    
    

end
