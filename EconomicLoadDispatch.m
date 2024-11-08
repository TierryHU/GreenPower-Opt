function f = funcELD(x, D, costCoefficients, powerLimits)
% funcELD: Economic Load Dispatch (ELD) Optimization
%
% This function calculates the fuel cost for an economic load dispatch problem.
% The objective is to minimize generation costs while satisfying demand and constraints.
%
% Inputs:
%   x - Array of generator outputs [P1, P2, ..., Pn] (in MW)
%   D - Total demand to be satisfied (in MW)
%   costCoefficients - Matrix of cost coefficients for each generator, where each row is [alpha, beta, gamma]
%   powerLimits - Matrix of power limits for each generator, where each row is [Pmin, Pmax]
%
% Outputs:
%   f - Total fuel cost for given generator outputs
%
% Example usage:
%   Define the generator outputs, demand, cost coefficients, and power limits as follows:
%     x = [100, 120, 150];                  % Example output power for each generator in MW
%     D = 370;                               % Example demand in MW
%     costCoefficients = [20 0.04 0.0001;    % Cost coefficients for each generator
%                         15 0.03 0.0002;
%                         30 0.05 0.00015];
%     powerLimits = [50 200;                % Power limits for each generator (Pmin, Pmax)
%                    30 150;
%                    40 180];
%   Call the function:
%     total_cost = funcELD(x, D, costCoefficients, powerLimits);

% Initialize total cost and define penalty for constraint violation
numGenerators = length(x);
totalCost = 0;
penalty = 1e6;  % Penalty for violating constraints

% Calculate the generation cost for each generator
for i = 1:numGenerators
    P = x(i);  % Power output of the i-th generator in MW
    
    % Enforce generator power limits
    if P < powerLimits(i, 1) || P > powerLimits(i, 2)
        f = penalty;  % Return a high cost if generator output is out of bounds
        return;
    end
    
    % Extract cost coefficients for this generator
    alpha = costCoefficients(i, 1);
    beta = costCoefficients(i, 2);
    gamma = costCoefficients(i, 3);
    
    % Calculate fuel cost using the quadratic cost function
    totalCost = totalCost + (alpha + beta * P + gamma * P^2);
end

% Enforce power balance constraint (total generation should match demand)
if abs(sum(x) - D) > 1e-3
    f = penalty;  % Return a high cost if total generation does not meet demand
else
    f = totalCost;  % Return the calculated total generation cost
end

end

% Example usage:
% Define generator outputs, demand, cost coefficients, and power limits
x = [100, 120, 150];  % Power output for each generator in MW
D = 370;              % Total demand in MW
costCoefficients = [20 0.04 0.0001; 15 0.03 0.0002; 30 0.05 0.00015];
powerLimits = [50 200; 30 150; 40 180];

% Calculate total cost for given outputs
total_cost = funcELD(x, D, costCoefficients, powerLimits);
fprintf('Total cost for given generator outputs: %.2f\n', total_cost);
