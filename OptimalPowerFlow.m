function total_cost = funcOPF(Pgen, Qgen, V, theta, demand, G, B, costCoeffs, limits)
% funcOPF: Optimal Power Flow (OPF) Optimization Function
%
% This function calculates the total fuel cost for an OPF problem by determining
% the ideal active and reactive power outputs of each generator and distribution
% across the network.
%
% Inputs:
%   Pgen - Array of active power outputs of each generator (MW)
%   Qgen - Array of reactive power outputs of each generator (MVAR)
%   V - Array of voltage magnitudes at each bus (pu)
%   theta - Array of voltage angle differences between buses (radians)
%   demand - Struct containing active and reactive demand for each bus (MW and MVAR)
%   G - Conductance matrix between buses
%   B - Susceptance matrix between buses
%   costCoeffs - Matrix of cost coefficients for each generator [a, b, c]
%   limits - Struct containing the following fields:
%       - PgenMin, PgenMax: Active power limits for generators (MW)
%       - QgenMin, QgenMax: Reactive power limits for generators (MVAR)
%       - Vmin, Vmax: Voltage limits for each bus (pu)
%       - Smax: Apparent power limits for each line (MVA)
%
% Outputs:
%   total_cost - Total economic cost of the optimal power flow solution
%
% Example usage:
%   % Define system parameters for a 3-bus, 3-generator system
%   Pgen = [80; 100; 120];               % Active power outputs (MW)
%   Qgen = [30; 40; 20];                 % Reactive power outputs (MVAR)
%   V = [1.05; 1.01; 1.0];               % Voltage magnitudes at buses (pu)
%   theta = [0; 0.05; -0.03];            % Voltage angle differences (radians)
%   demand.P = [70; 60; 50];             % Active power demand (MW)
%   demand.Q = [25; 20; 15];             % Reactive power demand (MVAR)
%   G = [0.1, -0.02, 0; -0.02, 0.1, -0.03; 0, -0.03, 0.05]; % Conductance matrix
%   B = [0.15, -0.05, 0; -0.05, 0.2, -0.1; 0, -0.1, 0.1];   % Susceptance matrix
%   costCoeffs = [30, 0.2, 0.01; 28, 0.18, 0.015; 35, 0.25, 0.02]; % Cost coefficients
%   limits.PgenMin = [50; 80; 90];
%   limits.PgenMax = [150; 200; 180];
%   limits.QgenMin = [10; 20; 15];
%   limits.QgenMax = [40; 50; 30];
%   limits.Vmin = [0.95; 0.98; 0.97];
%   limits.Vmax = [1.05; 1.02; 1.01];
%   limits.Smax = [100; 80; 90];  % Apparent power limits (MVA)
%   % Calculate total cost
%   total_cost = funcOPF(Pgen, Qgen, V, theta, demand, G, B, costCoeffs, limits);

% Initialization
numBuses = length(V);
numGens = length(Pgen);
total_cost = 0;   % Initialize total cost
penalty = 1e6;    % Penalty for constraint violations

% Calculate active and reactive power balance
for i = 1:numBuses
    % Active power balance (Equation limit9)
    P_balance = Pgen(i) - V(i) * sum(V .* (G(i, :) .* cos(theta(i) - theta) + B(i, :) .* sin(theta(i) - theta)));
    if abs(P_balance - demand.P(i)) > 1e-3
        total_cost = penalty; % Apply penalty if balance is not met
        return;
    end

    % Reactive power balance (Equation limit10)
    Q_balance = Qgen(i) - V(i) * sum(V .* (G(i, :) .* sin(theta(i) - theta) + B(i, :) .* cos(theta(i) - theta)));
    if abs(Q_balance - demand.Q(i)) > 1e-3
        total_cost = penalty; % Apply penalty if balance is not met
        return;
    end
end

% Calculate total cost from generator outputs
for i = 1:numGens
    % Check active and reactive power limits
    if Pgen(i) < limits.PgenMin(i) || Pgen(i) > limits.PgenMax(i) || ...
       Qgen(i) < limits.QgenMin(i) || Qgen(i) > limits.QgenMax(i)
        total_cost = penalty; % Apply penalty if power output is out of limits
        return;
    end
    
    % Fuel cost calculation (Equation obj3)
    a = costCoeffs(i, 1);
    b = costCoeffs(i, 2);
    c = costCoeffs(i, 3);
    total_cost = total_cost + (a + b * Pgen(i) + c * Pgen(i)^2);
end

% Check voltage magnitude constraints for each bus
if any(V < limits.Vmin) || any(V > limits.Vmax)
    total_cost = penalty; % Apply penalty if voltage magnitude is out of limits
    return;
end

% Apparent power limit for each line (Equation limit14)
for k = 1:numBuses
    S_k = V(k) * abs(G(k, :) * V .* cos(theta(k) - theta) + B(k, :) * V .* sin(theta(k) - theta));
    if S_k > limits.Smax(k)
        total_cost = penalty; % Apply penalty if apparent power exceeds limit
        return;
    end
end

end
