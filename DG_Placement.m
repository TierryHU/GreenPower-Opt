function [bestLoss, bestSite, bestSize] = optimalDGPlacement(nBuses, demand, G, B, V_limits, Smax_limits, fobj)
% optimalDGPlacement: Function to find the optimal site and size for DG placement.
%
% This function optimizes the placement and sizing of distributed generation (DG)
% units to minimize power losses within a distribution network.
%
% Inputs:
%   nBuses      - Number of buses in the system
%   demand      - Struct with fields P and Q for each bus representing power demand
%   G, B        - Conductance and susceptance matrices for power network
%   V_limits    - Voltage limits for each bus [Vmin, Vmax] (pu)
%   Smax_limits - Maximum apparent power limit for each line (MVA)
%   fobj        - Function handle for the objective function
%
% Outputs:
%   bestLoss    - Minimum power loss obtained
%   bestSite    - Optimal bus location for DG placement
%   bestSize    - Optimal size (MW) of the DG unit at bestSite
%
% Example usage:
%   nBuses = 5;
%   demand.P = [0; 10; 15; 20; 10];       % Power demand in MW
%   demand.Q = [0; 5; 7; 10; 5];          % Reactive demand in MVAR
%   G = [0.1, -0.01, -0.01, 0, 0;         % Conductance matrix
%        -0.01, 0.1, -0.01, -0.01, 0;
%        -0.01, -0.01, 0.1, -0.01, -0.01;
%        0, -0.01, -0.01, 0.1, -0.01;
%        0, 0, -0.01, -0.01, 0.1];
%   B = [0.05, -0.005, -0.005, 0, 0;      % Susceptance matrix
%        -0.005, 0.05, -0.005, -0.005, 0;
%        -0.005, -0.005, 0.05, -0.005, -0.005;
%        0, -0.005, -0.005, 0.05, -0.005;
%        0, 0, -0.005, -0.005, 0.05];
%   V_limits = [0.95, 1.05];              % Voltage limits
%   Smax_limits = [100; 80; 90; 75; 80];  % Apparent power limits (MVA)
%   fobj = @(DG) calcPowerLoss(DG, demand, G, B, V_limits, Smax_limits);
%   [bestLoss, bestSite, bestSize] = optimalDGPlacement(nBuses, demand, G, B, V_limits, Smax_limits, fobj);

% Initialize search variables
minLoss = Inf;
bestSite = 0;
bestSize = 0;

% Loop through potential DG placements
for i = 2:nBuses  % Assuming bus 1 is slack and cannot host DG
    for size = 0.01:0.01:0.63  % DG size ranges from 0 to 0.63 MW
        % Evaluate power loss for this site and size
        DG = struct('site', i, 'size', size);
        Ploss = fobj(DG);
        
        % Update if a better configuration is found
        if Ploss < minLoss
            minLoss = Ploss;
            bestSite = i;
            bestSize = size;
        end
    end
end

% Return best results
bestLoss = minLoss;

end

function Ploss = calcPowerLoss(DG, demand, G, B, V_limits, Smax_limits)
% calcPowerLoss: Calculates power loss given a DG configuration.
%
% This helper function computes the total system power loss for a given DG placement
% and sizing, based on network conductance and susceptance matrices.
%
% Inputs:
%   DG         - Struct with fields site and size representing DG placement and size
%   demand     - Struct with fields P and Q for each bus
%   G, B       - Conductance and susceptance matrices
%   V_limits   - Voltage limits [Vmin, Vmax] (pu)
%   Smax_limits- Apparent power limits for each line (MVA)
%
% Outputs:
%   Ploss      - Calculated power loss for the given DG placement
%
% Example usage in optimalDGPlacement function

% Calculate active and reactive power injection with DG at site i
P = demand.P;
Q = demand.Q;
P(DG.site) = P(DG.site) - DG.size;  % Update demand with DG generation

% Initialize loss
Ploss = 0;

% Calculate line losses based on conductance matrix
for i = 1:length(P)
    for j = i+1:length(P)
        if G(i, j) ~= 0
            % Calculate power on line i-j
            P_line = P(i) - P(j);
            Ploss = Ploss + P_line^2 * G(i, j);  % Add to total loss
        end
    end
end

% Check voltage and apparent power constraints (optional)
% This code assumes penalty-based constraint handling is implemented in a higher-level function

end
