function total_cost = funcUC(P, U, demand, spinning_reserve, costCoefficients, powerLimits, minUpDownTime, startUpCost)
% funcUC: Unit Commitment (UC) Optimization Function
%
% This function calculates the total economic cost for a Unit Commitment (UC) problem.
% It minimizes fuel and startup costs for a given set of generators and time slots.
%
% Inputs:
%   P - Matrix of generator outputs for each time slot (size: [n, T] where n is the number of generators and T is the number of time slots)
%   U - Binary matrix indicating generator on/off status for each time slot (size: [n, T])
%   demand - Array of demand values for each time slot (size: [1, T])
%   spinning_reserve - Array of spinning reserve requirements for each time slot (size: [1, T])
%   costCoefficients - Matrix of cost coefficients for each generator [alpha, beta, gamma] (size: [n, 3])
%   powerLimits - Matrix of power limits [Pmin, Pmax] for each generator (size: [n, 2])
%   minUpDownTime - Matrix of minimum up/down times for each generator [MUT, MDT] (size: [n, 2])
%   startUpCost - Matrix of start-up costs [Hot, Cold] for each generator (size: [n, 2])
%
% Outputs:
%   total_cost - Total economic cost of the unit commitment solution
%
% Example usage:
%   % Define inputs for a 3-generator, 24-hour system
%   n = 3; T = 24;
%   P = randi([0, 100], n, T);                   % Random generator outputs (MW)
%   U = randi([0, 1], n, T);                     % Random on/off status (0 or 1)
%   demand = 200 + 50 * sin(1:T);                % Example demand profile (MW)
%   spinning_reserve = 0.1 * demand;             % Example spinning reserve
%   costCoefficients = [20, 0.04, 0.0001; 15, 0.03, 0.0002; 30, 0.05, 0.00015]; % Cost coefficients
%   powerLimits = [50, 200; 30, 150; 40, 180];   % Power limits (MW)
%   minUpDownTime = [3, 2; 4, 3; 5, 3];         % Minimum up/down time (hours)
%   startUpCost = [500, 1000; 300, 600; 400, 900]; % Startup costs [Hot, Cold] (currency)
%   % Calculate total cost
%   total_cost = funcUC(P, U, demand, spinning_reserve, costCoefficients, powerLimits, minUpDownTime, startUpCost);

% Initialization
[n, T] = size(P); % Number of generators and time slots
total_cost = 0;   % Initialize total cost
penalty = 1e6;    % Penalty for constraint violations

% Loop through each time slot
for t = 1:T
    % Power balance and spinning reserve constraint
    total_power = sum(P(:, t) .* U(:, t)); % Sum of online generators' output
    if total_power < demand(t) + spinning_reserve(t)
        total_cost = penalty; % Apply penalty if demand + reserve is not met
        return;
    end
    
    % Loop through each generator
    for j = 1:n
        % Check generator power limits
        if P(j, t) < powerLimits(j, 1) * U(j, t) || P(j, t) > powerLimits(j, 2) * U(j, t)
            total_cost = penalty; % Apply penalty if output is outside limits
            return;
        end
        
        % Calculate fuel cost
        alpha = costCoefficients(j, 1);
        beta = costCoefficients(j, 2);
        gamma = costCoefficients(j, 3);
        total_cost = total_cost + (alpha + beta * P(j, t) + gamma * P(j, t)^2) * U(j, t);
        
        % Calculate startup cost
        if U(j, t) == 1 && (t == 1 || U(j, t-1) == 0) % Generator is starting up
            if t == 1 || U(j, t-1) == 0 % If previous state was off
                % Determine hot or cold start cost based on previous downtime
                downtime = sum(U(j, max(1, t - minUpDownTime(j, 2)):t-1) == 0);
                if downtime >= minUpDownTime(j, 2)
                    total_cost = total_cost + startUpCost(j, 2); % Cold start
                else
                    total_cost = total_cost + startUpCost(j, 1); % Hot start
                end
            end
        end
    end
end

% Check minimum up/down time constraints
for j = 1:n
    % Minimum up time
    on_times = find(diff([0, U(j, :)]) == 1); % Indices where generator turns on
    for k = 1:length(on_times)
        if k < length(on_times)
            next_off = on_times(k) + minUpDownTime(j, 1);
            if next_off <= T && sum(U(j, on_times(k):next_off-1)) < minUpDownTime(j, 1)
                total_cost = penalty; % Apply penalty if up time is too short
                return;
            end
        end
    end
    % Minimum down time
    off_times = find(diff([U(j, :), 0]) == -1); % Indices where generator turns off
    for k = 1:length(off_times)
        if k < length(off_times)
            next_on = off_times(k) + minUpDownTime(j, 2);
            if next_on <= T && sum(U(j, off_times(k):next_on-1) == 0) < minUpDownTime(j, 2)
                total_cost = penalty; % Apply penalty if down time is too short
                return;
            end
        end
    end
end

end
