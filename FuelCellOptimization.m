function mse = pemFuelCellOptimization(params, experimentalData)
% pemFuelCellOptimization: Identifies optimal parameters for PEM fuel cell model
%
% This function calculates the Mean Square Error (MSE) between simulated
% and experimental output voltages for a PEM fuel cell model.
%
% Inputs:
%   params - Array of model parameters [A, l, R_C, xi1, xi2, xi3, xi4, J_n, J_min, lambda, b]
%   experimentalData - Struct containing fields for experimental voltage (Uq),
%                      current density (J), and temperature (T)
%
% Outputs:
%   mse - Mean Square Error between experimental and simulated voltage

% Unpack experimental data
Uq = experimentalData.Uq; % Experimental voltage data
J = experimentalData.J;   % Current density data
T = experimentalData.T;   % Temperature data
Q = length(Uq);           % Number of data points

% Define constants
F = 96485;            % Faraday constant (C/mol)
R = 8.314;            % Universal gas constant (J/(mol·K))
PH2 = 1.5;            % Partial pressure of hydrogen (atm)
PO2 = 0.5;            % Partial pressure of oxygen (atm)
n = 1;                % Number of cells in series
T_ref = 298.15;       % Reference temperature (K)

% Extract parameters
A = params(1);
l = params(2);
R_C = params(3);
xi1 = params(4);
xi2 = params(5);
xi3 = params(6);
xi4 = params(7);
J_n = params(8);
J_min = params(9);
lambda = params(10);
b = params(11);

% Initialize error accumulator
errorSum = 0;

% Loop through each data point to calculate simulated voltage and error
for q = 1:Q
    % Compute variables and terms based on equations provided
    
    % Calculate E_Nernst
    deltaG = -237.13;  % Change in Gibbs free energy (J/mol) at standard state
    deltaS = -163.5;   % Change in entropy (J/(mol·K))
    E_nernst = (deltaG / (2 * F)) + (deltaS / (2 * F) * (T(q) - T_ref)) + (R * T(q) / (2 * F) * log(PH2 * sqrt(PO2)));

    % Calculate membrane resistivity (rho_M)
    rho_M = (181.6 * (1 + 0.03 * (J(q) / A) + 0.062 * (T(q) / 303)^2 * (J(q) / A)^2.5)) / ...
        ((lambda - 0.634 - 3 * (J(q) / A)) * exp(4.18 * (T(q) - 303) / T(q)));
    
    % Calculate R_M (membrane resistance)
    R_M = rho_M * l / A;

    % Oxygen concentration at the cathode
    CO2 = PO2 / (5.08 * 10^6 * exp(-498 / T(q)));

    % Activation, ohmic, and concentration voltage drops
    V_act = xi1 * T(q) + xi2 * T(q) * log(CO2) + xi3 * T(q) * log(J(q));
    V_ohmic = J(q) * (R_M + R_C);
    V_con = b * log(1 - J(q) / J_n);

    % Simulated voltage
    Vs = n * (E_nernst - V_act - V_ohmic - V_con);

    % Calculate the squared error
    errorSum = errorSum + (Uq(q) - Vs)^2;
end

% Compute MSE
mse = errorSum / Q;

end
