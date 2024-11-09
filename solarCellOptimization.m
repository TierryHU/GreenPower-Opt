function mse = solarCellOptimization(params, experimentalData)
% solarCellOptimization: Identifies optimal parameters for Solar Cell Model
%
% This function calculates the Mean Square Error (MSE) for identifying the
% parameters of a photovoltaic (PV) cell based on the objective function.
%
% Inputs:
%   params - Array of model parameters [a, Rs, Rp, T]
%   experimentalData - Struct containing experimental values for Isc, Voc, and Pmp
%
% Outputs:
%   mse - Objective function value representing the accumulated error

% Extract experimental data
Isc_exp = experimentalData.Isc;  % Short-circuit current (experimental)
Voc_exp = experimentalData.Voc;  % Open-circuit voltage (experimental)
Pmp_exp = experimentalData.Pmp;  % Maximum power point (experimental)

% Define constants
q = 1.602e-19;      % Charge of an electron (Coulombs)
k = 1.38e-23;       % Boltzmann constant (J/K)

% Extract parameters from input
a = params(1);    % Ideality factor
Rs = params(2);   % Series resistance
Rp = params(3);   % Parallel resistance
T = params(4);    % Temperature (in Kelvin)

% Calculate thermal voltage (Vt)
Vt = (k * T) / q;

% Calculate simulated values
Isc_sim = Isc(params, T);     % Simulated short-circuit current
Voc_sim = Voc(params, Vt);    % Simulated open-circuit voltage
Pmp_sim = Pmp(params, Vt);    % Simulated maximum power point

% Calculate error terms
f_Isc = (Isc_sim - Isc_exp) / Isc_exp;
f_Voc = (Voc_sim - Voc_exp) / Voc_exp;
f_Pmp = (Pmp_sim - Pmp_exp) / Pmp_exp;

% Objective function: Sum of absolute errors
mse = abs(f_Isc) + abs(f_Voc) + abs(f_Pmp);

end

% Additional helper functions for calculating Isc, Voc, and Pmp

function Isc_sim = Isc(params, T)
% Simulates the short-circuit current
% Inputs:
%   params - Array [a, Rs, Rp, T]
%   T - Operating temperature (K)
% Output:
%   Isc_sim - Simulated short-circuit current

    % Constants
    k = 1.38e-23;  % Boltzmann constant (J/K)
    q = 1.602e-19; % Electron charge (C)
    a = params(1);
    
    % Calculate short-circuit current (as a function of temperature and ideality factor)
    Isc_sim = T * exp(-q / (a * k * T));  % Simplified for illustration
end

function Voc_sim = Voc(params, Vt)
% Simulates the open-circuit voltage
% Inputs:
%   params - Array [a, Rs, Rp, T]
%   Vt - Thermal voltage (V)
% Output:
%   Voc_sim - Simulated open-circuit voltage

    a = params(1);
    Isc = Isc(params, Vt);
    Voc_sim = Vt * log((Isc + 1) / (1 - Isc));  % Simplified calculation
end

function Pmp_sim = Pmp(params, Vt)
% Simulates the maximum power point
% Inputs:
%   params - Array [a, Rs, Rp, T]
%   Vt - Thermal voltage (V)
% Output:
%   Pmp_sim - Simulated maximum power point

    Isc = Isc(params, Vt);
    Voc = Voc(params, Vt);
    Pmp_sim = Isc * Voc * 0.9; % Estimated at 90% of Voc and Isc for illustration
end
