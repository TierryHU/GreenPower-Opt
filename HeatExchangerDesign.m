function totalCost = heatExchangerOptimization(designVars)
% heatExchangerOptimization: Optimizes the design of a shell-and-tube heat exchanger.
%
% This function calculates the total cost of a shell-and-tube heat exchanger (STHX) design,
% including capital and operating costs, based on the specified design variables.
%
% Inputs:
%   designVars - Array of design variables:
%                [Dt, Gt, Nt, Re, L, Ntp, T_inlet, eta, i_R, TL, ec, TP]
%
% Outputs:
%   totalCost - Total cost of the heat exchanger for the given design parameters.

% Extract design variables
Dt      = designVars(1);  % Tube diameter
Gt      = designVars(2);  % Tube mass flux
Nt      = designVars(3);  % Number of tubes
Re      = designVars(4);  % Reynolds number
L       = designVars(5);  % Length of heat exchanger
Ntp     = designVars(6);  % Number of tube passes
T_inlet = designVars(7);  % Inlet temperature
eta     = designVars(8);  % Pump efficiency
i_R     = designVars(9);  % Interest rate
TL      = designVars(10); % Technical lifespan
ec      = designVars(11); % Energy unit cost
TP      = designVars(12); % Annual operating period

% Physical and material constants (example values, adjust as needed)
mu = 0.001;           % Dynamic viscosity (Pa·s)
rho = 1000;           % Density (kg/m^3)
Pr = 7;               % Prandtl number (dimensionless)
lambda = 0.6;         % Thermal conductivity (W/m·K)
Gz = Re * Pr;         % Graetz number

% Step 1: Heat Transfer Coefficient Calculation
if Re >= 2100 && Re <= 10000  % Transitional flow
    ht = transitionalHeatTransfer(Re, lambda, Dt, Pr, mu);
elseif Re > 10000              % Turbulent flow
    ht = 0.23 * lambda / Dt * Re^0.8 * Pr^(1/3) * (mu / mu)^0.14;
else                           % Laminar flow
    ht = laminarHeatTransfer(Gz, lambda, Dt, Pr, mu);
end

% Step 2: Pressure Drop Calculation
ft = 0.0035 + 0.264 / Re^0.42;  % Friction factor
deltaPt = 2 * ft * Gt^2 * L * Ntp / (Dt * rho * (mu / mu)^0.14);

% Shell-side pressure drop (simplified for illustration)
deltaPs = calculateShellPressureDrop(Ntp, rho, L, eta, Gt);

% Step 3: Capital and Operating Costs
% Capital cost
Cc = @(A) 100 * A^0.85;  % Capital cost coefficient and exponent (example)

% Operating cost
E_delta_p = (Gt * deltaPt) / (eta * rho) + (Gt * deltaPs) / (eta * rho);
operatingCost = TP * ec * E_delta_p;
capitalCost = Cc(L * Dt) * (i_R * (1 + i_R)^TL) / ((1 + i_R)^TL - 1);

% Total cost
totalCost = operatingCost + capitalCost;

end

function ht = transitionalHeatTransfer(Re, lambda, Dt, Pr, mu)
    % transitionalHeatTransfer: Computes heat transfer coefficient in transitional flow
    h_prime = 3.66 + (0.085 * Re * Pr / (1 + 0.047 * (Re * Pr)^(2/3))) * (mu / mu)^0.14;
    ht = h_prime + ((Re - 2100) / (10000 - 2100)) * (0.23 * lambda / Dt * Re^0.8 * Pr^(1/3) * (mu / mu)^0.14 - h_prime);
end

function ht = laminarHeatTransfer(Gz, lambda, Dt, Pr, mu)
    % laminarHeatTransfer: Computes heat transfer coefficient in laminar flow
    ht = (3.66 + (0.085 * Gz / (1 + 0.047 * Gz^(2/3))) * (mu / mu)^0.14) * lambda / Dt;
end

function deltaPs = calculateShellPressureDrop(Ntp, rho, L, eta, Gt)
    % calculateShellPressureDrop: Computes shell-side pressure drop (simplified)
    deltaP_c = 1.5 * Ntp * (Gt^2) / (rho * eta);
    deltaP_w = 0.8 * Ntp * (Gt^2) / (rho * eta);
    deltaP_e = 0.5 * Ntp * (Gt^2) / (rho * eta);
    deltaPs = deltaP_c + deltaP_w + deltaP_e;
end
