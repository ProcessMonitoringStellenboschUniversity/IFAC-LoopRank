% Controller Maintenance Prioritisation

% Copyright (C) 2019 JT McCoy and L Auret

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.

% Controller maintenance prioritisation with LoopRank for a milling and flotation plant
% submitted to IFAC MMM 2019 Conference

% JT McCoy and L Auret
% Department of Process Engineering
% University of Stellenbosch

function [L_next, F_prod, F_under] = ...
    TankModel(L, F_in, V, k_under, F_prod_spec, L_min)

% Model of a tank with a flow in (F_in), a product flow out (F_prod) which
% is under level control (F_prod_spec is controller setting of flow rate),
% underflow determed by gravity (F_under = k_under*L), volume (V), and level
% fraction L, between 0 and 1, with a minimum level L_min.

% Note time step assumed to be 1 min.

% Inputs:
    % L - current level, fraction
    % F_in - current flow of liquid into tank, [m^3/h]
    % V - tank volume, [m^3]
    % k_under - flow coefficient for underflow, F = k*L, [m^3/h]
    % F_prod_spec - product outflow, specified by controller, [m^3/h]
    % L_min - minimum tank level, fraction

% Outputs:
    % L_next - resulting level, fraction
    % F_prod - resulting product flow, [m^3/h]
    % F_under - resulting underflow, [m^3/h]

%% Calculate flow out due to gravity
F_under_temp = k_under*L;

%% Calculate next level
% time step is 1 min, so divide by 60 to convert to hours:
L_next_temp = L + (1/V)*(F_in - F_under_temp - F_prod_spec)/60;

%% Check for overflow
if L_next_temp > 1
    % increase F_prod_spec to maintain level at 1:
    F_prod = F_prod_spec + V*(L_next_temp - 1)*60;
    
    F_under = F_under_temp;
    
    % specify L_next:
    L_next = 1;

    %% Check for empty tank
elseif L_next_temp < L_min
    % gravity flow is uncontrolled, check if sufficient level for underflow
    % as calculated in F_under_temp, setting F_prod = 0:
    L_next_empty = L + (1/V)*(F_in - F_under_temp)/60;
    
    if L_next_empty < L_min
        % insufficient level for calculated gravity flow in F_under_temp
        % set underflow so L_next = L_min
        F_under = F_in + F_under_temp - V*(0 - L_next_empty)*60;
        F_prod = 0;
        L_next = L_min;
    else
        % sufficient level for calculated gravity flow in F_under_temp
        F_under = F_under_temp;
        % set F_prod so L_next = L_min
        F_prod = V*(L_next_empty - L_min)*60;
        L_next = L_min;
    end
else
    L_next = L_next_temp;
    F_prod = F_prod_spec;
    F_under = F_under_temp;
end

%% Ensure flow rates positive:
F_prod = max(0, F_prod);
F_under = max(0, F_under);