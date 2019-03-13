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

function [L_meas, F_in_meas, F_prod_meas, F_under_meas, F_feed_meas, Kc, tau_i, SS] = ...
    SimulateTankNetwork(V, L_SP, k_under, L_min, F_feed_SS, ...
    Feed_flows, Level_control_flows, under_flows, N, tc)

% Function to simulate a network of tanks under level control, with
% disturbances in the feed to the network. Measurement noise is added to
% true values, as normally distributed noise with a standard deviation of
% 0.05 by default.

%% Inputs specifying flows, volumes, etc:
% V, [m^3], tank volumes
% L_SP, fraction, tank level setpoints
% k_under, [m^3/h], underflow parameters
% L_min, fraction, minimum tank level
% F_feed_SS, [m^3/h], steady state feed to network

%% Inputs specifying tank connectivity:
% Rows are source, columns are destinations
% Each row must add up to 1, allowing splits in the flows to multiple tanks

% Feed_flows, specifying the tank to which feed stream connects

% Level_control_flows, specifying the connections under level control; all
% zeros in a row means the tank product flow, under level control, leaves
% the system.

% under_flows, specifying the connections of the underflows; all zeros in a
% row means this tank does not have an underflow stream.

%% Other inputs:
% N, timesteps, 1 min each
% tc, tunable parameter for direct synthesis controller tuning

%% Outputs:
% L_meas, fraction, measured tank levels
% F_in_meas, [m^3/h], measured flow to each tank
% F_prod_meas, [m^3/h], measured product flows from each tank
% F_under_meas, [m^3/h], measured underflows from each tank
% F_feed_meas, [m^3/h], measured flow to network
% Kc, controller proportional gains
% tau_i, controller integral gains
% SS, structure with the steady state input, product and output flows

%% begin function
T = length(V); % number of tanks

% ensure each row of Feed_flows, Level_control_flows, and under_flows sums
% to 1 (mass balance requirement):
Feed_flows = Feed_flows/sum(Feed_flows);
for i = 1:T
    if (sum(Level_control_flows(i,:)) ~= 0)
        Level_control_flows(i,:) = Level_control_flows(i,:)/sum(Level_control_flows(i,:));
    end
    
    if (sum(under_flows(i,:)) ~= 0)
        under_flows(i,:) = under_flows(i,:)/sum(under_flows(i,:));
    end
end

% update k_under for tanks with no underflow:
for i = 1:T
    k_under(i) = sum(under_flows(i,:))*k_under(i);
end

% calculate steady state underflows:
F_under_SS = zeros(1, T);
for i = 1:T
    F_under_SS(i) = sum(under_flows(i,:))*k_under(i)*L_SP(i);
end
SS.F_under = F_under_SS;

% calculate steady state feed and product flows:
F_in_SS = zeros(1, T);
F_prod_SS = zeros(1, T);

for i = 1:T
    F_in_SS(i) = F_feed_SS*Feed_flows(i) + F_prod_SS*Level_control_flows(:,i) + F_under_SS*under_flows(:,i);
    F_prod_SS(i) = F_in_SS(i) - F_under_SS(i);
end
SS.F_in = F_in_SS;
SS.F_prod = F_prod_SS;

% tune controllers using direct synthesis:
tau = V./F_in_SS/60; % process residence time based on feed flows
K = -V./F_prod_SS/60; % process gains based on product flows

Kc = 1./K.*tau/tc; % controller proportional gains
tau_i = tau; % controller integral time constant

%% Run tank models with disturbance in feed:

% define matrices to store actual values of variables:
F_feed = zeros(N,1); % [m^3/h]
F_in = zeros(N,T); % [m^3/h]
F_prod = zeros(N,T); % [m^3/h]
F_under = zeros(N,T); % [m^3/h]
L = zeros(N,T); % fraction

F_prod_spec = zeros(1,T); % [m^3/h]

% define matrices to store measured values of variables:
F_feed_meas = zeros(N,1); % [m^3/h]
F_in_meas = zeros(N,T); % [m^3/h]
F_prod_meas = zeros(N,T); % [m^3/h]
F_under_meas = zeros(N,T); % [m^3/h]
L_meas = zeros(N,T); % fraction

% initialise at steady state:
F_feed(1) = F_feed_SS;
F_in(1,:) = F_in_SS;
F_prod(1,:) = F_prod_SS;
F_under(1,:) = F_under_SS;
L(1,:) = L_SP;

% add noise:
F_feed_meas(1) = AddNoise(F_feed_SS);
F_in_meas(1,:) = AddNoise(F_in_SS);
F_prod_meas(1,:) = AddNoise(F_prod_SS);
F_under_meas(1,:) = AddNoise(F_under_SS);
L_meas(1,:) = AddNoise(L_SP);

for i = 2:N
    F_feed(i) = 0.97*F_feed(i-1) + normrnd(0,0.15) + F_feed_SS*(1-0.97);
    F_feed_meas(i) = AddNoise(F_feed(i));
    
    % calculate controller outputs:
    for j = 1:T
        F_prod_spec(j) = TankController(F_prod_meas(i-1,j), L_meas(i-1,j), L_SP(j), Kc(j), tau_i(j));
    end
    
    % run tank model for each tank:
    for j = 1:T
        F_in(i,j) = F_feed(i-1)*Feed_flows(j) + F_prod(i-1,:)*Level_control_flows(:,j) + F_under(i-1,:)*under_flows(:,j);
        F_in_meas(i,j) = AddNoise(F_in(i,j));
        
        [L(i,j), F_prod(i,j), F_under(i,j)] = ...
            TankModel(L(i-1,j), F_in(i,j), V(j), k_under(j), F_prod_spec(j), L_min);
        L_meas(i,j) = AddNoise(L(i,j));
        F_prod_meas(i,j) = AddNoise(F_prod(i,j));
        F_under_meas(i,j) = AddNoise(F_under(i,j));
    end
end