% Main code for Controller Maintenance Prioritisation,
% running a larger tank network.

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

close all; clear; clc;

%% Hyperparameter settings

% set seed for random number generators:
seed = 0;
rng(seed);

% tank network hyperparameters
T = 15; % specify number of tanks
N = 10000; % timesteps, 1 min each
controller_tuning = 5; % tunable parameter for direct synthesis tuning

% LoopRank hyperparameter:
LoopRank_factor = 0.15; % weighting factor in LoopRank

% Granger causality hyperparameters:
global GCplots kMax lagvars
alpha_gc = 0.01; % level of significance for GC connections
GCplots = 3; % set to 1 to view AIC plot and connectivity plot
lagvars = 10; % consecutive increasing values of AIC for early stopping of model order search

%% Setup for multitank model

% specify feed to system:
F_feed_SS = 5; % [m^3/h]

%% Automatically generate a tank network
% assume two banks of tanks, with underflow connections from the first bank
% to the second. The tanks in the first bank are on average twice the
% volume of the tanks in the second bank.

i_mid = ceil(T/2); % index of the last tank in the first bank

% specify tank volumes, levels and underflow parameters:
V(1:i_mid) = randi([5 7],[1,i_mid]); % [m^3], tank volumes
V(i_mid+1:T) = randi([3 5],[1,T-i_mid]); % [m^3], tank volumes
L_SP = 0.5*ones(1,T); % fraction, setpoints
L_min = 0.05; % minimum tank level, fraction
k_under = 0.5*sqrt(V); % [m^3/h], underflow parameters

% specify connections:
% (rows are source, columns are destinations)
% Feed flows is to first tank:
Feed_flows = zeros(1,T);
Feed_flows(1) = 1;

% (for Level_control_flows, all zeros in a row means this tank product
% flow, under level control, leaves the system)
% level controlled flows between each tank in the first bank:
Level_control_flows = zeros(T);
for i = 1:i_mid-1
    Level_control_flows(i,i+1) = 1;
end
% level controlled flows between some of the tanks in the second bank;
% other tanks product flows leave the system:
for i = i_mid+1:T-1
    Level_control_flows(i,i+1) = round(rand(1));
end

% (for under_flows, all zeros in a row means this tank does not have an
% underflow stream)
% underflows from each tank in the first bank to 
under_flows = zeros(T);
if i_mid == T/2
    % even number of tanks, end loop below at:
    loopend = i_mid;
else
    % odd number of loops, end loop below at:
    loopend = i_mid-1;
    % last tank in first bank also has underflow to last tank in second
    % bank:
    under_flows(loopend+1,loopend+i_mid) = 1;
end
for i = 1:loopend
    under_flows(i,i+i_mid) = 1;
end
% some of the tanks in the second bank have underflows which are recycled
% to the first tank in the first bank:
for i = i_mid+1:T-1
    under_flows(i,1) = round(rand(1));
end

% tune controllers using direct synthesis:
tc = controller_tuning; % tunable parameter

%% Some housekeeping to ensure the specifications are correct:

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

%% Run multitank network model:
[L_meas, F_in_meas, F_prod_meas, F_under_meas, F_feed_meas, Kc, tau_i, SS] = ...
    SimulateTankNetwork(V, L_SP, k_under, L_min, F_feed_SS, ...
    Feed_flows, Level_control_flows, under_flows, N, tc);

%% plot results
% figure;
% subplot(2, T+1, 1)
% plot(F_feed_meas)
% ylabel('Feed flow [m^3/hr]')
% for i = 1:T
%     subplot(2, T+1, i+1)
%     plot(L_meas(:,i))
%     ylabel(strcat({'Tank '}, int2str(i) ,' level [-]'))
%     q = gca;
%     set(q, 'YLim', [0, 1])
%     
%     subplot(2, T+1, i + T + 2)
%     hold on;
%     plot(F_in_meas(:,i))
%     plot(F_prod_meas(:,i))
%     plot(F_under_meas(:,i))
%     hold off;
%     ylabel(strcat({'Tank '}, int2str(i) ,' flows [m^3/hr]'))
% end
% subplot(2,T+1,T+2)  % merge two subplots together
% plot(0,0,  0,0,  0,0)
% axis off
% legend({'Inlet flow','Product flow','Underflow'});
% Multitank = gcf;
% Multitank.Position(3) = Multitank.Position(3)*2;
% % print(gcf,'-r200','-dtiff','Multitank')
% 
% figure;
% plot(L_meas)
% ylabel('Tank levels [-]')
% leg_cell = cell(T,1);
% for i = 1:T
%     str_cell = strcat({'Tank '}, {int2str(i)});
%     leg_cell{i} = str_cell{1};
% end
% legend(leg_cell)

%% Run LoopRank, using manually derived links
% determine links in a similar manner to that described in Bryan & Leise
% (2006): the importance of a loop is the weighted sum of the importances
% of all loops which influence that loop. The weighting is determined by
% the number of outlets a tank has.

% determine number of outputs from a given tank:
num_outputs = zeros(1,T);
for i = 1:T
    num_outputs(i) = sum(Level_control_flows(i,:) > 0) + sum(under_flows(i,:) > 0);
end

CV_link_physical = zeros(T);
for i = 1:T
    for j = 1:T
        if num_outputs(i) ~= 0
            CV_link_physical(j,i) = (Level_control_flows(i,j) + under_flows(i,j))/num_outputs(i);
        end
    end
end

ranking_physical = LoopRank(CV_link_physical, LoopRank_factor);

%% Run LoopRank, using manually derived links
% determine links in a similar manner to that described in Bryan & Leise
% (2006): the importance of a loop is the weighted sum of the importances
% of all loops which influence that loop. The weighting is determined by
% the steady state flow rates of the streams leaving each tank (but not
% leaving the system).

CV_link_flow = zeros(T);
for i = 1:T
    for j = 1:T
        denom = SS.F_prod(i) + SS.F_under(i);
        if denom ~= 0
            weight_level_control = SS.F_prod(i)/denom;
            weight_underflow = SS.F_under(i)/denom;
            CV_link_flow(j,i) = Level_control_flows(i,j)*weight_level_control + under_flows(i,j)*weight_underflow;
        end
    end
end

ranking_flow = LoopRank(CV_link_flow, LoopRank_factor);

%% Run LoopRank, using partial correlation

CV_link_pcorr = partialcorr(L_meas);

ranking_pcorr = LoopRank(CV_link_pcorr, LoopRank_factor);

%% Run LoopRank, using Granger causality

% process lag is a representative system time constant
ProcLag = mean(V)/F_feed_SS*60; % [min]; 

% maximum lag is a multiple of process lag:
kMax = round(5*ProcLag);

% note that within GC, the model order is limited to kMax, but the
% algorithm stops checking AIC after lagvars increasing values of AIC, to
% speed up large data sets

% Significance testing
alpha = alpha_gc;

dataMatrix = L_meas;
dataMatrix = zscore(dataMatrix);

gc_AllData = GrangerCausality();
[cm, pval,fStat,K,aic] = gc_AllData.grangerCaus(dataMatrix);

% refine connectivity matrix, using alpha as level of significance:
gc_AllData.refine(cm,pval,alpha);
cmRef = gc_AllData.NewCM;

names = cell(1,T);
for i = 1:T
    names{i} = int2str(i);
end

if GCplots < 3
    connMap_GC = ConnectivityMap();
    connMap_GC.drawMap(cmRef,names,'Granger Causality');
end

CV_link_GC = cmRef;

ranking_GC = LoopRank(CV_link_GC, LoopRank_factor);

%% Plot the rankings from the various methods:

figure;
hold on
plot(ranking_physical, 'o:')
plot(ranking_flow, 'd:')
plot(ranking_pcorr, '+:')
plot(ranking_GC, '*:')
hold off
ylabel('Importance score')
xlabel('Variable number')
legend('Physical', 'Flows', 'Partial correlation', 'Granger causality')

figure;
plot([0, max(ranking_flow)],[0, max(ranking_flow)],'k',ranking_flow,ranking_pcorr,'+',ranking_flow,ranking_GC,'*')
legend('Parity','Partial correlation', 'Granger causality')

RMSE_PC = sqrt(sum((ranking_flow-ranking_pcorr).^2)/T)
RMSE_GC = sqrt(sum((ranking_flow-ranking_GC).^2)/T)