% Main code for Controller Maintenance Prioritisation,
% running multiple replicates of larger tank network.

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

close all; clear all; clc;

% specify the vector of network sizes:
Tvec = 5:15;
% specify the vector of random number generator seeds:
seedvec = 0:15;

GC_rmse = zeros(length(seedvec),length(Tvec));
PC_rmse = GC_rmse;

for j = length(Tvec)
    j
    for i = 1:length(seedvec)
        i
        Tin = Tvec(j);
        seedin = seedvec(i);
        [RMSE_PC,RMSE_GC] = jtmgrid(Tin,seedin);
        GC_rmse(i,j) = RMSE_GC;
        PC_rmse(i,j) = RMSE_PC;
    end
end

GC_mean = mean(GC_rmse);
PC_mean = mean(PC_rmse);

GC_std = std(GC_rmse);
PC_std = std(PC_rmse);

figure;
hold on
errorbar(Tvec,PC_mean,PC_std, 'o:')
errorbar(Tvec,GC_mean,GC_std, 'x:')
hold off
xlabel('Number of tanks in network')
ylabel('Root mean square error')
legend('Partial correlation','Granger causality')
% print(gcf,'-r250','-dtiff','PCvsGC')