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

function ranking = LoopRank(CV_link_matrix, m)
% Function to return the ranking of a set of control loops, with the
% connectivity/link matrix of the controlled values in each loop given in
% CV_link_matrix. The value of m, [0, 1], determines the weighting for the
% weighted average described in Bryan & Leise (2006); for m in (0, 1], the
% method is guaranteed to return a unique ranking. For m = 0, if there are
% sub-networks within CV_link_matrix, the method may return ambiguous
% rankings. For m = 0, the method is the same as described in Farenzana &
% Trierweiler (2009).

%% Inputs:
% CV_link_matrix, an N by N matrix in which element ij indicates the
% strength of connection from variable i to variable j. CV_link_matrix can
% be determined manually, from correlation methods, or from causal methods.
% N is the number of variables.
% m, weighting factor. default 0.15 as in Bryan & Leise (2006).
if nargin < 2
    m = 0.15;
end

%% Outputs:
% ranking, a column vector with N elements. The value in element i reflects
% the relative importance of the i-th control loop. Values are scaled to
% sum to 1; larger values indicate more important loops.

%% begin function
% preprocess CV_link_matrix to ensure it meets the requirements of the
% PageRank/LoopRank method:

[N, ~] = size(CV_link_matrix);

% take absolute value, as strength of connection is important:
CV_link_matrix = abs(CV_link_matrix);

% set diagonal elements to 0, as a loop is not connected to itself:
for i = 1:N
    CV_link_matrix(i,i) = 0;
end

% normalise each column, to provide correct weightings:
for i = 1:N
    if sum(CV_link_matrix(i,:)) ~= 0
        CV_link_matrix(:,i) = CV_link_matrix(i,:)/sum(CV_link_matrix(i,:));
    end
end

% use the weighted average trick from Bryan and Leise (2006):
S = ones(N)/N; % matrix with each entry 1/N
CV_link_matrix = (1-m)*CV_link_matrix + m*S;

% find eigenvalues and eigenvectors of CV_link_matrix:
[V,D] = eig(CV_link_matrix,'vector');

% find the eigenvalue closest to 1:
[~, eig_index] = min((D - 1).^2);

% find corresponding eigenvector:
eig_vec = V(:,eig_index);

% scale ranking and return a column vector:
ranking = eig_vec/sum(eig_vec);
ranking = ranking(:);