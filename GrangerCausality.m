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

% Acknowledgements to Brian Lindner for supplying Granger causality code.
% For open source version, see: 
% https://github.com/ProcessMonitoringStellenboschUniversity/IFAC-Causality-Analysis

classdef GrangerCausality<handle
    %GRANGERCAUSALITY class for calculation of Granger Causality
    
    properties
        
        ConnMatrix
        NewCM
        PVals
        Parameters        
        
    end
    
    methods
        
        function this = GrangerCausality()
            %Granger Causality Class constructor
            
        end
        
        function this = calculate(this,data)
            
            [cm, pval,~,k,~] = this.grangerCaus(data);
            this.refine(cm,pval,0.05)
            % Package results
            this.ConnMatrix = cm;
            this.PVals = pval;
%             this.Parameters.KMax = kMax;
            this.Parameters.K = k;
            
        end
        
        function [cm, pval,fStat,K,aic] = grangerCaus(this,X,varargin)
            %GRANGERCAUS Granger causality 
            global GCplots kMax lagvars
            [N,M] = size(X);
            
            if nargin>2
                % if model order is specified
                K = varargin{1};
                aic = [];
            else
                % Model order estimation
                aic = zeros(kMax,1) + Inf;
                SIGMAfull = cell(kMax,1);
                E = cell(kMax,1);
                
%                 for k = 1:kMax
%                     % Get residuals covariance matrix
%                     [~,SIGMAfull{k},E{k}] = this.dataToVar(X,k);
%                     % Using the AIC from Duan (2014) (Eq 4.18)
% %                     fprintf('Model order = %d\n',k);
%                     % working out aic with just sigma from X to Y
%                     aic(k) = log(det(SIGMAfull{k})) + (2*k*M^2)/N;
% %                     aic(k) = (N-k)*log(det(SIGMAfull{k}(1))) + (N-k)*(2*k*M^2)/((N-k)-(k*M^2)-1);%Toolbox uses a corrected version
%                 end
                
                % generate first few data points of aic:
                for k = 1:lagvars+1
                    % Get residuals covariance matrix
                    [~,SIGMAfull{k},E{k}] = this.dataToVar(X,k);
                    % Using the AIC from Duan (2014) (Eq 4.18)
%                     fprintf('Model order = %d\n',k);
                    % working out aic with just sigma from X to Y
                    aic(k) = log(det(SIGMAfull{k})) + (2*k*M^2)/N;
%                     aic(k) = (N-k)*log(det(SIGMAfull{k}(1))) + (N-k)*(2*k*M^2)/((N-k)-(k*M^2)-1);%Toolbox uses a corrected version
                end
                
                % continue generating aic until upward trend occurs
                for k = lagvars+2:kMax
%                     display(k/kMax)
                    gtaic = aic(k-lagvars:k-1) > aic(k-lagvars-1:k-2); % vector of booleans                    
                    if sum(gtaic) >= lagvars
                        break
                    else
                        % Get residuals covariance matrix
                        [~,SIGMAfull{k},E{k}] = this.dataToVar(X,k);
                        % Using the AIC from Duan (2014) (Eq 4.18)
                        %                     fprintf('Model order = %d\n',k);
                        % working out aic with just sigma from X to Y
                        aic(k) = log(det(SIGMAfull{k})) + (2*k*M^2)/N;
                        %                     aic(k) = (N-k)*log(det(SIGMAfull{k}(1))) + (N-k)*(2*k*M^2)/((N-k)-(k*M^2)-1);%Toolbox uses a corrected version
                    end
                end
                
                % Select k with lowest AIC value:
                [~,K] = min(real(aic));                
%                 fprintf('Chosen model order: K = %d\n',K);
                
                if GCplots < 2
                    figure;
                    plot(1:kMax,aic,'bo-')
                    title('Model Order Estimation: Akaike Information Criteria')
                    xlabel('Model Order')
                    ylabel('AIC')
                end
                GCplots = GCplots + 1;
            end
            % Full regression
            [~,SIGMAfull,E] = this.dataToVar(X,K);
            rssF = sum(E.^2);

            % Reduced regression and F calculation:
            vars = 1:M;
            cm = zeros(M,M);
            fStat = zeros(M,M);
            pval = zeros(M,M);
            lSigFull = log(diag(SIGMAfull));
            for r = 1:M
                % Regression including all variables except variable No. r
                omitVar = vars(vars~=r); % omit variable No. r
                [~, sigma,eR] = this.dataToVar(X(:,omitVar),K);
                lSig = log(diag(sigma));
                rssR = sum(eR.^2);
                for ind = 1:M-1
                    c = omitVar(ind);
                    %influence omitted variable has on variable no. c
                    cm(r,c) = lSig(ind)-lSigFull(c);
                    fStat(r,c) = (rssR(ind) - rssF(c))*(N-2*K-1)/(rssF(c)*K);
                    pval(r,c) = 1-fcdf(fStat(r,c),K,(N-2*K-1));
                end
            end
        end
                
        function [cm, pVal,fStat,K] = grangerCausWrap(this,X,varargin)
            %GRANGERCAUS Pairwise granger causality (nonconditional)
            
            [~,M] = size(X);  
            cm = zeros(M,M);
            pVal = zeros(M,M);
            fStat = zeros(M,M);
            for r = 1:M
                for c = 1:M
                    if c == r
                        % Skipping diagonals
                    else
                        if nargin>2
                            % if model order is specified
                            K = varargin{1};
                            [gccm, pval,fstat,~,~] = this.grangerCaus([X(:,r) X(:,c)],K);
                        else
                            [gccm, pval,fstat,K(r,c),~] = this.grangerCaus([X(:,r) X(:,c)]);
                        end
                        cm(r,c) = gccm(2,1);
%                         cm(c,r) = gccm(1,2);
                        pVal(r,c) = pval(2,1);
%                         pVal(c,r) = pval(1,2);
                        fStat(r,c) = fstat(2,1);
%                         fStat(c,r) = fstat(1,2);
                    end
                end
            end
        end
                
        function [gcCM, pVal,fStat,K,sampleSize] = sampleSizeAnalysis(this,X,Y,nMin)
            N = length(X);
           
            nRange = nMin:N;
            parfor n = 1:length(nRange)
                disp(n)
                [cm, pval,fstat,K(n),~] = this.grangerCaus([X(1:nRange(n)),Y(1:nRange(n))]);
                sampleSize(n) = length(X(1:nRange(n)));
               gcCM(n) = cm(2,1);
               pVal(n) = pval(2,1);
               fStat(n) = fstat(2,1);               
            end
        end
        
        function refine(this,cm,pVals,alpha)
            %REFINE Remove non-significant values in the connectivity
            %matrix
            
            %Unpack
%             cm = this.ConnMatrix;
%             pVals = this.PVals;
            
            %Remove values where p>alpha
            sigInd = zeros(size(cm));
            sigInd(pVals<alpha) = 1;
            newCM = cm;
            newCM(~sigInd) = 0;
            
            %Package results
            this.Parameters.Alpha = alpha;
            this.NewCM = newCM;
                        
        end
        
    end
    
    methods (Static)
        
        function [ aic, aicOrder ] = orderEstimation(X,kMax)
            %ORDERESTIMATION Estimation of model order using Akaike information
            %criteria
            
            [N,M] = size(X);
            aic = zeros(kMax,1);
            for k = 1:kMax
                % Get residuals covariance matrix
                [~,SIGMAfull,~] = dataToVar(X,k);
                % Using the AIC from Duan (2014) (Eq 4.18)
                fprintf('Model order = %d\n',k);
                aic(k) = log(det(SIGMAfull)) + (2*k*M^2)/(N-k);%Toolbox uses a corrected version
            end
            [~,aicOrder] = min(aic);
            figure
            plot(1:kMax,aic,'bo-')
            title('Model Order Estimation: Akaike Information Criteria')
            xlabel('Model Order')
            ylabel('AIC')
        end
        
        function [A,SIGMA,E] = dataToVar(X,K)
            %
            %   Detailed explanation goes here
            
            [N,M] = size(X);
            K1 = K+1;
            X0 = X(K1:end,:);
            XL = zeros(N-K,M,K); % lagged X values
            for k = 1:K
                XL(:,:,k) = X(K1-k:N-k,:);
            end
            XL = reshape(XL,N-K,M*K);
            A = X0'/XL';
            E = X0'-A*XL';
            SIGMA = (E*E')/(N-K-1);
            E = E';
            A = reshape(A,M,M,K);
            
        end
        
        
    end
    
end

