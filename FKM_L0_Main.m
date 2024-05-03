

clear all
clc
close all

%% Generate data

% clear
% N=115; T=2; C=4;
% unità 111-115 non chiaramente assegnabili ad alcun cluster
rng(123);
X1=randn(25,2)-5;
X2=randn(30,2);
X3=randn(35,2)+5;
X4=[randn(20,1)+5,randn(20,1)-5];
X12=randn(2,2)-2.5;
X23=randn(3,2)+2.5;
X=[X1;X2;X3;X4;X12;X23];
x1=X(:,1);
x2=X(:,2);
% plot(x1(1:25),x2(1:25),'*',x1(26:55),x2(26:55),'*',x1(56:90),x2(56:90),'*',x1(91:110),x2(91:110),'*',x1(111:112),x2(111:112),'o',x1(113:115),x2(113:115),'o');

% size(X)

%% Set and run FKM

m = 2;            % Fuzzyfing parameter (>1)
conv = 1e-16;     % Convergency treshold
Max_C = 20;
Max_iter = 1000;
stand = 1;
alpha = 1;

[N, T] = size(X); 

%% Choose number of clusters

% validity_values = zeros(1, Max_C);
num_random_starts = 100; % Numero di partenze casuali

validity_values = zeros(num_random_starts, Max_C);
parfor i = 1:num_random_starts
    for k = 2:Max_C
    
    [U,P,J] = FKM(X,m,k,conv);
    
    % Xie & Beni index
    XB_result = XB(X, U, P, m);
    validity_values(i, k) = XB_result;
    end
end

% Find the minimum of the validity index
min_validity_values = min(validity_values(:, 2:end), [], 1);

[~, optimal_C_index] = min(min_validity_values);
optimal_C = optimal_C_index + 1;

fprintf('Optimal number of clusters using Xie & Beni: %d\n', optimal_C);


%% Choose penalty parameter

C = optimal_C;
[U,P,J] = FKM(X,m,C,conv);  

% Specify lambda values
lambda_values = 0:(max(J)/(100000)):(max(J)/100);

% Initialize variables to store Xie & Beni index values
silhouette_values = zeros(size(lambda_values));

% Loop over lambda values
for lambda_idx = 1:length(lambda_values)
    lambda = lambda_values(lambda_idx);

    rng('default');
    % Call the FKM_L0_Lambda function to obtain the optimal U and P
    [U, P, ~] = FKM_L0_Lambda(X, C, m, lambda, conv, Max_iter, stand, alpha);

    % Compute the Xie & Beni index
    silhouette_values(lambda_idx) = XB(X, U, P, m);
end

% Find the lambda value that minimizes the index
[optimal_silhouette, optimal_lambda_idx] = min(silhouette_values);
optimal_lambda_L0 = lambda_values(optimal_lambda_idx);


[UL0, PL0, JL0] = FKM_L0(X, C, m, optimal_lambda_L0, conv, Max_iter, stand);  % Compute FKM-L0  "Non-Exhaustive"

% [UL0P, PL0P, JL0P] = FKM_L0_P(X, C, m, optimal_lambda_L0, conv, Max_iter,
% stand);  % Compute FKM-L0 "Exhaustive"


[~, Out] = max(UL0, [], 2);

% save Results.mat;
