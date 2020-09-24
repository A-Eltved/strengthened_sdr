close all
clear all
clc

addpath('..')
check_cvx
% Example 1 from paper
% minimize    -(x1^2 + x2^2) - x1 - x2
% subject to  ||x|| <= 1
%             ||x|| <= -x1 -x2 -1

% Create instance 
n = 2;
R = 1; 
r = 0;
a = -1;
b = -ones(n,1);
c = 0*ones(n,1);

H = -eye(n);
g = [-0.5; -0.5];

options = struct();
prob = create_problem(H,g,r,R,a,b,c,options);

% Manually set interior point
prob.data.xhat = zeros(n,1);

plot_prob(prob); 

% Shor+KSOC separation
fprintf('Solve Shor+KSOC relaxation...\n')
prob = solve_sdr_shor_ksoc_cone(prob);
fprintf('Objective value: %.4f\n',prob.sol.val)
fprintf('NOTE: The constant 1 should be added to obtain the actual objective, since constants are not considered.\n')

