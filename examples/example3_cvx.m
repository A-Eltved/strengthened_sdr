close all
clear all
clc

addpath('..')
check_cvx
% Example 3 from paper
% minimize    -(x1^2 + x2^2) - 1.1*x1 - x2
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
g = [-0.55; -0.5];

options = struct();

prob = create_problem(H,g,r,R,a,b,c,options);

% Manually set interior point
prob.data.xhat = zeros(n,1);

plot_prob(prob); 
% Solve Shor relaxation
pshor = solve_sdr_shor_cone(prob);
if pshor.sol.Y.rat < pshor.options.rattol
    fprintf('Running Shor separation\n')
    pshor = compute_cmax(pshor);
    pshor.baseSol = pshor.sol;
    cuts = 0;
    fprintf('Cuts: ')
    while cuts == pshor.cuts.count
        fprintf(' %i,',cuts)
        pshor = shor_separation(pshor);
        pshor = solve_sdr_shor_cone_with_cuts(pshor);
        if pshor.sol.Y.rat >= pshor.options.rattol
            fprintf('\nRank-1 solution found!\n')
            break
        end
        cuts = cuts + 1;
    end
    fprintf('\n')
end

prob = solve_sdr_shor_ksoc_cone(prob);

% Shor+KSOC separation
if prob.sol.Y.rat < prob.options.rattol
    fprintf('Running Shor+KSOC separation\n')
    prob = compute_cmax(prob);
    prob.baseSol = prob.sol;
    cuts = 0;
    fprintf('Cuts: ')
    while cuts == prob.cuts.count
        fprintf(' %i,',cuts)
        prob = shor_ksoc_separation(prob);
        prob = solve_sdr_shor_ksoc_cone_with_cuts(prob);
        if prob.sol.Y.rat >= prob.options.rattol
            fprintf('\nRank-1 solution found!\n')
            break
        end
        cuts = cuts + 1;
    end
    fprintf('\n')
end
