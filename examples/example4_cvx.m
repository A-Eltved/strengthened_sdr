clc
clear all
close all

addpath('..')
check_cvx

a = -0.77;
b = zeros(2,1);
c = [-0.38; 0.18];
H = [-1.32, 0.64; -0.22, -0.81];
H = (H + H')/2;
g = [-0.25; 0.05];
r = 0;
R = 1;

p = create_problem(H,g,r,R,a,b,c);

p.data.xhat = [-0.667118085937451; 0.733646700682918];
    
p = solve_ksoc_with_cuts(p);

fprintf('Cut:')
fprintf('fl:')
disp(p.cuts.fl{1})
fprintf('gl:\n')
disp(p.cuts.gl{1})
fprintf('[q+l]_{min}:')
disp(p.cuts.qlmin{1})
fprintf('fq:')
disp(p.cuts.fq{1})
fprintf('gq:\n')
disp(p.cuts.gq{1})
fprintf('Hq:\n')
disp(p.cuts.Hq{1})

fprintf('Shor+KSOC solution: %.4f\n',p.baseSol.val)
fprintf('Shor+KSOC+Cut solution: %.4f\n',p.sol.val)

fprintf('Engenvalue ratio: %.2e\n',p.sol.Y.rat)
fprintf('Y(x,X):\n')
disp(p.sol.Y.val)

plot_prob(p);
