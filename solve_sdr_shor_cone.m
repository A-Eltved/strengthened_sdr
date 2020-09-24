function prob = solve_sdr_shor_cone(prob)
% SYNTAX: prob = solve_sdr_shor_cone(prob)
%
% DESCRIPTION:
% Solve Shor semidefinite relaxation in cone representation. 

% Local variables 
n  = prob.data.n;
H  = prob.data.H;
gT = prob.data.g';
r  = prob.data.r;
rsq= r*r;
R  = prob.data.R;
Rsq= R*R;
a  = prob.data.a;
b  = prob.data.b;
bT = prob.data.b';
c  = prob.data.c;
cT = prob.data.c';
asq= a*a;
bbT= prob.data.b*prob.data.b';
cTc= prob.data.c'*prob.data.c; 

A = cell(4);
A{1} = [-rsq, zeros(1,n); zeros(n,1), eye(n)];
A{2} = [Rsq, zeros(1,n); zeros(n,1), -eye(n)];
A{3} = [asq - cTc, cT - a*bT; c - a*b, bbT - eye(n)];
A{4} = [-a, bT/2; b/2, zeros(n)];

if prob.options.verbose
    cvx_begin
else
    cvx_begin quiet
end
if strcmp(prob.options.solver,'sedumi')
    cvx_solver sedumi
end
if strcmp(prob.options.solver,'mosek')
    cvx_solver mosek
end
try % added later so older problems might not have it and should be solved with high precision
    if prob.options.highprecision
        if prob.options.verbose
            fprintf('Using high precision')
        end
        cvx_precision high
    end
catch
    if prob.options.verbose
        fprintf('Using high precision')
    end
    cvx_precision high
end
    variable x(n)
    variable X(n,n) semidefinite

    minimize trace(H*X) +  2*gT*x 
    subject to 
        for i = 1:4
            trace(A{i}*[1,x'; x, X]) >= 0
        end
        [1,x'; x, X] <In> semidefinite(n+1)
cvx_end
prob.sol.status = cvx_status;
Y = [1,x';x,X];
evals = sort(eig(Y),'descend');

prob.sol.Y.val = Y;
prob.sol.Y.x = x;
prob.sol.Y.X = X;
prob.sol.Y.rat = evals(1)/abs(evals(2)); %putting absolute value, since eigenvalue can be slighly negative 
prob.sol.Y.evals = evals;
prob.sol.val = cvx_optval;
