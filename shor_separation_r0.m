function prob = shor_separation(prob)
%SYNTAX: prob = shor_separation(prob)
%
%DESCRIPTION: 
%Find nonnegative functions l(x) and q(x) for cut by solving the "primal".
%use r=0 as stated in corrollary 1 (paper version 2020-08-09)

n = prob.data.n;
x = prob.sol.Y.x;
X = prob.sol.Y.X;
R = prob.data.R;
Rsq = R*R;
r = prob.data.r;
a = prob.data.a;
b = prob.data.b;
c = prob.data.c;
cmax = prob.data.cmax;
xhat = prob.data.xhat;

rsq= r*r;
bT = prob.data.b';
cT = prob.data.c';
asq= a*a;
bbT= prob.data.b*prob.data.b';
cTc= prob.data.c'*prob.data.c; 

A = cell(4);
A{1} = [-rsq, zeros(1,n); zeros(n,1), eye(n)];
A{2} = [Rsq, zeros(1,n); zeros(n,1), -eye(n)];
A{3} = [asq - cTc, cT - a*bT; c - a*b, bbT - eye(n)];
A{4} = [-a, bT/2; b/2, zeros(n)];

Cq = Rsq*[1,x';x,X];

flterm = R*(b'*x - a) + c'*x + cmax*R;
glterm = R*(X*b - a*x) + X*c + cmax*R*x;
Cl = [flterm,glterm';glterm,zeros(n)];

Yhat = [1;xhat]*[1;xhat]';

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
    variable zl(4) nonnegative
    variable zq(4) nonnegative
    variable zql(4) nonnegative
    variable Jl(n+1,n+1) symmetric
    variable Jq(n+1,n+1) symmetric
    variable Zl(n+1,n+1) semidefinite
    variable Zq(n+1,n+1) semidefinite
    variable Zql(n+1,n+1) semidefinite
    variable qlmin(1) nonnegative

    minimize trace(Cq*Jq) + trace(Cl*Jl) - trace(X)*qlmin
    subject to 
        Jq == Zq + zq(1)*A{1} + zq(2)*A{2} +zq(3)*A{3} +zq(4)*A{4} 
        Jl == Zl + zl(1)*A{1} + zl(2)*A{2} +zl(3)*A{3} +zl(4)*A{4} 
        Jl + Jq - [qlmin, zeros(1,n); zeros(n,1), zeros(n)] == Zql + zql(1)*A{1} + zql(2)*A{2} +zql(3)*A{3} +zql(4)*A{4} 
        Jl(2:end,2:end) == zeros(n)
        trace(Yhat*Jq) <= 1
        trace(Yhat*Jl) <= 1
cvx_end

% disp(Jl)
% disp(Jq)
% disp(Jl + Jq - [qlmin, zeros(1,n); zeros(n,1), zeros(n)])

primVal = cvx_optval;
if primVal < -prob.options.cuttol
    % Extract function matrices and vectors
    fl = Jl(1,1);
    gl = Jl(2:end,1);
    fq = Jq(1,1);
    gq = Jq(2:end,1);
    Hq = Jq(2:end,2:end);
    prob.r0cuts.count = prob.r0cuts.count + 1;
    prob.r0cuts.Jl{prob.r0cuts.count} = Jl;
    prob.r0cuts.Jq{prob.r0cuts.count} = Jq;
    prob.r0cuts.fl{prob.r0cuts.count} = fl;
    prob.r0cuts.gl{prob.r0cuts.count} = gl;
    prob.r0cuts.fq{prob.r0cuts.count} = fq;
    prob.r0cuts.gq{prob.r0cuts.count} = gq;
    prob.r0cuts.Hq{prob.r0cuts.count} = Hq;
    prob.r0cuts.qlmin{prob.r0cuts.count} = qlmin;
    prob.r0cuts.primVal{prob.r0cuts.count} = primVal;
    prob.r0cuts.notSeparated = 0;
else
    if prob.options.verbose
        disp('Primal value not negative')
    end
    prob.r0cuts.notSeparated = 1;
end

