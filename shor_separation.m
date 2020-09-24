function prob = shor_separation(prob)
%SYNTAX: prob = shor_separation(prob)
%
%DESCRIPTION: 
%Find nonnegative functions l(x) and q(x) for cut by solving the "primal"

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

flterm = (r+R)*(b'*x - a) - r*R + c'*x + cmax*R;
glterm = (r+R)*(X*b - a*x) - r*R*x + X*c + cmax*R*x;
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
    prob.cuts.count = prob.cuts.count + 1;
    prob.cuts.Jl{prob.cuts.count} = Jl;
    prob.cuts.Jq{prob.cuts.count} = Jq;
    prob.cuts.fl{prob.cuts.count} = fl;
    prob.cuts.gl{prob.cuts.count} = gl;
    prob.cuts.fq{prob.cuts.count} = fq;
    prob.cuts.gq{prob.cuts.count} = gq;
    prob.cuts.Hq{prob.cuts.count} = Hq;
    prob.cuts.qlmin{prob.cuts.count} = qlmin;
    prob.cuts.primVal{prob.cuts.count} = primVal;
    prob.cuts.notSeparated = 0;
    prob.cuts.status{prob.cuts.count} = cvx_status;
    prob.cuts.notSeparated = 0;
    prob.cuts.method{prob.cuts.count} = 'Primal';
else
    prob.cuts.notSeparated = 1;
    prob.cuts.notSeparatedStatus = cvx_status;
    if strcmp(cvx_status,'Failed')
        warning('Optimization in separation failed. Try dual.')
    end
    if prob.options.verbose
        fprintf('Primal (%.3e) value not negative\n', primVal)
    end
end

