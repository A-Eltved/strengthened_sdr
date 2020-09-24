function prob = shor_separation_dual(prob)
%SYNTAX: prob = shor_separation_dual(prob)
%
%DESCRIPTION: 
%Find nonnegative functions l(x) and q(x) for cut
% by dual of the dual approach. For now, let's see if they can help.
% Address how to recover the functions for the problem from the dual later

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
    variable tl(1) nonnegative
    variable tq(1) nonnegative
    variable tql(1) nonnegative
    % variable Yl(n+1,n+1) semidefinite
    % variable Yq(n+1,n+1) semidefinite
    % variable Yql(n+1,n+1) semidefinite
    variable Yl(n+1,n+1) symmetric
    dual variable Zl
    dual variable al{4}
    variable Yq(n+1,n+1) symmetric
    dual variable Zq
    dual variable aq{4}
    variable Yql(n+1,n+1) symmetric
    dual variable Zql
    dual variable aql{4}
    variable Z(n,n) symmetric

    maximize -tq - tl
    subject to 
        Cq == Yq + Yql - tq*Yhat
        Cl == Yl + [0,zeros(1,n); zeros(n,1), Z] + Yql - tl*Yhat
        -trace(X) == -Yql(1,1) + tql
        for i = 1:4
            al{i} : trace(A{i}*Yl) >= 0
            aq{i} : trace(A{i}*Yq) >= 0
            aql{i} : trace(A{i}*Yql) >= 0
        end
        Zl : Yl <In> semidefinite(n+1)
        Zq : Yq <In> semidefinite(n+1)
        Zql : Yql <In> semidefinite(n+1)
cvx_end
Jl = Zl + al{1}*A{1} + al{2}*A{2} + al{3}*A{3} + al{4}*A{4};
Jq = Zq + aq{1}*A{1} + aq{2}*A{2} + aq{3}*A{3} + aq{4}*A{4};
Jql = Zql + aql{1}*A{1} + aql{2}*A{2} + aql{3}*A{3} + aql{4}*A{4};

% Extract function matrices and vectors
fl = Jl(1,1);
gl = Jl(2:end,1);
fq = Jq(1,1);
gq = Jq(2:end,1);
Hq = Jq(2:end,2:end);

qlmin = fl + fq - Jql(1,1);

primVal = trace(Cq*Jq) + trace(Cl*Jl) - trace(X)*qlmin;
if primVal < -prob.options.cuttol
    prob.cuts.count = prob.cuts.count + 1;
    prob.cuts.Jl{prob.cuts.count} = Jl;
    prob.cuts.Jq{prob.cuts.count} = Jq;
    prob.cuts.fl{prob.cuts.count} = fl;
    prob.cuts.gl{prob.cuts.count} = gl;
    prob.cuts.fq{prob.cuts.count} = fq;
    prob.cuts.gq{prob.cuts.count} = gq;
    prob.cuts.Hq{prob.cuts.count} = Hq;
    prob.cuts.qlmin{prob.cuts.count} = qlmin;
    prob.cuts.dualVal{prob.cuts.count} = cvx_optval;
    prob.cuts.primVal{prob.cuts.count} = primVal;
    prob.cuts.status{prob.cuts.count} = cvx_status;
    prob.cuts.notSeparated = 0;
    prob.cuts.method{prob.cuts.count} = 'Dual';
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

