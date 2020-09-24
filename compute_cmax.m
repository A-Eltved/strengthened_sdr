function prob = compute_cmax(prob)
% SYNTAX: prob = compute_cmax(prob)
%
% DESCRIPTION
% Computes [c]_{max} by bisection

epsi = 1e-5;

% Copy data to local variables
n = prob.data.n;
r = prob.data.r;
R = prob.data.R;
a = prob.data.a;
b = prob.data.b;
c = prob.data.c;
bT = b';
cT = c';

l = 0; 
u = norm(c,2);

while u-l > epsi
    mu = l + (u-l)/2;
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
        variable x(n,1)
        minimize mu*sum_square_abs(x) - r*cT*x
        subject to
            norm(x,2) <= R
            norm(x-c,2) <= bT*x-a
    cvx_end
    v = cvx_optval;
    if v < 0
        l = mu;
    else
        u = mu;
    end
end

prob.data.cmax = u;
