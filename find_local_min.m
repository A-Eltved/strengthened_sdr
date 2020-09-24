function prob = find_local_min(prob)
% SYNTAX: prob = find_local_min(prob)
%
% DESCRIPTION
% Find local minimum.
options = optimoptions('fmincon','Display','off');
nlocal = 100;

objfun = @(x) objective_fun(x,prob);
cons = @(x) nonlcons(x,prob);
for i = 1:nlocal
    x0 = randn(prob.data.n,1);
    [xs(:,i),fvals(i),EXITFLAG] = fmincon(objfun,x0,[],[],[],[],[],[],cons,options);
    if EXITFLAG ~= 1
        fvals(i) = NaN;
    end
end
[fval,id] = min(fvals);

prob.local.x = xs(:,id);
prob.local.fval = fval;

end

function [f,g] = objective_fun(x,prob)

f = x'*prob.data.H*x + 2*prob.data.g'*x;

if nargout > 1 % gradient required
    g = 2*prob.data.H*x + 2*prob.data.g;
end
end

function [c,ceq] = nonlcons(x,prob)
c = [norm(x) - prob.data.R; 
    prob.data.r - norm(x);  
    norm(x - prob.data.c) - prob.data.b'*x + prob.data.a];
ceq = [];
end
