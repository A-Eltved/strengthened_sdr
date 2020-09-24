function prob = instance_generator_TTRS(n,options)
% SYNTAX: prob = instance_generator_TTRS(n,options)
%
% DESCRIPTION:
% Generates TTRS instance (b=r=0) with non-empty interior

if nargin < 2
    options = struct();
end

beta = 1; % Width of distribution for a
R = 1; % We might as well fix this
r = 0; 
% r = rand*R;
xt = randn(n,1);
rt = r + (R-r)*rand;
xhat = xt/norm(xt,2)*rt;
% b = randn(n,1);
b = zeros(n,1);
c = randn(n,1);
amax = b'*xhat - norm(xhat - c,2);
a = amax - beta + beta*rand;

H = randn(n);
H = (H + H')/2;
g = randn(n,1);

prob = create_problem(H,g,r,R,a,b,c,options);

prob.data.xhat = xhat;


