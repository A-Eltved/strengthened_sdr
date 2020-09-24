function prob = instance_generator_r0_c0_a0(n,options)
% SYNTAX: prob = instance_generator(n,options)
%
% DESCRIPTION:
% Generates instance with non-empty interior

if nargin < 2
    options = struct();
end

beta = 1; % Width of distribution for a
R = 1; % We might as well fix this
r = 0; 
b = ones(n,1)*(1/sqrt(n) + rand*2*n);
c = zeros(n,1);
amax = 0;
a = 0;

xhat = 1/2*b/norm(b,2); 

H = randn(n);
H = (H + H')/2;
g = randn(n,1);

prob = create_problem(H,g,r,R,a,b,c,options);

prob.data.xhat = xhat;


