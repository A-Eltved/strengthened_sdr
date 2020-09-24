function prob = create_problem(H,g,r,R,a,b,c,options)
% SYNTAX: prob = create_problem(H,g,r,R,a,b,c,options)
%
% DESCRIPTION:
% Creates structure with problem data for the problem
%   minimize    x^T H x + 2 g^T x
%   subject to  r <= ||x|| <= R
%               ||x - c|| <= b^T x - a
%
% INPUT: 
%   nxn matrices:   H
%   nx1 vectors:    g, c, b
%   scalars:         r, R, a

prob = struct();
% Info about problem
prob.data.n = length(g);
% Objective
prob.data.H = H;
prob.data.g = g;
% Constraints
prob.data.r = r;
prob.data.R = R;
prob.data.a = a;
prob.data.b = b;
prob.data.c = c;
% Fields for cuts
prob.cuts.count = 0;
prob.cuts.fl = [];
prob.cuts.gl = [];
prob.cuts.fq = [];
prob.cuts.gq = [];
prob.cuts.Hq = [];
prob.cuts.qlmin = [];
prob.cuts.dualVal = [];
prob.cuts.primVal = [];
% Fields for r=0 cuts
prob.r0cuts.count = 0;
prob.r0cuts.fl = [];
prob.r0cuts.gl = [];
prob.r0cuts.fq = [];
prob.r0cuts.gq = [];
prob.r0cuts.Hq = [];
prob.r0cuts.qlmin = [];
prob.r0cuts.dualVal = [];
prob.r0cuts.primVal = [];
% Options
% default options
prob.options.use_kron = true;
prob.options.verbose = false;
prob.options.rattol = 1e4;
prob.options.cuttol = 1e-6;
prob.options.solver = 'mosek';
prob.options.highprecision = false;
% Adjust to user options
if nargin == 8
    if isfield(options,'verbose'); prob.options.verbose = options.verbose; end
    if isfield(options,'rattol'); prob.options.rattol = options.rattol; end
    if isfield(options,'cuttol'); prob.options.cuttol = options.cuttol; end
    if isfield(options,'solver'); prob.options.solver = options.solver; end
    if isfield(options,'highprecision'); prob.options.highprecision = options.highprecision; end
end
end
