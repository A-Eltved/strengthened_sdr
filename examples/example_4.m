clear;

n = 2;
a = -0.77;
b = zeros(n, 1);
c = [-0.38; 0.18];
r = 0;
R = 1;
H = [-1.32, 0.21; 0.21, -0.81];
g = [-0.25; 0.05];
f = 0;
cmax = 0;

bootstrap = 4;

setup_shor_ksoc;

% Solve Shor \cap KSOC

nost = solvesdp(con, obj, mysettings);
if nost.problem ~= 0
    fprintf('Warning! YALMIP return code ~= 0\n');
end

double(obj)

% double(Y)

% Add cut

Hq = -4.9035 * eye(n);
gq = [-1.8633; 0.8826];
fq = 2.0403;

gl = [1.8633; -0.8826];
fl = 4.1236;

qlmin = 1.2604;

add_cut;

% Re-solve

nost = solvesdp(con, obj, mysettings);
if nost.problem ~= 0
    fprintf('Warning! YALMIP return code ~= 0\n');
end

double(obj)
double(Y)
