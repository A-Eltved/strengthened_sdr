clear;

n = 2;
a = -1;
b = -ones(n, 1);
c = zeros(n, 1);
r = 0;
R = 1;
H = -eye(n);
g = [-0.55; -0.5];
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
double(Y)

% Add first cut

Hq = -0.3812 * eye(n);
gq = [-0.5578; -0.5531];
fq = 0.8563;

gl = [0.3462; 0.3608];
fl = 1;

qlmin = 1.42;

add_cut;

% Add second cut

Hq = [-0.7065, 0.1719; 0.1719, -0.4368];
gq = [-0.7808; -0.7278];
fq = 1;

gl = [0.3442; 0.3626];
fl = 1;
qlmin = 1.155;

add_cut;

% Add third cut

Hq = [-0.6296, 0.2398; 0.2398, -0.4512];
gq = [-0.7868; -0.7580];
fq = 1;

gl = [0.3479; 0.3591];
fl = 1;

qlmin = 1.149;

add_cut;

% Re-solve

nost = solvesdp(con, obj, mysettings);
if nost.problem ~= 0
    fprintf('Warning! YALMIP return code ~= 0\n');
end

double(obj)
double(Y)
