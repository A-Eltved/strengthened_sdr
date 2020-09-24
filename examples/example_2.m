clear;

n = 2;
a = 0;
b = ones(n, 1);
c = zeros(n, 1);
r = 0;
R = 1;
H = -eye(n) - [1; 0] * [1; 1]';
g = 0.5 * ([1; 1] + [1; 0]);
f = 0;

bootstrap = 4;

setup_shor_ksoc;

solvesdp(con, obj, mysettings)

double(obj)
