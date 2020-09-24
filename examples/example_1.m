clear;

n = 2;
a = -1;
b = [-1; -1];
c = zeros(n, 1);
r = 0;
R = 1;
H = -eye(n);
g = 0.5 * b;
f = -a;

bootstrap = 4;

setup_shor_ksoc;

solvesdp(con, obj, mysettings)

double(obj)
