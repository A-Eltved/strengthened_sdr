% Run the "verify_chenetal" code to setup symbolic quantities, which
% will be made concrete below

verify_chenetal

% Set specific parameter values to test

valL11 =  0; valU11 = 1;
valL12 = -1; valU12 = 1;
valL22 =  0; valU22 = 1;

% Or uncomment the next few lines to generate random values for the
% parameters

% valL11 =  rand; valU11 = rand;
% valL12 = -rand; valU12 = rand;
% valL22 =  rand; valU22 = rand;
%
% if valL11 > valU11;
%     tmp = valU11;
%     valU11 = valL11;
%     valL11 = tmp;
% end
%
% if valL22 > valU22;
%     tmp = valU22;
%     valU22 = valL22;
%     valL22 = tmp;
% end

% Setup structures to do the symbolic substitutions below

old = {L11, L12, L22, U11, U12, U22};
new = {valL11, valL12, valL22, valU11, valU12, valU22};

% First get numerical values of b1 and b2

b1 = double(subs(valb(1), old, new));
b2 = double(subs(valb(2), old, new));

% Next get nuemerical values of p0 through pi4

pi0 = double(subs(pi0, old, new));
pi1 = double(subs(pi1, old, new));
pi2 = double(subs(pi2, old, new));
pi3 = double(subs(pi3, old, new));
pi4 = double(subs(pi4, old, new));

% Finally, replace the main parameter values. Everything should be
% concrete values now, not symbols

L11 = double(subs(L11, valL11));
U11 = double(subs(U11, valU11));
L12 = double(subs(L12, valL12));
U12 = double(subs(U12, valU12));
L22 = double(subs(L22, valL22));
U22 = double(subs(U22, valU22));

% Setup constants to build model. Compared to the paper, this code is a
% little confusing because, in this example, the norm is not taken with
% respect to all variables, just the first 2 out of 3

% Overall dimension is 3, but broken into two subsets, J and K

n = 3;
J = [1, 2];
K = [3];

% Setup identiy and zero matrices (different dimensions from one another)

I = eye(2);
Z = zeros(3);

% Setup parameters to match paper. Note here the relevant dimension is 2

r = sqrt(L11);
R = sqrt(U11);
a = 0;
b = [b1; b2];
c = [0; 0];

% Setup SDP variables

x = sdpvar(n,1);
X = sdpvar(n,n,'symm');
Y = [1, x'; x, X];

% Set Shor constraints. This includes the bounds on x(3)

con = [];
con = [con; Y >= 0];
con = [con; r^2 <= trace(X(J,J)); trace(X(J,J)) <= R^2];
con = [con; trace(X(J,J)) - 2*c'*x(J) + c'*c <= b'*X(J,J)*b - 2*a*b'*x(J) + a^2];
con = [con; 0 <= b'*x(J) - a];
con = [con; sqrt(L22) <= x(K); x(K) <= sqrt(U22)];

% Add RLT constraint coming from the bounds on x(3)

con = [con; X(K,K) + sqrt(L22) * sqrt(U22) <= (sqrt(L22) + sqrt(U22)) * x(K)];

% Add SOCRLT constraints, which link x(1:2) with x(3). I hope this is right!

con = [con; norm(sqrt(U22) * x(J) - X(J,K)) <= R * (sqrt(U22) - x(K))];
con = [con; norm(sqrt(U22) * x(J) - X(J,K)) <= b' * (sqrt(U22) * x(J) - X(J,K))];

con = [con; norm(X(J,K) - sqrt(L22) * x(J)) <= R * (x(K) - sqrt(L22))];
con = [con; norm(X(J,K) - sqrt(L22) * x(J)) <= b' * (X(J,K) - sqrt(L22) * x(J))];

% Add KSOC constraint

M11 = R * [ b'*x(J) - a, (x(J) - c)';
            x(J) - c,    (b'*x(J) - a) * I ];

Mx1 = [ b'*X(J,1) - a*x(1), (X(J,1) - x(1)*c)';
        X(J,1) - x(1)*c,    (b'*X(J,1) - a*x(1)) * I ];

Mx2 = [ b'*X(J,2) - a*x(2), (X(J,2) - x(2)*c)';
        X(J,2) - x(2)*c,    (b'*X(J,2) - a*x(2)) * I ];

M = [ M11, Mx1, Mx2;
      Mx1, M11, Z;
      Mx2, Z,   M11];

con = [con; M >= 0];

% Setup objective. See paper

obj = (pi0 + U11 * U22) + ...
    (pi1 - U22) * (X(1,1) + X(2,2)) + ...
    (pi2 - U11) * X(3,3) + ...
    pi3 * X(1,3) + ...
    pi4 * X(2,3);

% Solve SDP. A negative objective value means that the Chen-et-al cut is
% not captured by Shor \cap KSOC

nost = solvesdp(con, obj)
fprintf('KSOC: mosek = %f\n', double(obj));
