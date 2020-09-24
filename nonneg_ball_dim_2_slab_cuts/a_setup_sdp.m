% Do not put 'clear' at the top of this file!

% Option 'bootstrap' is set by calling code

% Set settings for YALMIP

mysettings = sdpsettings('verbose', 0);

% Set constants (do I need all of these?)

n = 2;
e = ones(n, 1);
z = zeros(n, 1);
I = eye(n);
Z = zeros(n+1);

% Here is the problem in our notation:
%
%     min   x'Hx + 2g'x
%     st.   ||x|| <= 1
%           ||x|| <= e'x

% Set our parameters

r = 0;
R = 1;
a = 0;
b = e;
c = z;

% Set variables

x = sdpvar(n,1);
X = sdpvar(n,n,'symm');
Y = [1, x'; x, X];

% Set Shor constraints

con = [];
con = [con; Y >= 0];
con = [con; trace(X) <= R^2];
con = [con; trace(X) - 2*c'*x + c'*c <= b'*X*b - 2*a*b'*x + a^2];
con = [con; 0 <= b'*x - a];

% Set RLT and SOCRLT constraints (these lines assume n=2)

if bootstrap == 2 | bootstrap == 3
    con = [con; X(1,2) >= 0];
end

if bootstrap == 3
    con = [con; norm(X(1,:)) <= x(1)];
    con = [con; norm(X(2,:)) <= x(2)];
end

% Setup KSOC constraint

if bootstrap == 4

    M11 = R * [ b'*x - a, (x - c)';
                x - c,    (b'*x - a) * I ];

    Mx1 = [ b'*X(:,1) - a*x(1), (X(:,1) - x(1)*c)';
            X(:,1) - x(1)*c,    (b'*X(:,1) - a*x(1)) * I ];

    Mx2 = [ b'*X(:,2) - a*x(2), (X(:,2) - x(2)*c)';
            X(:,2) - x(2)*c,    (b'*X(:,2) - a*x(2)) * I ];

    M = [ M11, Mx1, Mx2;
          Mx1, M11, Z;
          Mx2, Z,   M11];

    con = [con; M >= 0];

end
