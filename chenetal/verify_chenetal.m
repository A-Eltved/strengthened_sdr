% Setup symbolic parameters (both constants and variables)

syms L11 L22 L12 U11 U22 U12 % note: constants defined by Chen et al
syms X11 X22 X33 X13 X23     % note: there is no X12
syms b1 b2                   % note: these are our b1 and b2

% Setup the derivated constants established by Chen at al. These should
% match their notation exactly.

pi0 = -sqrt(L11)*sqrt(L22)*sqrt(U11)*sqrt(U22);
pi1 = -sqrt(L22)*sqrt(U22);
pi2 = -sqrt(L11)*sqrt(U11);

fL12 = (sqrt(1 + L12^2) - 1)/L12;
fU12 = (sqrt(1 + U12^2) - 1)/U12;

pi3 = (sqrt(L11) + sqrt(U11))*(sqrt(L22) + sqrt(U22))*(1 - fL12*fU12)/(1 + fL12*fU12);
pi4 = (sqrt(L11) + sqrt(U11))*(sqrt(L22) + sqrt(U22))*(fL12 + fU12)/(1 + fL12*fU12);

% Calculate the values of b (as a vector) according to our Proposition 4
% (as of 2020-07-29).

valb = inv([1, L12; 1, U12])*[sqrt(1 + L12^2); sqrt(1 + U12^2)];

%%% Part 1: Upper bound cut

% Setup the expression of our cut over several steps. We have expr >= 0.
% As of 2020-07-29, this explicit expression is found in the proof of
% Theorem 2 in Section 3, although it fundamentally comes from Theorem 1
% in Section 2.

expr = (U22 - X33) * sqrt(U11) + (sqrt(L22) + sqrt(U22)) * ( b1 * X13 + b2 * X23 );
expr = expr * ( sqrt(L11) + sqrt(U11) );
expr = expr - ( U22 + sqrt(L22) * sqrt(U22) ) * ( X11 + X22 + sqrt(L11) * sqrt(U11) );

% Substitute the values of b1 and b2 into expr.

expr = subs(expr, {b1, b2}, {valb(1), valb(2)});

% Create the "upper bound" cut by Chen et al symbollically. The form is
% chen >= 0.

chen = (pi0 + U11 * U22) + (pi1 - U22) * (X11 + X22) + (pi2 - U11) * X33;
chen = chen + pi3 * X13 + pi4 * X23;

% Simplify chen - expr. Should get 0

simplify(chen - expr)

%%% Part 2: Lower bound cut

% Setup the expression of our cut over several steps. We have expr >= 0.
% As of 2020-07-29, this explicit expression is found in the proof of
% Theorem 3 in Section 3, although it fundamentally comes from Proposition 2
% in Section 2.

expr = (sqrt(L22) + sqrt(U22)) * ( b1 * X13 + b2 * X23 ) - (X33 - L22) * sqrt(L11);
expr = expr * ( sqrt(L11) + sqrt(U11) );
expr = expr - ( L22 + sqrt(L22) * sqrt(U22) ) * ( X11 + X22 + sqrt(L11) * sqrt(U11) );

% Substitute the values of b1 and b2 into expr.

expr = subs(expr, {b1, b2}, {valb(1), valb(2)});

% Create the "lower bound" cut by Chen et al symbollically. The form is
% chen >= 0.

chen = (pi0 + L11 * L22) + (pi1 - L22) * (X11 + X22) + (pi2 - L11) * X33;
chen = chen + pi3 * X13 + pi4 * X23;

% Simplify chen - expr. Should get 0

simplify(chen - expr)