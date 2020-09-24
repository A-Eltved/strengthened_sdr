clear

% Set the bootstrap relaxation (in search mode, want KSOC)

% bootstrap = 1; % Shor
% bootstrap = 2; % Shor + RLT
% bootstrap = 3; % Shor + RLT + SOCRLT
bootstrap = 4; % Shor + KSOC

% Setup constants

tau_sep = -1.0e-8;
tau_rank = 10^5;
max_cuts = 25;
max_iter = 100;

% Setup SDP 

a_setup_sdp

% Initialize array to hold indices of subproblems *not* solved by KSOC

if isfile('b_search_for_instances.mat')
     load('b_search_for_instances.mat')
else
    furthest_searched = 0;
    not_solved_by_KSOC = [];
end

upper_limit = furthest_searched + max_iter;

for iter = (furthest_searched + 1) : upper_limit

    % Print progress

    if mod(iter, 100) == 0
        fprintf('iter = %d (out of %d)\n', iter, upper_limit);
    end

    % Set seed

    rng(iter);

    % Setup objective

    % This code assumes n=2
    tmp = randn(5,1);
    tmp = tmp/norm(tmp);
    H = zeros(n); g = zeros(n,1);
    H(1,1) = tmp(1);
    H(1,2) = tmp(2);
    H(2,1) = tmp(2);
    H(2,2) = tmp(3);
    g(1) = tmp(4);
    g(2) = tmp(5);

    % H = randn(n); H = 0.5*(H + H');
    % g = randn(n, 1);

    obj = H(:)'*X(:) + 2*g'*x;

    % obj = e'*x - X(1,:)*e + x(1) - trace(X);

    % Solve SDP

    nost = solvesdp(con, obj, mysettings);
    if nost.problem ~= 0
        fprint('Mosek not happy\n');
        return
    else
        % fprintf('obj = %.3e\n', double(obj));
        evals = sort(eig(double(Y)), 'descend');
        if evals(1) / abs(evals(2)) > tau_rank
            % Solved
        else
            % Not solved
            not_solved_by_KSOC = [not_solved_by_KSOC; iter];
        end
    end

    furthest_searched = iter;
    save('b_search_for_instances.mat', 'not_solved_by_KSOC', 'furthest_searched');

end
