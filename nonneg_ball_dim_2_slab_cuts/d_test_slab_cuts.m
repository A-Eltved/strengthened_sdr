clear

% Initialize values

tau_rank = 1.0e4;
tau_sep = -1.0e-5;
max_cuts = 25;

% Set the bootstrap relaxation

% bootstrap = 1; % Shor
% bootstrap = 2; % Shor + RLT
% bootstrap = 3; % Shor + RLT + SOCRLT
bootstrap = 4; % Shor + KSOC

% Setup SDP 

a_setup_sdp

% Load instances *not* solved by KSOC (saved by random number seed)

files = dir('mat_files/*.mat');

problems_solved = 0;

rank_meas_0 = [];
rank_meas_1 = [];

for iter = 1:min(100, length(files))

    load(strcat('mat_files/', files(iter).name));
    % Set seed

    % rng(not_solved_by_KSOC(iter));

    % Setup objective

    % This code assumes n=2
    % tmp = randn(5,1);
    % tmp = tmp/norm(tmp);
    % H = zeros(n); g = zeros(n,1);
    % H(1,1) = tmp(1);
    % H(1,2) = tmp(2);
    % H(2,1) = tmp(2);
    % H(2,2) = tmp(3);
    % g(1) = tmp(4);
    % g(2) = tmp(5);

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
        fprintf('obj = %.3e\n', double(obj));
        evals = sort(eig(double(Y)), 'descend');
        rm = evals(1) / abs(evals(2));
        rank_meas_0 = [rank_meas_0; rm];
        if rm > tau_rank
            % Problem solved
        else
            % evals(1) / abs(evals(2))
        end
    end

    dX = double(X);
    dx = double(x);
    y = dX*e - dx;
    II = find(y >= 0);
    JJ = setdiff(1:n, II)';
    viol_II = e'*dx - norm(y(II)) - trace(dX); 
    viol_JJ = 1.0 - norm(y(JJ)) - trace(dX);

    cuts_added = 0;
    cuts = [];

    while cuts_added <= max_cuts & (viol_II < tau_sep | viol_JJ < tau_sep)

        if viol_II < tau_sep
            s = z;
            s(II) = y(II) / norm(y(II));
            cuts = [cuts; trace(X) <= e'*x - s'*(X*e - x)];
            cuts_added = cuts_added + 1;
        end

        if viol_JJ < tau_sep
            s = z;
            s(JJ) = - y(JJ) / norm(y(JJ));
            cuts = [cuts; trace(X) <= 1.0 + s'*(X*e - x)];
            cuts_added = cuts_added + 1;
        end

        nost = solvesdp([con, cuts], obj, mysettings);
        if nost.problem ~= 0
            fprint('Mosek not happy\n');
        end

        dX = double(X);
        dx = double(x);
        y = dX*e - dx;
        II = find(y >= 0);
        JJ = setdiff(1:n, II)';
        viol_II = e'*dx - norm(y(II)) - trace(dX);
        viol_JJ = 1.0 - norm(y(JJ)) - trace(dX);

    end

    fprintf('us  = %.3e\n', double(obj));
    evals = sort(eig(double(Y)), 'descend');
    rm = evals(1) / abs(evals(2));
    rank_meas_1 = [rank_meas_1; rm];
    if rm > tau_rank
        fprintf('problem solved (cuts added = %d)\n', cuts_added);
        problems_solved = problems_solved + 1;
    else
        % evals(1) / abs(evals(2))
    end

end

