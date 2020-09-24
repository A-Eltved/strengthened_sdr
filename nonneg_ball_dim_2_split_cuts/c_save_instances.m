clear

n = 2;

% Load instances *not* solved by KSOC (saved by random number seed)

if ~isfile('b_search_for_instances.mat')
    fprintf('Error! First search for instances!\n');
end

load('b_search_for_instances.mat');

for iter = 1:length(not_solved_by_KSOC)

    % Set seed

    rng(not_solved_by_KSOC(iter));

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

    mystr = sprintf('mat_files/instance_seed_%010d.mat', not_solved_by_KSOC(iter));
    save(mystr, 'H', 'g');

end

