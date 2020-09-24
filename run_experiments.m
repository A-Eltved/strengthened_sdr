close all
clear all
clc

check_cvx; % add CVX to path

seed = 1331;
N = 15000; % number of experiments
n = 3; % dimension of problems

disp(n)

filename = sprintf('data/simulations/general__seed_%i__n_%i__N_%i',seed,n,N);

rng(seed)


% Create instances

options = struct();
probs = cell(N,1);
fprintf('Creating instances...\n')
for i = 1:N
    prob = instance_generator(n,options);
    
    probs{i} = prob;
end

fprintf('Solving shor relaxation...\n')
shor_exact = false(N,1);
for i = 1:N
    fprintf('%i/%i\n',i,N)
    tic
    probs{i} = solve_sdr_shor_cone(probs{i});
    probs{i}.shor_time = toc;
    if probs{i}.sol.Y.rat > probs{i}.options.rattol
        shor_exact(i) = true;
    end
end
nprobs = sum(~shor_exact);
fprintf('Pruning... Number of inexact problems: %i \n',nprobs)
probs = probs(~shor_exact); % Only keep problems that are not exact
fprintf('Solving Shor+KSOC and adding cuts...\n')
for i = 1:nprobs
    fprintf('%i/%i\n',i,nprobs)

    % Solve Shor relaxation
    tic
    probs{i}.pshor = solve_shor_with_cuts(probs{i});
    probs{i}.pshor.time = toc;

    %Solve Shor+KSOC
    tic
    probs{i}.pksoc = solve_ksoc_with_cuts(probs{i});
    probs{i}.pksoc.time = toc;
end

% Save to file
num = 0;
savefile = strcat(filename,'_',num2str(num),'.mat');
while exist(savefile)
    num = num+1;
    savefile = strcat(filename,'_',num2str(num),'.mat');
end
fprintf('Saving as: %s\n',savefile)
save(savefile,'probs','shor_exact')

