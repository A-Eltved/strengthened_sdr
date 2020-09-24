function prob = solve_ksoc_with_cuts(prob)
% SYNTAX: prob = solve_ksoc_with_cuts(prob)
%
% DESCRIPTION: 
% Solve Shor+KSOC relaxation with cuts added. 

prob = solve_sdr_shor_ksoc_cone(prob);
if prob.sol.Y.rat >= prob.options.rattol
    fprintf('Shor+KSOC is exact.\n')
    prob.isExact = 1;
    return
end
prob.isExact = 0;
prob.baseSol = prob.sol;
prob = compute_cmax(prob);

totalCuts = 0;
fprintf('Cuts: ')
while prob.sol.Y.rat < prob.options.rattol
    fprintf('%i, ',totalCuts)
    prob = shor_ksoc_separation(prob);

    if prob.cuts.notSeparated
        if strcmp(prob.cuts.notSeparatedStatus,'Failed')
            prob = shor_ksoc_separation_dual(prob);
        end
    end
    if prob.cuts.notSeparated
        prob = shor_ksoc_separation_r0(prob);

        if prob.r0cuts.notSeparated
            fprintf('\nPoint not separated.\n')
            prob = find_local_min(prob);
            prob.gapClosure = 1 - (prob.local.fval - prob.sol.val)/(prob.local.fval - prob.baseSol.val);
            prob.closed = false;
            return
        end
        fprintf('(r=0)->')
    end
    totalCuts = totalCuts + 1;

    prob = solve_sdr_shor_ksoc_cone_with_cuts(prob);
end
prob.closed = true;
prob.gapClosure = 1;
fprintf('\nRank-1 solution found!.\n')
