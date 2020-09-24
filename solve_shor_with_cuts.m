function prob = solve_shor_with_cuts(prob)
% SYNTAX: prob = solve_shor_with_cuts(prob)
%
% DESCRIPTION: 
% Solve Shor+Cuts. 

prob = solve_sdr_shor_cone(prob);
if prob.sol.Y.rat >= prob.options.rattol
    fprintf('Shor is exact.\n')
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
    try
        prob = shor_separation(prob);
    catch
        warning('Error in shor separation. Returning...')
        prob.separation_failed = true;
        fprintf('\nPoint not separated.\n')
        prob = find_local_min(prob);
        prob.gapClosure = 1 - (prob.local.fval - prob.sol.val)/(prob.local.fval - prob.baseSol.val);
        prob.closed = 0;
        return
    end

    if prob.cuts.notSeparated
        if strcmp(prob.cuts.notSeparatedStatus,'Failed')   
            try
                prob = shor_separation_dual(prob);
            catch
                warning('Error in shor dual separation. Returning...')
                prob.separation_failed = true;
                fprintf('\nPoint not separated.\n')
                prob = find_local_min(prob);
                prob.gapClosure = 1 - (prob.local.fval - prob.sol.val)/(prob.local.fval - prob.baseSol.val);
                prob.closed = 0;
                return
            end
        end
    end
    if prob.cuts.notSeparated
        prob = shor_separation_r0(prob);

        if prob.r0cuts.notSeparated
            fprintf('\nPoint not separated.\n')
            prob = find_local_min(prob);
            prob.gapClosure = 1 - (prob.local.fval - prob.sol.val)/(prob.local.fval - prob.baseSol.val);
            prob.closed = 0;
            return
        end
        fprintf('<-(r=0)')
    end
    totalCuts = totalCuts + 1;

    prob = solve_sdr_shor_cone_with_cuts(prob);
end
prob.closed = 1;
prob.gapClosure = 1;
fprintf('\nRank-1 solution found!.\n')
