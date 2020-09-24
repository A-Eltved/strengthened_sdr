% Analyze cell of problems

nprobs = length(probs);
n = probs{1}.data.n;
fprintf('Dimension: %i\n',n)
fprintf('Number of problems: %i\n',nprobs)
% Exactness
shorExact = false(nprobs,1);
ksocExact = false(nprobs,1);
shorCutExact = true(nprobs,1);
ksocCutExact = true(nprobs,1);
% Number of cuts
ncutsShor = zeros(nprobs,1);
ncutsKSOC = zeros(nprobs,1);
% Gap improvement
objShor = zeros(nprobs,1);
objShorCuts = zeros(nprobs,1);
objKSOC = zeros(nprobs,1);
objKSOCCuts = zeros(nprobs,1);
localMin = zeros(nprobs,1);
gapClosureShor = zeros(nprobs,1);
gapClosureKSOC = zeros(nprobs,1);

idclosedKSOC = false(nprobs,1);
idclosedShor = false(nprobs,1);


for i = 1:nprobs 
    if isfield(probs{i}.pshor,'isExact')
        if probs{i}.pshor.isExact
            objShor(i) = probs{i}.pshor.sol.val;
            shorExact(i) = true;
            ksocExact(i) = true;
        else
            baseobjShor(i) = probs{i}.pshor.baseSol.val;
            objShor(i) = probs{i}.pshor.sol.val;
            if isfield(probs{i}.pshor,'gapClosure')
                if ~isnan(probs{i}.pshor.gapClosure)
                    gapClosureShor(i) = probs{i}.pshor.gapClosure;
                else
                    % disp('Separation not solved. Using KSOC local to calculate gap...')
                    if probs{i}.pksoc.isExact
                        loc = probs{i}.pksoc.sol.val;
                    else
                        loc = probs{i}.pksoc.local.fval;
                    end
                    probs{i}.pshor.gapClosure = 1 - (loc - probs{i}.pshor.sol.val)/(loc - probs{i}.pshor.baseSol.val);
                    probs{i}.pshor.closed = 0;
                end
            end
                
            ncutsShor(i) = probs{i}.pshor.cuts.count;
            if probs{i}.pshor.r0cuts.count > 0
                ncuts(i) = ncuts(i) + probs{i}.pshor.r0cuts.count;
            end
            if probs{i}.pshor.sol.Y.rat >= probs{i}.options.rattol
                idclosedShor(i) = true;
            else
                shorCutExact(i) = false;
            end
        end
    end
    if probs{i}.pksoc.isExact
        baseobjKSOC(i) = probs{i}.pksoc.sol.val;
        objKSOC(i) = probs{i}.pksoc.sol.val;
        ksocExact(i) = true;
    else
        if ~isfield(probs{i}.pksoc,'gapClosure')
            probs{i}.pksoc = find_local_min(probs{i}.pksoc);
            probs{i}.pksoc.gapClosure = 1 - (probs{i}.pksoc.local.fval - probs{i}.pksoc.sol.val)/(probs{i}.pksoc.local.fval - probs{i}.pksoc.baseSol.val);
            probs{i}.pksoc.closed = false;
        end
        gapClosureKSOC(i) = probs{i}.pksoc.gapClosure;
        baseobjKSOC(i) = probs{i}.pksoc.baseSol.val;
        objKSOC(i) = probs{i}.pksoc.sol.val;
        if probs{i}.pksoc.sol.Y.rat >= probs{i}.pksoc.options.rattol
            idclosedKSOC(i) = true;
        else
            ksocCutExact(i) = false;
        end
        ncutsKSOC(i) = probs{i}.pksoc.cuts.count;
    end
end


% Extract some relevant problems
ninex = sum(~ksocCutExact);
pinex = cell(ninex,1);
tmp = 1;
for i = find(~ksocCutExact)'
    pinex{tmp} = probs{i};
    tmp = tmp +1;
end
nimp = sum(ncutsKSOC > 0);
pimp = cell(nimp,1);
tmp = 1;
for i = find(ncutsKSOC > 0)'
    pimp{tmp} = probs{i};
    tmp = tmp +1;
end
nnc = sum(~idclosedKSOC & ~ksocExact);
pnc = cell(nnc,1);
tmp = 1;
for i = find(~idclosedKSOC & ~ksocExact)'
    pnc{tmp} = probs{i};
    tmp = tmp +1;
end


tab{1,1} = nprobs - sum(shorExact);
tab{1,2} = sum(ncutsShor>0) - sum(idclosedShor);
tab{1,3} = mean(ncutsShor(ncutsShor > 0 & ~idclosedShor));
tab{1,4} = mean(gapClosureShor(ncutsShor > 0 & ~idclosedShor));
tab{1,5} = sum(idclosedShor);
tab{1,6} = mean(ncutsShor(idclosedShor));

fprintf('Number of inexact Shor relaxations: %i\n',tab{1,1})
fprintf('Number of problems where cuts help Shor: %i\n',tab{1,2})
fprintf('Average number of cuts (Shor): %.2f\n',tab{1,3})
fprintf('Average gap clusure (Shor): %.2f\n',tab{1,4})
fprintf('Number of closed gaps (Shor): %i\n',tab{1,5})
fprintf('Average number of cuts (Shor): %.2f\n',tab{1,6})   

tab2{1,1} = nprobs - sum(ksocExact);
tab2{1,2} = sum(ncutsKSOC>0) - sum(idclosedKSOC);
tab2{1,3} = mean(ncutsKSOC(ncutsKSOC > 0 & ~idclosedKSOC));
tab2{1,4} = mean(gapClosureKSOC(ncutsKSOC > 0 & ~idclosedKSOC));
tab2{1,5} = sum(idclosedKSOC);
tab2{1,6} = mean(ncutsKSOC(idclosedKSOC));

fprintf('Number of inexact KSOC relaxations: %i\n',tab2{1,1})
fprintf('Number of problems where cuts help KSOC: %i\n',tab2{1,2})
fprintf('Average number of cuts (KSOC): %.2f\n',tab2{1,3})
fprintf('Average gap clusure (KSOC): %.2f\n',tab2{1,4})
fprintf('Number of closed gaps (KSOC): %i\n',tab2{1,5})
fprintf('Average number of cuts (KSOC): %.2f\n',tab2{1,6})
