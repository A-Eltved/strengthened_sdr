function prob = solve_sdr_shor_cone(prob)
% SYNTAX: prob = solve_sdr_shor_cone(prob)
%
% DESCRIPTION:
% Solve Shor semidefinite relaxation in cone representation. 



% Local variables 
n  = prob.data.n;
H  = prob.data.H;
gT = prob.data.g';
r  = prob.data.r;
rsq= r*r;
R  = prob.data.R;
Rsq= R*R;
a  = prob.data.a;
b  = prob.data.b;
bT = prob.data.b';
c  = prob.data.c;
% if sum(abs(c))>0
%     warning("c =/= 0 not tested properly.")
% end
cT = prob.data.c';
asq= a*a;
bbT= prob.data.b*prob.data.b';
cTc= prob.data.c'*prob.data.c; 

A = cell(4);
A{1} = [-rsq, zeros(1,n); zeros(n,1), eye(n)];
A{2} = [Rsq, zeros(1,n); zeros(n,1), -eye(n)];
A{3} = [asq - cTc, cT - a*bT; c - a*b, bbT - eye(n)];
A{4} = [-a, bT/2; b/2, zeros(n)];

% Build large matrix 
Bij = cell((n+1)^2);
% fill in zero matrices first
for i = 1:(n+1)^2
    for j = 1:(n+1)^2 % Only fill lower triangle, since it's symmetric
        Bij{i,j} = zeros(n+1);
    end
end
% Fill diagonal
for i = 1:(n+1)^2
    % fprintf('Diagonal, i = %i\n',i)
    Bij{i,i} = R/2*[-2*a, bT; b, zeros(n)];
    % disp(Bij{i,i})
end
% Fill column of block diagonal
for k = 0:n
    for j = 1:n
    % fprintf('Row of block diagonal, k = %i, j = %i\n',k,j)
    % disp([k*(n+1)+1+j, k*(n+1)+1])
        Bij{k*(n+1)+1+j, k*(n+1)+1} = 1/2*[0, -a*evec(j,n)';...
        -a*evec(j,n), evec(j,n)*bT + b*evec(j,n)'];
        Bij{k*(n+1)+1, k*(n+1)+1+j} = 1/2*[0, -a*evec(j,n)';...
        -a*evec(j,n), evec(j,n)*bT + b*evec(j,n)'];
        % disp(Bij{k*(n+1)+1, k*(n+1)+1+j})
    end
end
% Fill diagonal of column blocks
for k = 1:n
    for j = 1:(n+1)
        % fprintf('Diagonal of row matrices, k = %i, j = %i\n',k,j)
        % disp([k*(n+1) + j, j])
        Bij{k*(n+1) + j, j} = R/2*[-2*c(k), evec(k,n)'; evec(k,n), zeros(n)];
        Bij{j, k*(n+1) + j} = R/2*[-2*c(k), evec(k,n)'; evec(k,n), zeros(n)];
    end
end
% Fill first column of column blocks
for k = 1:n
    for j = 1:n
        % fprintf('Row of row matrices, k = %i, j = %i\n',k,j)
        % disp([k*(n+1)+1+j, 1])
        % disp(1/2*[0, -c(k)*evec(j,n)'; -c(k)*evec(j,n),...
        %                     evec(j,n)*evec(k,n)' + evec(k,n)*evec(j,n)'])
        Bij{k*(n+1)+1+j, 1} = 1/2*[0, -c(k)*evec(j,n)'; -c(k)*evec(j,n),...
                            evec(j,n)*evec(k,n)' + evec(k,n)*evec(j,n)'];
        Bij{k*(n+1)+1, 1+j} = 1/2*[0, -c(k)*evec(j,n)'; -c(k)*evec(j,n),...
                            evec(j,n)*evec(k,n)' + evec(k,n)*evec(j,n)'];
        Bij{1 ,k*(n+1)+1+j} = 1/2*[0, -c(k)*evec(j,n)'; -c(k)*evec(j,n),...
                            evec(j,n)*evec(k,n)' + evec(k,n)*evec(j,n)'];
        Bij{1+j,k*(n+1)+1} = 1/2*[0, -c(k)*evec(j,n)'; -c(k)*evec(j,n),...
                            evec(j,n)*evec(k,n)' + evec(k,n)*evec(j,n)'];
    end
end

if prob.options.verbose
    cvx_begin
else
    cvx_begin quiet
end
if strcmp(prob.options.solver,'sedumi')
    cvx_solver sedumi
end
if strcmp(prob.options.solver,'mosek')
    cvx_solver mosek
end
try % added later so older problems might not have it and should be solved with high precision
    if prob.options.highprecision
        if prob.options.verbose
            fprintf('Using high precision')
        end
        cvx_precision high
    end
catch
    if prob.options.verbose
        fprintf('Using high precision')
    end
    cvx_precision high
end
    variable x(n)
    variable X(n,n) symmetric
    variable Y(n+1,n+1) semidefinite
    % KSOC
    variable B(n+1,n+1,n+1)
    variable Bmat((n+1)^2,(n+1)^2) semidefinite

    minimize trace(H*X) +  2*gT*x 
    subject to 
        Y == [1,x'; x, X]
        for i = 1:4
            trace(A{i}*Y) >= 0
        end
        % KSOC
        % B(1,1,1) == R*( bT*Y(2:end,1) - a*Y(1,1) )
        % for i = 2:n+1
        %     B(i,i,1) == R*( bT*Y(2:end,1) - a*Y(1,1) )
        % end
        % for i = 2:n
        %     B(i,i+1,1) == 0
        %     B(i+1,i,1) == 0
        % end
        % for i = 1:n
        %     for j = 1:n+1
        %         B(j,j,i+1) == R*( Y(i+1,1) - c(i)*Y(1,1) ) 
        %     end
        % end
        % for i = 1:n
        %     for j = 2:n
        %         B(j+1,j,i+1) == 0
        %         B(j,j+1,i+1) == 0
        %     end
        % end
        % % first row
        % for i = 1:n
        %     B(1,i+1,1) == bT*Y(2:end,i+1) - a*Y(i+1,1)
        %     B(i+1,1,1) == bT*Y(2:end,i+1) - a*Y(i+1,1)
        % end
        % for i = 1:n
        %     for j = 1:n
        %         B(1,j+1,i+1) == Y(i+1,j+1) - c(i)*Y(j+1,1)
        %         B(j+1,1,i+1) == Y(i+1,j+1) - c(i)*Y(j+1,1)
        %     end
        % end
        % % populate block diagonal
        % for i = 1:(n+1)
        %     Bmat((i-1)*(n+1)+1:i*(n+1),(i-1)*(n+1)+1:i*(n+1)) == B(:,:,1)
        % end
        % % populate row
        % for i = 1:n
        %     Bmat(i*(n+1)+1:(i+1)*(n+1),1:(n+1)) == B(:,:,i+1)
        % end
        % % populate zeros
        % for i = 1:n
        %     for j = (i+1):n 
        %         Bmat(i*(n+1)+1:(i+1)*(n+1),j*(n+1)+1:(j+1)*(n+1)) == zeros(n+1) 
        %     end
        % end
        for i = 1:(n+1)^2
            for j = 1:(n+1)^2 
                % fprintf('i = %i, j = %i\n',i,j)
                % disp(Bij{i,j})
                Bmat(i,j) == trace(Bij{i,j}*Y);
            end
        end
cvx_end
prob.sol.status = cvx_status;
evals = sort(eig(Y),'descend');

prob.sol.Y.val = Y;
prob.sol.Y.x = x;
prob.sol.Y.X = X;
prob.sol.Y.rat = evals(1)/abs(evals(2)); %putting absolute value, since eigenvalue can be slighly negative 
prob.sol.Y.evals = evals;
prob.sol.val = cvx_optval;
end

function ei = evec(i,n)
ei = zeros(n,1);
ei(i) = 1;
end
