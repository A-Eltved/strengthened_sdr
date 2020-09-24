function prob = shor_ksoc_separation_dual(prob)
%SYNTAX: prob = shor_ksoc_separation_dual(prob)
%
%DESCRIPTION: 
%Find nonnegative functions l(x) and q(x) for cut
% by dual of the dual approach. 
% Recovery is ad hoc and not proven. 


n = prob.data.n;
x = prob.sol.Y.x;
X = prob.sol.Y.X;
R = prob.data.R;
Rsq = R*R;
r = prob.data.r;
a = prob.data.a;
b = prob.data.b;
c = prob.data.c;
cmax = prob.data.cmax;
xhat = prob.data.xhat;

rsq= r*r;
bT = prob.data.b';
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

Cq = Rsq*[1,x';x,X];

flterm = (r+R)*(b'*x - a) - r*R + c'*x + cmax*R;
glterm = (r+R)*(X*b - a*x) - r*R*x + X*c + cmax*R*x;
Cl = [flterm,glterm';glterm,zeros(n)];

Yhat = [1;xhat]*[1;xhat]';

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
    variable tl(1) nonnegative
    variable tq(1) nonnegative
    variable tql(1) nonnegative
    variable Yl(n+1,n+1) symmetric
    dual variable Zl
    dual variable al{4}
    variable Yq(n+1,n+1) symmetric
    dual variable Zq
    dual variable aq{4}
    variable Yql(n+1,n+1) symmetric
    dual variable Zql
    dual variable aql{4}
    variable Z(n,n) symmetric
    % KSOC
    variable Bq(n+1,n+1,n+1)
    variable Bqmat((n+1)^2,(n+1)^2) symmetric
    dual variable Wq
    variable Bl(n+1,n+1,n+1)
    variable Blmat((n+1)^2,(n+1)^2) symmetric
    dual variable Wl
    variable Bql(n+1,n+1,n+1)
    variable Bqlmat((n+1)^2,(n+1)^2) symmetric
    dual variable Wql

    maximize -tq - tl
    subject to 
        Cq == Yq + Yql - tq*Yhat
        Cl == Yl + [0,zeros(1,n); zeros(n,1), Z] + Yql - tl*Yhat
        -trace(X) == -Yql(1,1) + tql
        for i = 1:4
            al{i} : trace(A{i}*Yl) >= 0
            aq{i} : trace(A{i}*Yq) >= 0
            aql{i} : trace(A{i}*Yql) >= 0
        end
        Zl : Yl <In> semidefinite(n+1)
        Zq : Yq <In> semidefinite(n+1)
        Zql : Yql <In> semidefinite(n+1)
        % KSOC
        Wq : Bqmat(:,:) <In> semidefinite((n+1)^2)
        Wl : Blmat(:,:) <In> semidefinite((n+1)^2)
        Wql : Bqlmat(:,:) <In> semidefinite((n+1)^2)
        for i = 1:(n+1)^2
            for j = 1:(n+1)^2 
                Bqmat(i,j) == trace(Bij{i,j}*Yq);
                Blmat(i,j) == trace(Bij{i,j}*Yl);
                Bqlmat(i,j) == trace(Bij{i,j}*Yql);
            end
        end
        % Bq(1,1,1) == R*( bT*Yq(2:end,1) - a*Yq(1,1) )
        % Bl(1,1,1) == R*( bT*Yl(2:end,1) - a*Yl(1,1) )
        % Bql(1,1,1) == R*( bT*Yql(2:end,1) - a*Yql(1,1) )
        % for i = 2:n+1
        %     Bq(i,i,1) == R*( bT*Yq(2:end,1) - a*Yq(1,1) )
        %     Bl(i,i,1) == R*( bT*Yl(2:end,1) - a*Yl(1,1) )
        %     Bql(i,i,1) == R*( bT*Yql(2:end,1) - a*Yql(1,1) )
        % end
        % for i = 2:n
        %     Bq(i,i+1,1) == 0
        %     Bq(i+1,i,1) == 0
        %     Bl(i,i+1,1) == 0
        %     Bl(i+1,i,1) == 0
        %     Bql(i,i+1,1) == 0
        %     Bql(i+1,i,1) == 0
        % end
        % for i = 1:n
        %     for j = 1:n+1
        %         Bq(j,j,i+1) == R*( Yq(i+1,1) - c(i)*Yq(1,1) ) 
        %         Bl(j,j,i+1) == R*( Yl(i+1,1) - c(i)*Yl(1,1) ) 
        %         Bql(j,j,i+1) == R*( Yql(i+1,1) - c(i)*Yql(1,1) ) 
        %     end
        % end
        % for i = 1:n
        %     for j = 2:n
        %         Bq(j+1,j,i+1) == 0
        %         Bq(j,j+1,i+1) == 0
        %         Bl(j+1,j,i+1) == 0
        %         Bl(j,j+1,i+1) == 0
        %         Bql(j+1,j,i+1) == 0
        %         Bql(j,j+1,i+1) == 0
        %     end
        % end
        % % first row
        % for i = 1:n
        %     Bq(1,i+1,1) == bT*Yq(2:end,i+1) - a*Yq(i+1,1)
        %     Bq(i+1,1,1) == bT*Yq(2:end,i+1) - a*Yq(i+1,1)
        %     Bl(1,i+1,1) == bT*Yl(2:end,i+1) - a*Yl(i+1,1)
        %     Bl(i+1,1,1) == bT*Yl(2:end,i+1) - a*Yl(i+1,1)
        %     Bql(1,i+1,1) == bT*Yql(2:end,i+1) - a*Yql(i+1,1)
        %     Bql(i+1,1,1) == bT*Yql(2:end,i+1) - a*Yql(i+1,1)
        % end
        % for i = 1:n
        %     for j = 1:n
        %         Bq(1,j+1,i+1) == Yq(i+1,j+1) -   c(i)* Yq(j+1,1)
        %         Bq(j+1,1,i+1) == Yq(i+1,j+1) -   c(i)* Yq(j+1,1)
        %         Bl(1,j+1,i+1) == Yl(i+1,j+1) -   c(i)* Yl(j+1,1)
        %         Bl(j+1,1,i+1) == Yl(i+1,j+1) -   c(i)* Yl(j+1,1)
        %         Bql(1,j+1,i+1) == Yql(i+1,j+1) - c(i)*Yql(j+1,1)
        %         Bql(j+1,1,i+1) == Yql(i+1,j+1) - c(i)*Yql(j+1,1)
        %     end
        % end
        % % populate block diagonal
        % for i = 1:(n+1)
        %     Bqmat((i-1)*(n+1)+1:i*(n+1),(i-1)*(n+1)+1:i*(n+1)) == Bq(:,:,1)
        %     Blmat((i-1)*(n+1)+1:i*(n+1),(i-1)*(n+1)+1:i*(n+1)) == Bl(:,:,1)
        %     Bqlmat((i-1)*(n+1)+1:i*(n+1),(i-1)*(n+1)+1:i*(n+1)) == Bql(:,:,1)
        % end
        % % populate row
        % for i = 1:n
        %     Bqmat(i*(n+1)+1:(i+1)*(n+1),1:(n+1)) == Bq(:,:,i+1)
        %     Blmat(i*(n+1)+1:(i+1)*(n+1),1:(n+1)) == Bl(:,:,i+1)
        %     Bqlmat(i*(n+1)+1:(i+1)*(n+1),1:(n+1)) == Bql(:,:,i+1)
        % end
        % % populate zeros
        % for i = 1:n
        %     for j = (i+1):n 
        %         Bqmat(i*(n+1)+1:(i+1)*(n+1),j*(n+1)+1:(j+1)*(n+1)) == zeros(n+1) 
        %         Blmat(i*(n+1)+1:(i+1)*(n+1),j*(n+1)+1:(j+1)*(n+1)) == zeros(n+1) 
        %         Bqlmat(i*(n+1)+1:(i+1)*(n+1),j*(n+1)+1:(j+1)*(n+1)) == zeros(n+1) 
        %     end
        % end
cvx_end

% % Matrices for KSOC
% MB0diag = [-a*R, R/2*bT; R/2*b, zeros(n)];
% MB0r = cell(n,1);
% MBdiag = cell(n,1);
% MBr = cell(n,1);
% for i = 1:n
%     MB0r{i} = [R*c(i), R/2*evec(i,n)'; R/2*evec(i,n), zeros(n)];
%     MBdiag{i} = 1/2*[0, -a*evec(i,n)'; -a*evec(i,n), b*evec(i,n)' + evec(i,n)*bT];
%     MBr{i} = cell(n,1); 
%     for j = 1:n
%         MBr{i}{j} = 1/2*[0, -c(j)*evec(i,n)';
%         -c(j)*evec(i,n), evec(j,n)*evec(i,n)' + evec(i,n)*evec(j,n)'];
%     end
% end


% ad hoc recovery of primal variables
Jl = Zl + al{1}*A{1} + al{2}*A{2} + al{3}*A{3} + al{4}*A{4};
for i = 1:(n+1)^2
    for j = 1:(n+1)^2
        Jl = Jl + Wl(i,j)*Bij{i,j};
    end
end

Jq = Zq + aq{1}*A{1} + aq{2}*A{2} + aq{3}*A{3} + aq{4}*A{4};
for i = 1:(n+1)^2
    for j = 1:(n+1)^2
        Jq = Jq + Wq(i,j)*Bij{i,j};
    end
end

Jql = Zql + aql{1}*A{1} + aql{2}*A{2} + aql{3}*A{3} + aql{4}*A{4};
for i = 1:(n+1)^2
    for j = 1:(n+1)^2
        Jql = Jql + Wql(i,j)*Bij{i,j};
    end
end

% Extract function matrices and vectors
fl = Jl(1,1);
gl = Jl(2:end,1);
fq = Jq(1,1);
gq = Jq(2:end,1);
Hq = Jq(2:end,2:end);

qlmin = fl + fq - Jql(1,1);

primVal = trace(Cq*Jq) + trace(Cl*Jl) - trace(X)*qlmin;
if primVal < -prob.options.cuttol && primVal >= cvx_optval - 1e-1
    prob.cuts.count = prob.cuts.count + 1;
    prob.cuts.Jl{prob.cuts.count} = Jl;
    prob.cuts.Jq{prob.cuts.count} = Jq;
    prob.cuts.fl{prob.cuts.count} = fl;
    prob.cuts.gl{prob.cuts.count} = gl;
    prob.cuts.fq{prob.cuts.count} = fq;
    prob.cuts.gq{prob.cuts.count} = gq;
    prob.cuts.Hq{prob.cuts.count} = Hq;
    prob.cuts.qlmin{prob.cuts.count} = qlmin;
    prob.cuts.dualVal{prob.cuts.count} = cvx_optval;
    prob.cuts.primVal{prob.cuts.count} = primVal;
    prob.cuts.status{prob.cuts.count} = cvx_status;
    prob.cuts.notSeparated = 0;
    prob.cuts.method{prob.cuts.count} = 'Dual';
else
    prob.cuts.notSeparated = 1;
    prob.cuts.notSeparatedStatus = cvx_status;
    if strcmp(cvx_status,'Failed')
        warning('Optimization in separation failed. Try dual.')
    end
    if prob.options.verbose
        fprintf('Primal (%.3e) value not negative\n', primVal)
    end
end

end

function ei = evec(i,n)
ei = zeros(n,1);
ei(i) = 1;
end
