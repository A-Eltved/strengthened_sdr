function check_cvx(cvx_path)
% Checks if cvx is loaded
% adds cvx_path if not, and runs startup
if ~exist('cvx_begin')
    disp('Loading CVX...')
    if nargin < 1
        cvx_path = '~/Documents/MATLAB/cvx';
    end
    addpath(cvx_path)
    cvx_startup
end

