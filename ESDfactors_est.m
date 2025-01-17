function est = ESDfactors_est(Y,W,Xy,Z,Xz,Ry,Rz,x0,itr_times,corr_times,message,options)
%ESDfactors_estimation
%   Inputs:
%       Y: n*(T+1) matrix of outcomes
%       W: n*n*(T+1) spatial weight matrix.
%       Xy: n*(T+1)*ky matrix. Independent variables in main outcome equation
%       Z: n*p*(T+1) matrix. Dependent variable in auxiliary equation
%       Xz: n*(T+1)*max_kzp*p matrix. Independent variables in auxiliary
%           equation. The third dimension equals to the maximum number of
%           independent variable for p. For those with lest than max_kzp,
%           keep NaN in data.
%       Ry: number of factors in y
%       Rz: number of factors in z
%       x0: initial guess for the optimization problem. If one does not
%           want to assign initial guess, keep this argument as [].
%       itr_times:  default value is 9. If x0 is [], then the algorithm will
%                   generate "itr_times" number of start points to avoid 
%                   local optimal solution
%       corr_times: default value is 1. Times for bias correction. After 
%                   each bias correction, the estimator could be used to 
%                   calculate more precise asymptotic variance, and based 
%                   on this asymptotic variance, the route can update the 
%                   bias-corrected estimator.
%       message: Ture or flase. If Ture, the route will generate message 
%                notifying the user which step is on.
%       options: optimal algorithm options.
%
%
%
%   Outputs: est - struct with fields: param, nb, bc, exitflag
%       est.param: saves numbers of factors, n, T, numbers of parameters in
%                  the model.
%       est.nb: saves estimators before bias correction.
%       est.bc: saves estimators after bias correction.
%       est.exitflag: exit condition, same as for fminunc.

arguments
    Y
    W
    Xy
    Z
    Xz
    Ry              double {mustBeInteger}
    Rz              double {mustBeInteger}
    x0 = []
    itr_times       double {mustBeInteger} = 9
    corr_times      double {mustBeInteger} = 1
    message                         = isempty(getCurrentTask())
    options.TolX = 1e-9                % The termination tolerance for x
    options.MaxFunEvals = 100000       % bound on the number of function evaluations
    options.MaxIter = 100000           % bound on the number of solver iterations
    options.TolFun = 1e-9
    options.Display = 'notify'         % displays output only if the function does not converge.
    options.PlotFcn = []  %{@optimplotx,@optimplotfval,@optimplotfirstorderopt}
    options.UseParallel = isempty(getCurrentTask())  % use parallel by default, but if doing simulation, close.
end




%% Get rid of the first observations.
Y0 = Y(:,1);
W0 = W(:,:,1);
Z0 = Z(:,:,1);

Y(:,1) = [];
W(:,:,1) = [];
Xy(:,1,:) = [];
Z(:,:,1) = [];
Xz(:,1,:,:) = [];

%% get model parameters
[n, T, ky] = size(Xy);
[~,~,max_kzp,p] = size(Xz);
J = p*(p+1)/2;

kzp = nan(p,1);
for j = 1 : p
    if sum(isnan(Xz(:,:,:,j)),'all') == 0
        kzp(j) = max_kzp;
    else
        indices = find( isnan(Xz(:,:,:,j)) );
        [~, ~, temp_c] = ind2sub(size(Xz(:,:,:,j)), indices);
        kzp(j) = max(temp_c)-1;
    end
end
kz = sum(kzp);
K = [ky, 3, p, kz, p^2, 1, J];



est.param.Ry = Ry;
est.param.Rz = Rz;
est.param.n = n;
est.param.T = T;
est.param.ky = ky;
est.param.max_kzp = max_kzp;
est.param.p = p;
est.param.J = J;
est.param.kzp = kzp;
est.param.kz = kz;
est.param.K = K;



%% estimate
opt = optimoptions('fminunc',   ...
    'TolX',options.TolX,...
    'MaxFunEvals',options.MaxFunEvals, ...
    'MaxIter', options.MaxIter, ...
    'TolFun', options.TolFun, ...
    'Display', options.Display, ...
    'Algorithm', 'quasi-newton',...
    'PlotFcn', options.PlotFcn,  ...
    'UseParallel', options.UseParallel);


if isempty(x0)

    fval = inf(itr_times,1);
    x = nan(sum(est.param.K),itr_times);
    exitflag=nan(itr_times,1);

    if message
        fprintf('Estimating...  0%%');
    end

    for itr = 1 : itr_times
        x0 = randn(sum(est.param.K),1);
        
        while true
            try

                percentDone = 100 * itr / itr_times ;

                if message
                    fprintf('\b\b\b\b%3.0f%%', percentDone);  % Display progress as percentage
                end

                [temp_x, fval_temp, temp_exitflag, ~] ...
                    = fminunc(@(arg) obj_Q(arg,Y,W,Xy,Z,Xz,Y0,Z0,W0, Ry, Rz),x0,opt);

                if fval_temp < fval(itr)
                    fval(itr) = fval_temp;
                    x(:,itr) = temp_x;
                    exitflag(itr) = temp_exitflag;
                end

                break

            catch

                warning('initial guess not feasible, try another one...')
            
            end
        end
    end

    if message
        fprintf('\b\b\b\b%3.0f%%', 100);  % Display progress as percentage
        fprintf('\n');  % in-Command Wondow waitbar ends
    end

    [est.nb.fval, index] = min(fval);
    final_x = x(:,index);
    est.exitflag = exitflag(index);

else
    if message
        fprintf('Estimating...  \n');
    end
    [final_x, est.nb.fval, est.exitflag, ~] ...
            = fminunc(@(arg) obj_Q(arg,Y,W,Xy,Z,Xz,Y0,Z0,W0, Ry, Rz),x0,opt);

end

%%% adjust some estimators
final_x(ky+1:ky+3) = tanh(final_x(ky+1:ky+3));
final_x(ky+3+p+kz+p^2+1) = exp(final_x(ky+3+p+kz+p^2+1)); % alpha_xi
% for alpha, since Sigma_epsilon^(-1) = Sigma_epsilon_epsilon^(-1/2)*Sigma_epsilon_epsilon^(-1/2)
% where Sigma_epsilon_epsilon^(-1/2) make differs in sign, so we use
% sqrtm() function to make sure the diagonal elements are positive.

temp_alpha = final_x(ky+3+p+kz+p^2+1+1 : ky+3+p+kz+p^2+1+J);
Sigma_epsilon_neg_half_temp = zeros(p,p);
Sigma_epsilon_neg_half_temp(tril(true(p,p))) = temp_alpha;
Sigma_epsilon_neg_half_temp = Sigma_epsilon_neg_half_temp + tril(Sigma_epsilon_neg_half_temp,-1)';
Sigma_epsilon_temp = (Sigma_epsilon_neg_half_temp)^(-2);
Sigma_epsilon_neg_half_temp = sqrtm(inv(Sigma_epsilon_temp));

final_x(ky+3+p+kz+p^2+1+1 : ky+3+p+kz+p^2+1+J) = nonzeros(tril(Sigma_epsilon_neg_half_temp));
est.nb.theta = final_x;

clear temp_alpha Sigma_epsilon_neg_half_temp Sigma_epsilon_temp final_x

% retrive estimators
est.nb = theta_to_params(est.nb, est.param);

% check stationary condtion
[c, ~] = stability_cond(est.nb.theta, W, est);
if sum(c > 0) ~= 0  % constrains not satisfied, we use fmincon
    warning('Stability condition not satisfied for unbias correced estimators.')
end
clear c



%% estimate factors and factor loadings and define the projection matrix
[Hy, Hz]  = H(est.nb.theta,Y,W,Xy,Z,Xz,Y0,Z0,W0,est);
est.nb = factors_and_proj_matrix(est.nb, Hz, Hy, est.param);



%% original higher moments of epsilon and xi, variance and bias corrector
if message
    fprintf('Calculating original asymptotic variance...    ');
end
est.nb = moments_and_residuals(est.nb,est.param,Y,W,Xy,Z,Xz,Y0,Z0,W0);

est.nb = asymp_var(est.nb, est.param, Y, W, Xy, Z, Xz, Y0, Z0, W0, message);
est.nb = rmfield(est.nb, {'P_Fz', 'M_Fz', 'P_Gz', 'M_Gz', 'P_Fy', 'M_Fy', 'P_Gy', 'M_Gy'});
est.nb.bcorr = (est.nb.Sigma \ est.nb.phi) / sqrt(n*T); % bias-corrector



%% bias corrected estimators
est.bc.theta = est.nb.theta - est.nb.bcorr;
est.bc = theta_to_params(est.bc, est.param);

% varaince matrix of est.bc.theta
[Hy, Hz]  = H(est.bc.theta, Y,W,Xy,Z,Xz,Y0,Z0,W0,est);
est.bc = factors_and_proj_matrix(est.bc, Hz, Hy, est.param);

if message
    fprintf('Calculating bias corrected asymptotic variance (%d/%d)...    ', 1, corr_times);
end

est.bc = moments_and_residuals(est.bc,est.param,Y,W,Xy,Z,Xz,Y0,Z0,W0);
est.bc = asymp_var(est.bc, est.param, Y, W, Xy, Z, Xz, Y0, Z0, W0, message);
est.bc = rmfield(est.bc, {'P_Fz', 'M_Fz', 'P_Gz', 'M_Gz', 'P_Fy', 'M_Fy', 'P_Gy', 'M_Gy'});


% update bias corrected estimator based on better phi
est.bc.bcorr = (est.bc.Sigma \ est.bc.phi) / sqrt(n*T); % bias corrector
est.bc.theta = est.nb.theta - est.bc.bcorr;
est.bc = theta_to_params(est.bc, est.param);



% we have already corrected theta above (for one time), so i start at 2
for i = 2 : corr_times

    [Hy, Hz]  = H(est.bc.theta, Y,W,Xy,Z,Xz,Y0,Z0,W0,est);
    est.bc = factors_and_proj_matrix(est.bc, Hz, Hy, est.param);


    if message
        fprintf('Calculating bias corrected asymptotic variance (%d/%d)...    ', i, corr_times);
    end
    est.bc = moments_and_residuals(est.bc,est.param,Y,W,Xy,Z,Xz,Y0,Z0,W0);
    est.bc = asymp_var(est.bc, est.param, Y, W, Xy, Z, Xz, Y0, Z0, W0, message);
    est.bc = rmfield(est.bc, {'P_Fz', 'M_Fz', 'P_Gz', 'M_Gz', 'P_Fy', 'M_Fy', 'P_Gy', 'M_Gy'});


    % update bias corrected estimator based on better phi
    est.bc.bcorr = (est.bc.Sigma \ est.bc.phi) / sqrt(n*T);
    est.bc.theta = est.nb.theta - est.bc.bcorr;
    est.bc = theta_to_params(est.bc, est.param);

end

est.bc.fval = obj_Q(est.bc.theta,Y,W,Xy,Z,Xz,Y0,Z0,W0,est.param.Ry,Rz);


end
