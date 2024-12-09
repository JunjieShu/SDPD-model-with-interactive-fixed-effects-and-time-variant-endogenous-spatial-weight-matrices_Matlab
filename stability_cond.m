function [c, ceq] = stability_cond(arg, W, est)
    % arg = est.nb.theta
    
    ky = est.param.ky;
    p = est.param.p;
    kz = est.param.kz;
    J = est.param.J;
    T = est.param.T;

    % no equality constrains
    ceq = [];

    gamma = arg(ky+1);
    rho = arg(ky+2);
    lambda = arg(ky+3);
    vect_Upsilon = arg(ky+3+p+kz+1 : ky+3+p+kz+p^2); 
    alpha_xi = arg(ky+3+p+kz+p^2+1);
    alpha =  arg(ky+3+p+kz+p^2+1+1 : ky+3+p+kz+p^2+1+J) ; 

    sigma_xi = 1 / alpha_xi;
    Upsilon = reshape(vect_Upsilon, p, p);
    Sigma_epsilon_neg_half = zeros(p,p);
    Sigma_epsilon_neg_half(tril(true(p,p))) = alpha;
    Sigma_epsilon_neg_half = Sigma_epsilon_neg_half + tril(Sigma_epsilon_neg_half,-1)';
    Sigma_epsilon = Sigma_epsilon_neg_half^(-2);
    delta = arg(ky+3+1 : ky+3+p);
    sigma_v_epsilon = Sigma_epsilon * delta;
    sigma_v = sqrt(sigma_xi^2 + sigma_v_epsilon' * delta);
    Sigma_v_epsilon = [sigma_v^2, sigma_v_epsilon'; 
                    sigma_v_epsilon, Sigma_epsilon];
    
    c_w = -Inf;
    for t = 1 : T
        if norm(W(:,:,t), Inf) > c_w
            c_w = norm(W(:,:,t), Inf);
        end
    end

    
    
    % multiple inequality constrains c(x) <= 0
    c = [ norm(Upsilon,1) - 1;  ...     % norm(Upsilon,1) < 1. Stable condition for Z
          abs(lambda)*c_w + abs(rho) + abs(gamma) * c_w - 1;  % Stable condition for y. See Qu et.al (2017)
        ];






end