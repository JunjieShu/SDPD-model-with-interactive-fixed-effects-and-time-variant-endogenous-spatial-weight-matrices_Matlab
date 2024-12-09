function s = theta_to_params(s, param_s)
arguments
    s struct
    param_s struct
end

ky = param_s.ky;
p = param_s.p;
kz = param_s.kz;
J = p*(p+1)/2;


s.beta_y = s.theta(1:ky);
s.gamma = s.theta(ky+1);
s.rho = s.theta(ky+2);
s.lambda = s.theta(ky+3);
s.delta = s.theta(ky+3+1 : ky+3+p);
s.beta_z = s.theta(ky+3+p+1 : ky+3+p+kz);
s.vect_Upsilon = s.theta(ky+3+p+kz+1 : ky+3+p+kz+p^2);
s.Upsilon = reshape(s.vect_Upsilon, p, p);
s.alpha_xi = s.theta(ky+3+p+kz+p^2+1);
s.sigma_xi = 1/s.alpha_xi;
s.alpha = s.theta(ky+3+p+kz+p^2+1+1 : ky+3+p+kz+p^2+1+J);

s.Sigma_epsilon_neg_half = zeros(p,p);
s.Sigma_epsilon_neg_half(tril(true(p,p))) = s.alpha;
s.Sigma_epsilon_neg_half = s.Sigma_epsilon_neg_half + tril(s.Sigma_epsilon_neg_half,-1)';
s.Sigma_epsilon = (s.Sigma_epsilon_neg_half)^(-2);
s.Sigma_epsilon_neg_half = sqrtm(inv(s.Sigma_epsilon));


% The below estimator can be retrived from the endogenous assuption
s.sigma_v_epsilon = s.Sigma_epsilon * s.delta;
s.sigma_v = sqrt(s.sigma_xi^2 + s.sigma_v_epsilon' * s.delta  );
s.Sigma_v_epsilon = [s.sigma_v,           s.sigma_v_epsilon';
                     s.sigma_v_epsilon ,  s.Sigma_epsilon ];


end