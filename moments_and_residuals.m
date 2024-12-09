function s = moments_and_residuals(s, param_s, Y,W,Xy,Z,Xz,Y0,Z0,W0)

% s = est.nb; param_s = est.param

arguments
    s           struct
    param_s     struct
    Y           double
    W           double
    Xy          double
    Z           double
    Xz          double
    Y0          double
    Z0          double
    W0          double
end


p = param_s.p;
n = param_s.n;
T = param_s.T;
ky = param_s.ky;
kzp = param_s.kzp;

In = speye(n);
Ip = speye(p);

s.error.Ec = nan(n*p,T);
s.error.xi = nan(n*T,1);


Last_Z = Z0;
last_Y = Y0;
Last_WY = W0 * Y0;

for t = 1 : T

    Xntz = [];
    for j = 1 : p
        Xntz = blkdiag(Xntz, reshape(Xz(:,t,1:kzp(j),j), n, kzp(j) )  );
    end

    epsilon_t = reshape(Z(:,:,t), n*p, 1) ...
        - kron(Ip, Last_Z) * s.vect_Upsilon ...
        - Xntz * s.beta_z ...
        - s.wt_Gz * s.Fz(t, :)';

    s.error.Ec(: , t) = epsilon_t;

    epsilon_np = reshape(epsilon_t, n, p);

    v = Y(:,t) - s.lambda * W(:,:,t) * Y(:,t) ...
        - s.gamma * last_Y - s.rho * Last_WY ...
        - reshape(Xy(:,t,:), n, ky) * s.beta_y ...
        - s.wt_Gy_wt_Fy(:,t);


    s.error.xi((t-1)*n+1 : t*n) = v - epsilon_np * s.delta;

    Last_Z = Z(:,:,t);
    last_Y = Y(:,t);
    Last_WY = W(:,:,t) * Y(:,t);

end

s.error.vec_xi_L = s.alpha_xi * s.error.xi;
s.moments.mu_3_vec_xi_L = moment(s.error.vec_xi_L, 3);
s.moments.mu_4_vec_xi_L = moment(s.error.vec_xi_L, 4);

s.error.vec_epsilon_Lp = vect( kron(s.Sigma_epsilon_neg_half, In) * s.error.Ec);
s.moments.mu_3_vec_epsilon_Lp = moment(s.error.vec_epsilon_Lp, 3);
s.moments.mu_4_vec_epsilon_Lp = moment(s.error.vec_epsilon_Lp, 4);





end