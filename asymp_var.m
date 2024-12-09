function s = asymp_var(s, param_s, Y, W, Xy, Z, Xz, Y0, Z0, W0, message)




% s = est.nb;   param_s = est.param;



%% parameters
K = param_s.K;
n = param_s.n;
T = param_s.T;
ky = param_s.ky;
p = param_s.p;
J = param_s.J;
L = n*T;
kzp = param_s.kzp;
kz = param_s.kz;



%% some matrices will be needed
In = speye(n);
Ip = speye(p);
I_T = speye(T);
K_pL = com_mat(p,L);
K_np = com_mat(n,p);
M_Fy_otimes_M_Gy = kron(s.M_Fy, s.M_Gy);

S = eye(n) - s.lambda * W;
G = pagemtimes(W, pageinv(S));

clear S

%% reshape eqn Y
Y_last_m = [Y0, Y(:, 1:end-1)];
WY_last_m = [ W0*Y0, cell2mat(arrayfun(@(t) W(:, :, t) * Y(:, t), 1:T-1, 'UniformOutput', false))]; % size: n * T
WY_m = [cell2mat(arrayfun(@(t) W(:, :, t) * Y(:, t), 1:T, 'UniformOutput', false))]; % size: n * T

blocks = arrayfun(@(t) G(:, :, t), 1:T, 'UniformOutput', false);
wt_G = sparse(blkdiag(blocks{:}));

blocks = arrayfun(@(t) W(:, :, t), 1:T, 'UniformOutput', false);
W1L = sparse(blkdiag(blocks{:}));

clear G blocks

W2L = [zeros(n,n*(T-1)) ,  zeros(n,n);
       speye(n*(T-1)) ,    zeros(n*(T-1), n) ];

W3L = W2L * W1L;

WL = s.lambda * W1L + s.gamma * W2L + s.rho * W3L;
SL = speye(n*T) - WL;

clear WL

G1L = full(W1L / SL);  % we can check: trace(G1L) = trace(wt_G), see Qu et. al. (2018), p.176 (3,4,5)
G2L = full(W2L / SL);  % trace(G2L) = 0;    % G2L = full(W2L) / full(SL)
G3L = full(W3L / SL);  % trace(G3L) = 0

% X_Ly = reshape(Xy, n*T, ky);
% ell_y = [s.gamma * Y0 + s.rho * W0 * Y0; 
%                     zeros(L-n,1)                 ];

clear W1L W3L SL



%% reshape eqn Z
S2L_U = speye(L*p) - kron(s.Upsilon', W2L);
G2L_U = kron(Ip, W2L) / S2L_U;

clear W2L S2L_U   


%% bias and others
phi = nan(sum(K),1);
a_Lp = zeros(L*p, sum(K));
b_L = zeros(L, sum(K));

A_Lp = cell(1,sum(K));
B_L = cell(1,sum(K));
D_L_Lp = cell(1,sum(K));
A_Lp(:) = {0};
B_L(:) = {0};
D_L_Lp(:) = {0};


% beta_y
for k = 1 : ky
    Z_y_k = Xy(:,:,k);
    b_L(:,k) = vect(s.alpha_xi * s.M_Gy * Z_y_k * s.M_Fy);
end
phi(1:ky) = 0;
clear Z_y_k


% gamma
% wt_B_L = M_Fy_otimes_M_Gy * G2L;
% b_L(:,ky+1) = s.alpha_xi * M_Fy_otimes_M_Gy ...
%     * (G2L * (X_Ly * s.beta_y + ell_y + vect(s.wt_Gy_wt_Fy)) + [Y0; zeros(L-n,1)]);
% B_L{ky+1} = (wt_B_L + wt_B_L') / 2;
% D_L_Lp{ky+1} = s.alpha_xi * M_Fy_otimes_M_Gy ...
%     * kron(s.delta', G2L) * K_pL ...
%     * kron(I_T, K_np * kron(s.Sigma_epsilon_neg_half^(-1), In));
% phi(ky+1) = trace(wt_B_L) / sqrt(n*T);
% clear G2L wt_B_L

b_L(:,ky+1) = s.alpha_xi * vect(s.M_Gy * Y_last_m * s.M_Fy);
phi(ky+1) = trace(M_Fy_otimes_M_Gy * G2L) / sqrt(n*T);
clear G2L 



% rho
% wt_B_L = M_Fy_otimes_M_Gy * G3L;
% b_L(:,ky+2) = s.alpha_xi * M_Fy_otimes_M_Gy ...
%     * (G3L * (X_Ly * s.beta_y + ell_y + vect(s.wt_Gy_wt_Fy)) + [W0*Y0; zeros(L-n,1)]);
% B_L{ky+2} = (wt_B_L + wt_B_L') / 2;
% D_L_Lp{ky+2} = s.alpha_xi * M_Fy_otimes_M_Gy ...
%     * kron(s.delta', G3L) * K_pL ...
%     * kron(I_T, K_np * kron(s.Sigma_epsilon_neg_half^(-1), In));
% phi(ky+2) = trace(wt_B_L) / sqrt(n*T);
% clear G3L wt_B_L

b_L(:,ky+2) = s.alpha_xi * vect(s.M_Gy * WY_last_m * s.M_Fy);
phi(ky+2) = trace(M_Fy_otimes_M_Gy * G3L) / sqrt(n*T);
clear G3L 



% lambda
% wt_B_L = M_Fy_otimes_M_Gy * G1L;
% b_L(:,ky+3) = s.alpha_xi * M_Fy_otimes_M_Gy ...
%     * G1L * (X_Ly * s.beta_y + ell_y + vect(s.wt_Gy_wt_Fy));
% B_L{ky+3} = (wt_B_L + wt_B_L') / 2;
% D_L_Lp{ky+3} = s.alpha_xi * M_Fy_otimes_M_Gy ...
%     * kron(s.delta', G1L) * K_pL ...
%     * kron(I_T, K_np * kron(s.Sigma_epsilon_neg_half^(-1), In));
% phi(ky+3) = (trace(wt_B_L) - trace(wt_G)) / sqrt(n*T);
% clear G1L wt_G wt_B_L ell_y

b_L(:,ky+3) = s.alpha_xi * vect(s.M_Gy * WY_m * s.M_Fy);
phi(ky+3) = (trace(M_Fy_otimes_M_Gy * G1L) - trace(wt_G)) / sqrt(n*T);
clear G1L wt_G



% delta
for k = 1 : p
    ek = zeros(p,1);
    ek(k) = 1;
    
    b_L(:,ky+3+k) =  vect(s.alpha_xi * s.M_Gy * kron(ek', In) * s.wt_Gz * s.Fz' * s.M_Fy);
    D_L_Lp{ky+3+k} = kron(s.alpha_xi * s.M_Fy, s.M_Gy * kron(ek' * s.Sigma_epsilon_neg_half^(-1), In));
end
phi(ky+3+1 : ky+3+p) = 0;
clear ek


%  beta_z
X_z_np_kz = nan(n*p,kz,T);
for t = 1 : T
    Xntz = [];
    for j = 1 : p
        Xntz = blkdiag(Xntz, reshape(Xz(:,t,1:kzp(j),j), n, kzp(j) )  );
    end
    X_z_np_kz(:,:,t) = Xntz;
end

for k = 1 : kz
    Z_z_k = reshape(X_z_np_kz(:,k,:), n*p, T);
    a_Lp(:,ky+3+p+k) = vect(s.M_Gz * kron(s.Sigma_epsilon_neg_half, In) ...
        * Z_z_k * s.M_Fz);
    b_L(:,ky+3+p+k) = - s.alpha_xi * vect(s.M_Gy * kron(s.delta', In) ...
        * Z_z_k * s.M_Fy);
end
phi(ky+3+p+1 : ky+3+p+kz) = 0;
clear X_z_np_kz Z_z_k


%vect_Upsilon
Q = zeros(n,L,T);
for t = 1 : T
    Q( : , (t-1)*n+1 : t*n ,t) = eye(n);
end

for k = 1 : p^2
    Ek = zeros(p,p);
    Ek(k) = 1;

    temp_E_Q = [];
    for t = 1 : T
        temp_E_Q = [      temp_E_Q; 
                     kron(Ek', Q(:,:,t)) ];
    end
    

    wt_A = kron(s.M_Fz, s.M_Gz * kron(s.Sigma_epsilon_neg_half, In)) ...
        * temp_E_Q * G2L_U * K_pL * kron(I_T, K_np * kron(s.Sigma_epsilon_neg_half^(-1), In));

    phi(ky+3+p+kz+k) = trace(wt_A) / sqrt(n*T);


    Z_z_k = nan(n*p, T);
    for t = 1 : T
        if t == 1
            tempZk = kron(Ip, Z0);
        else
            tempZk = kron(Ip, Z(:, :, t-1) );
        end
    
        Z_z_k(:,t) = tempZk(:,k);
    end

    a_Lp(:,ky+3+p+kz+k) = vect(s.M_Gz * kron(s.Sigma_epsilon_neg_half, In) * Z_z_k * s.M_Fz);
    b_L(:,ky+3+p+kz+k) = - s.alpha_xi * vect(s.M_Gy * kron(s.delta', In) ...
        * Z_z_k * s.M_Fy);


end
clear Q Ek G2L_U temp_E_Q Xc_Lz wt_A


% alpha_xi
B_L{ky+3+p+kz+p^2+1} = - s.sigma_xi * M_Fy_otimes_M_Gy;
phi(ky+3+p+kz+p^2+1) = (sqrt(n/T)+sqrt(T/n)) * param_s.Ry * s.sigma_xi ;
% phi(ky+3+p+kz+p^2+1) = sqrt(n*T) * s.sigma_xi - s.sigma_xi * trace(M_Fy_otimes_M_Gy) / sqrt(n*T) ;
clear M_Fy_otimes_M_Gy


% alpha
for k = 1 : J
    ek = zeros(J,1);
    ek(k) = 1;

    Ek = zeros(p,p);
    Ek(tril(true(p))) = ek;
    Ek = Ek + tril(Ek,-1)';

    temp_S_E = s.Sigma_epsilon_neg_half \ Ek;

    wt_A = - kron(s.M_Fz, kron(temp_S_E , In) * s.M_Gz);
    A_Lp{ky+3+p+kz+p^2+1+k} = (wt_A' + wt_A) / 2;
    phi(ky+3+p+kz+p^2+1+k) = ...
        sqrt(n/T) * param_s.Rz * trace(temp_S_E) ...
        + sqrt(T/n) * trace( kron(temp_S_E, In) * s.P_Gz);
end
clear wt_A temp_S_E Ek 


%% Variance of C_nT
%%% Sigma
s.Sigma = zeros(sum(K), sum(K));
for k1 = 1 : sum(K)
    % percentDone = 100 * k1 / sum(K) ;
    percentDone = 100 * (k1 * sum(K) - (k1 * (k1 - 1)) / 2) / (sum(K)*(sum(K)+1)/2 ) ;
    if message
        fprintf('\b\b\b\b%3.0f%%', percentDone);  % Display progress as percentage
    end
    
    s.Sigma(k1, k1:end) = ...
        + cellfun(@(c) sum(D_L_Lp{k1} .* c, 'all' ), D_L_Lp(k1:end) ) ...
        + 2 * cellfun(@(c) sum(A_Lp{k1} .* c, 'all' ), A_Lp(k1:end) ) ...
        + 2 * cellfun(@(c) sum(B_L{k1} .* c, 'all' ), B_L(k1:end) );
     
end

s.Sigma = s.Sigma + triu(s.Sigma,1)' ...
                    + a_Lp' * a_Lp ...
                    + b_L' * b_L;

diagA_Lp = nan(L*p, sum(K));
diagB_L = nan(L, sum(K));
for k = 1 : sum(K)
    diagA_Lp(:,k) =  diag(A_Lp{k});
    diagB_L(:,k) =  diag(B_L{k});
end


if message
    fprintf('\b\b\b\b%3.0f%%', 100);  % Display progress as percentage
    fprintf('\n');  % in-Command Wondow waitbar ends
end
 

s.Sigma = s.Sigma / (n*T) ;
s.phi = phi;

%% variance matrix of estimators

s.Vtheta = inv(s.Sigma) / (n*T)  ;



end


