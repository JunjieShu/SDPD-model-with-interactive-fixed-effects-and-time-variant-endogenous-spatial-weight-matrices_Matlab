function Q = obj_Q(arg,Y,W,Xy,Z,Xz,Y0,Z0,W0,Ry,Rz)

% arg = [theta_0];
% arg = rand(size(theta_0))

p = size(Z,2);
J = p*(p+1)/2;

[n, T] = size(Y);
ky = size(Xy,3);

max_kzp = size(Xz,3);
kzp = nan(p,1);
for j = 1 : p
    if sum(sum(sum(isnan(Xz(:,:,:,j))))) == 0
        kzp(j) = max_kzp;
    else
        indices = find(isnan(Xz(:,:,:,j)));
        [~, ~, temp_c] = ind2sub(size(Xz(:,:,:,j)),indices);
        kzp(j) = max(temp_c)-1;
    end
end
kz = sum(kzp);


% assign parameters
beta_y = arg(1:ky);
gamma = tanh(arg(ky+1));
rho = tanh(arg(ky+2));
lambda = tanh(arg(ky+3)); 
delta = arg(ky+3+1 : ky+3+p);
beta_z = arg(ky+3+p+1 : ky+3+p+kz);
vect_Upsilon = arg(ky+3+p+kz+1 : ky+3+p+kz+p^2); 
alpha_xi = exp(arg(ky+3+p+kz+p^2+1)); 
sigma_xi = sqrt(1 / alpha_xi^2); 

alpha =  arg(ky+3+p+kz+p^2+1+1 : ky+3+p+kz+p^2+1+J) ; 
Sigma_epsilon_neg_half = zeros(p,p);
Sigma_epsilon_neg_half(tril(true(p,p))) = alpha;
Sigma_epsilon_neg_half = Sigma_epsilon_neg_half + tril(Sigma_epsilon_neg_half,-1)';
Sigma_epsilon = Sigma_epsilon_neg_half^(-2);

In = eye(n);
Ip = eye(p);
S = In - lambda * W;



%% write variable in matrix form
Xy_m = Xy;  % size: n * T * ky
last_Y_m = [Y0, Y(:,1:end-1)]; % size: n * T
last_WY_m = [ W0*Y0, cell2mat(arrayfun(@(t) W(:, :, t) * Y(:, t), 1:T-1, 'UniformOutput', false))]; % size: n * T
SY_m = cell2mat(arrayfun(@(t) S(:, :, t) * Y(:, t), 1:T, 'UniformOutput', false)); % size: n * T


Z_m = reshape(Z,n*p,T);
% variable of Z contains p, hard to express directly, use loop
Xz_m = nan(n*p, T, kz);
last_Z_m = nan(n*p, T, p*p);

for t = 1 : T

    if t ~= 1
        last_Z_m(:,t,:) = kron(Ip, Z(:,:,t-1) );
    else
        last_Z_m(:,t,:) = kron(Ip, Z0);
    end

    Xntz = [];    
    for j = 1 : p
        Xntz = blkdiag(Xntz, reshape(Xz(:, t, 1:kzp(j), j), n, kzp(j) ));
    end

    Xz_m(:,t,:) = Xntz;

end



%% Dz
Xzsum = zeros(n*p,T);
for k = 1 : p^2
    Xzsum = Xzsum + last_Z_m(:,:,k) * vect_Upsilon(k);
end
for k = 1 : kz
    Xzsum = Xzsum + Xz_m(:,:,k) * beta_z(k);
end
% check 1:
% Z_m - Xzsum
% wt_G_nz * F_Tz(2:end,:)' + Ec(:,2:end)
% wt_G_nz * F_Tz(2:end,:)' + reshape(permute(reshape(epsilon(n+1:end,:),n,T,p),[1,3,2]),n*p,T)

Dz = kron(Sigma_epsilon_neg_half, In) * (Z_m - Xzsum);



%% Dy
Xysum = zeros(n,T);
Xysum = Xysum + gamma * last_Y_m + rho * last_WY_m;
for k = 1 : ky
    Xysum = Xysum + Xy_m(:,:,k) * beta_y(k);
end

% check 1:
% SY_m - Xysum
% wt_G_ny * wt_F_Ty(2:end,:)' + reshape(v(n+1:end, :),n,T)

% check 2:
% reshape(v(n+1:end, :),n,T)
% reshape( xi(n+1:end, :),n,T) +  kron(delta', In) * Ec(:,2:end)

% check 3:
% SY_m - Xysum - kron(delta', In) * (Z_m - Xzsum)
% wt_G_ny * wt_F_Ty(2:end,:)' - kron(delta', In) * (wt_G_nz * F_Tz(2:end,:)') + reshape( xi(n+1:end, :),n,T) 

% check 4: 
% sigma_xi
% std(xi(n+1:end, :))


Dy = 1 / sigma_xi * (SY_m - Xysum - kron(delta', In) * (Z_m - Xzsum));
% Dy = nan(n,T);
% for t = 1 : T
%    Dy(:,t) = 1 / sigma_xi * (SY_m(:,t) - Xysum(:,t) - reshape(Z_m(:,t) - Xzsum(:,t), n, p) * delta);
% end





%% components of objective function
sum_ln_det_Snt = sum(arrayfun(@(t) log(det(S(:, :, t))), 1:T));

ez = sort(eig(Dz*Dz'), 'descend');
Lz = sum(ez(Rz+1: end)) / (n*T);

ey = sort(eig(Dy*Dy'), 'descend');
Ly = sum(ey(Ry+1: end)) / (n*T) ;




%% Q

Q = - 1/2 * log(det(Sigma_epsilon)) ...
    - 1/2 * log(sigma_xi^2)   ...
    + sum_ln_det_Snt / (n*T)   ...
    - 1/2 * Lz - 1/2 * Ly       ;

Q = - Q; 


end

