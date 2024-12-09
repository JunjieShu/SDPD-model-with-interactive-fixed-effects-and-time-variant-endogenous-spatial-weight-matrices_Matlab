function [PC, IC, ER, GR] = ESDfactors_numer_of_factors(theta,Y,W,Xy,Z,Xz,Rz_max,Ry_max)

% theta = est.bc.theta; Rz_max = 8; Ry_max = 8;

arguments
    theta 
    Y 
    W 
    Xy 
    Z 
    Xz 
    Rz_max = min(size(Z,[1,3]) - [0,1]) - 1 
    Ry_max = min(size(Y) - [0,1]) - 1 
end

%% parameters and variables
Y0 = Y(:,1);
W0 = W(:,:,1);
Z0 = Z(:,:,1);


Y(:,1) = [];
W(:,:,1) = [];
Xy(:,1,:) = [];
Z(:,:,1) = [];
Xz(:,1,:,:) = [];


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
        indices = find(isnan(Xz(:,:,:,j)) == 1);
        [~, ~, temp_c] = ind2sub(size(Xz(:,:,:,j)),indices);
        kzp(j) = max(temp_c)-1;
    end
end
kz = sum(kzp);


%%% Rz and Ry must be smaller than n, T:
if Rz_max > n
    error('Rz_max must be smaller or equal than n.')
elseif Rz_max > T
    error('Rz_max must be smaller or equal than T.')
end

if Ry_max > n
    error('Ry_max must be smaller or equal than n.')
elseif Ry_max > T
    error('Ry_max must be smaller or equal than T.')
end





%% retrive parameters
beta_y = theta(1:ky);
gamma = theta(ky+1);
rho = theta(ky+2);
lambda = theta(ky+3);
delta = theta(ky+3+1 : ky+3+p);
beta_z = theta(ky+3+p+1 : ky+3+p+kz);
vect_Upsilon = theta(ky+3+p+kz+1 : ky+3+p+kz+p^2);
alpha_xi = theta(ky+3+p+kz+p^2+1);
sigma_xi = 1 / alpha_xi;
alpha =  theta(ky+3+p+kz+p^2+1+1 : ky+3+p+kz+p^2+1+J) ;


Sigma_epsilon_neg_half = zeros(p,p);
Sigma_epsilon_neg_half(tril(true(p,p))) = alpha;
Sigma_epsilon_neg_half = Sigma_epsilon_neg_half + tril(Sigma_epsilon_neg_half,-1)';
Sigma_epsilon = Sigma_epsilon_neg_half^(-2);


In = speye(n);
Ip = speye(p);
S = In - lambda * W;


%% calculate the squared error terms Dz, Dy.
Xy_m = Xy;  % size: n * T * ky
last_Y_m = [Y0, Y(:,1:end-1)]; % size: n * T
last_WY_m = [ W0*Y0, cell2mat(arrayfun(@(t) W(:, :, t) * Y(:, t), 1:T-1, 'UniformOutput', false))]; % size: n * T
SY_m = cell2mat(arrayfun(@(t) S(:, :, t) * Y(:, t), 1:T, 'UniformOutput', false)); % size: n * T


Z_m = reshape(Z,n*p,T);
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



%%% Dz
Xzsum = zeros(n*p,T);
for k = 1 : p^2
    Xzsum = Xzsum + last_Z_m(:,:,k) * vect_Upsilon(k);
end
for k = 1 : kz
    Xzsum = Xzsum + Xz_m(:,:,k) * beta_z(k);
end
Dz = kron(Sigma_epsilon_neg_half, In) * (Z_m - Xzsum);

%%% Dy
Xysum = zeros(n,T);
Xysum = Xysum + gamma * last_Y_m + rho * last_WY_m;
for k = 1 : ky
    Xysum = Xysum + Xy_m(:,:,k) * beta_y(k);
end
Dy = 1 / sigma_xi * (SY_m - Xysum - kron(delta', In) * (Z_m - Xzsum));

%%% eigenvalues
ez = sort(eig(Dz*Dz'), 'descend') / (n*T);
ey = sort(eig(Dy*Dy'), 'descend') / (n*T);




%% 2002 Bai and Ng PC, IC estimator
PCz = nan(Rz_max + 1, 3); % we also want to cover t=0 case, so + 1
ICz = nan(Rz_max + 1, 3);

PCy = nan(Ry_max + 1, 3);
ICy = nan(Ry_max + 1, 3);

% the first column is for z, the second one is for Y
PC = nan(3,2);
IC = nan(3,2);


% we calculate Rz first
for r = 0 : Rz_max

    Vz_r = sum(ez(r+1: end));

    PCz(r+1, :) = Vz_r + r * g(n,T)';  % p. 199
    ICz(r+1, :) = log(Vz_r) + r * g(n,T)';  % p. 200

end

[~, index] = min(PCz);
PC(:,1) = index' - 1;
[~, index] = min(ICz);
IC(:,1) = index' - 1;


% then Ry
for r = 0 : Ry_max

    Vy_r = sum(ey(r+1: end));

    PCy(r + 1, :) = Vy_r + r * g(n,T)';
    ICy(r + 1, :) = log(Vy_r) + r * g(n,T)';

end

[~, index] = min(PCy);
PC(:,2) = index' - 1;
[~, index] = min(ICy);
IC(:,2) = index' - 1;




%% 2013 Ahn and Alex, see p.1207 and p.1208 for the construction of
m = min(n,T);

% we'd like to cover r = 0 case, so + 1 below
ERz = nan(Rz_max + 1, 1);
GRz = nan(Rz_max + 1, 1);
ERy = nan(Ry_max + 1, 1);
GRy = nan(Ry_max + 1, 1);
% the first column is for Z, the second one is for Y
ER = nan(1,2);
GR = nan(1,2);


% estimators

% Rz first
for r = 0 : Rz_max
    
    if r ~= 0
        % eigenvalue ratio: % p.1207
        ERz(r + 1) = ez(r) / ez(r+1);  
    
        % growth ratio: % p.1207
        mu_star_r = ez(r) / sum(ez(r+1 : end));
        mu_star_r_plus_1 = ez(r+1) / sum(ez(r+1+1 : end));
        GRz(r + 1) = log(1 + mu_star_r) / log(1 + mu_star_r_plus_1); 

    else
        
        mu_star_0 = sum(ez(1 : m)) / log(m);  % p.1208
        mu_star_1 = ez(1) / sum(ez(1+1 : end));
        
        % eigenvalue ratio
        ERz(1) = mu_star_0 / mu_star_1;

        % growth ratio
        GRz(1) = log(1 + mu_star_0) / log(1 + mu_star_1);

    end

end

[~, index] = max(ERz);
ER(1,1) =  index - 1;
[~, index] = max(GRz);
GR(1,1) =  index - 1;


% the Ry
for r = 0 : Ry_max
    
    if r ~= 0
        ERy(r + 1) = ey(r) / ey(r+1);
       
        mu_star_r = ey(r) / sum(ey(r+1 : end));
        mu_star_r_plus_1 = ey(r+1) / sum(ey(r+1+1 : end));
        GRy(r + 1) = log(1 + mu_star_r) / log(1 + mu_star_r_plus_1);
        
    else
        mu_star_0 = sum(ey(1 : m)) / log(m);
        mu_star_1 = ey(1) / sum(ey(1+1 : end));

        ERy(1) = mu_star_0 / mu_star_1;
        GRy(1) = log(1 + mu_star_0) / log(1 + mu_star_1);

    end

end

[~, index] = max(ERy);
ER(1,2) =  index - 1;
[~, index] = max(GRy);
GR(1,2) =  index - 1;









end