function [Hy, Hz]  = H(theta,Y,W,Xy,Z,Xz,Y0,Z0,W0, est)
% theta = est.nb.theta; 

p = est.param.p;
J = est.param.J;
n = est.param.n;
T = est.param.T;
ky = est.param.ky;
kzp = est.param.kzp;
kz = est.param.kz;


In = speye(n);
Ip = speye(p);


% assign parameters
beta_y = theta(1:ky);
gamma = theta(ky+1);
rho = theta(ky+2);
lambda = theta(ky+3); 
delta = theta(ky+3+1 : ky+3+p);
beta_z = theta(ky+3+p+1 : ky+3+p+kz);
vect_Upsilon = theta(ky+3+p+kz+1 : ky+3+p+kz+p^2); 
alpha_xi = theta(ky+3+p+kz+p^2+1); 
alpha =  theta(ky+3+p+kz+p^2+1+1 : ky+3+p+kz+p^2+1+J) ; 


Sigma_epsilon_neg_half = zeros(p,p);
Sigma_epsilon_neg_half(tril(true(p,p))) = alpha;
Sigma_epsilon_neg_half = Sigma_epsilon_neg_half + tril(Sigma_epsilon_neg_half,-1)';



%% Lz and Ly
last_Z = Z0;
last_y = Y0;
last_Wy = W0*Y0;
S = nan(n,n,T);

Hz = nan(n*p,T);
Hy = nan(n,T);

for t = 1 : T

    S(:,:,t) = In - lambda * W(:,:,t);

    z_nt = vect(Z(:,:,t));
    y_nt = Y(:,t);
    
    Xnty = reshape(Xy(:,t,:), n, ky);

    Xntz = [];    
    for j = 1 : p
        Xntz = blkdiag(Xntz, reshape(Xz(:,t,1:kzp(j),j), n, kzp(j) ));
    end

    dz = kron(Sigma_epsilon_neg_half, In) ...
            * ( z_nt ...
                - kron(Ip, last_Z) * vect_Upsilon ...
                - Xntz * beta_z  );


    Hz(:,t) = dz;

    
    dy = alpha_xi * ...
           ( S(:,:,t) * y_nt ...
             - gamma * last_y ...
             - rho * last_Wy  ...
             - Xnty * beta_y  ...
             - kron(delta', In) ...
                * ( z_nt ...
                    - kron(Ip, last_Z) * vect_Upsilon ...
                    - Xntz * beta_z  )     )        ;  
    

    Hy(:,t) = dy;


    last_Z = reshape(Z(:,:,t), n, p);
    last_y = y_nt;
    last_Wy = W(:,:,t) * y_nt;
end



end