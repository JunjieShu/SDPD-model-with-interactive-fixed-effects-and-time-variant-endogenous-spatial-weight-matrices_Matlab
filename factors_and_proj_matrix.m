function s = factors_and_proj_matrix(s, Hz, Hy, param_struct)

%{
s = est.nb
param_struct = est.param
%}


arguments
    s struct
    Hz 
    Hy 
    param_struct struct
end


Rz = param_struct.Rz;
Ry = param_struct.Ry;
[n,T] = size(Hy);
p = size(Hz,1) / n;



%% Factors
[X,d] = eig(Hz' * Hz, 'vector');
[~, ind] = sort(d, 'descend');
X = X(:, ind);
s.Fz = X(:,1:Rz);

[X,d] = eig(Hz*Hz', 'vector');
[~, ind] = sort(d, 'descend');
X = X(:, ind);
s.Gz = X(:,1:Rz);




[X,d] = eig(Hy' * Hy, 'vector');
[~, ind] = sort(d, 'descend');
X = X(:, ind);
s.Fy = X(:,1:Ry);

[X,d] = eig(Hy*Hy', 'vector');
[~, ind] = sort(d, 'descend');
X = X(:, ind);
s.Gy = X(:,1:Ry);




s.wt_Gz = kron(s.Sigma_epsilon_neg_half, speye(n)) \ s.Gz;
s.wt_Gy_wt_Fy = s.sigma_xi * s.Gy * s.Fy' ...
    + kron(s.delta', speye(n)) * s.wt_Gz * s.Fz';


%% projection matrix
s.P_Fz = s.Fz * s.Fz';
s.M_Fz = eye(T,T) - s.P_Fz;
s.P_Gz = s.Gz * s.Gz';
s.M_Gz = eye(n*p,n*p) - s.P_Gz;

s.P_Fy = s.Fy * s.Fy';
s.M_Fy = eye(T,T) - s.P_Fy;
s.P_Gy = s.Gy * s.Gy';
s.M_Gy = eye(n,n) - s.P_Gy;



end