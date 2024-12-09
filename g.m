function result = g(n,T)
%G Three penalty functions with respect to n and T.
%   See Bai and Ng (2002, p.201) 

c_nT = min(sqrt(n), sqrt(T));   % definition of c_nt, see Bai and Ng (2002, p.198, thrm 1) 

g1 = (n+T)/(n*T) * log((n*T)/(n+T));
g2 = (n+T)/(n*T) * log(c_nT^2);
g3 = 2 / c_nT^2 * log(c_nT^2);  % I multiply 2 in original g3 in their paper. The requirements are still satisfied in thrm 2 (p.199)

result = [g1; g2; g3];
end