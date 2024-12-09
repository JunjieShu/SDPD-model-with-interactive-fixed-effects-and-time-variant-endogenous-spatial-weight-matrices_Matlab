function test = ESDfactors_test(theta, Vtheta, theta0, significance_level)

%{
theta = est.bc.theta;
Vtheta = est.bc.Vtheta;
theta0 = theta_0;
%}
arguments
    theta (:,1) double
    Vtheta (:,1) double
    theta0 (:,1) double
    significance_level = 0.05
end

test.p = 2 * ( 1 - normcdf( abs(theta - theta0) ./ sqrt(Vtheta) ) );

test.h = (test.p < significance_level);

end

