function [bl,bu_mu,bu_spd] = condNoBound(lambdal,e_exact,e_fin,mu,eigAmin,Aspd)


% upper bound
lamdbdal_mu = lambdal + mu;
eigaAmin_mu = eigAmin+mu;
if Aspd
    bu_spd = (lamdbdal_mu + e_exact + e_fin) * ...
        (1/lamdbdal_mu + (1+e_fin)/(eigaAmin_mu));
else
    bu_spd = NaN;
end

mu_e_fin = mu - e_fin;
if mu_e_fin <0
    bu_mu = NaN;
else
    bu_mu = 1 + (lambdal + e_exact + 2*e_fin)/mu_e_fin;
end

bl = (lambdal + mu_e_fin)/eigaAmin_mu;