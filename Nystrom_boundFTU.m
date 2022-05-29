function [bound_min,p_min] = Nystrom_boundFTU(lambda,n,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (real,int) <- (vct,int,int)
% Computes the error bound of Frangella, Tropp, Udell (2021), Proposition
% 2.2 for the randomised rank-k Nystrom approximation A_N of an (n x n)
% SPSD A:
% || A - A_N || <= min_p ((1 + 2(l-p)/(p-1)) lambda_{l-p+1} +
%                           + 2*e^2*l/(p^2 - 1)*(\sum_{j>l-p} lambda_j) ),
% where lambda_j is the jth largest eigenvalue of A and l>=4.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = sort(lambda,'descend');
p_min = 2;


c1 = 2*exp(2)*k;

for p = 2:(k-2)
    kmp = k-p; 
    bound = (1 + 2*kmp/(p-1))*lambda(k-p+1) + (c1/(p^2 -1))*sum(lambda(kmp+1:n));
    
    % set initial bound_min
    if p == 2
        bound_min = bound;
    end
    
    if bound < bound_min
        bound_min = bound;
        p_min = p;
    end
end
    
