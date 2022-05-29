function Px = Nystrom_Pinvsplit(x, U, Lambda, mu, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (vct) <- (vct,mtx,mtx,real, real)
% Returns vector Px obtained by multiplying x and the split factor of 
% inverse of Nystrom preconditioner P:
% P^{-1/2} = (alpha + mu)^{-1/2} * U * (Lambda + mu*Id)^{-1/2} * U^T 
%                               + (Id - U*U^T),
% where U*Lambda*U^T is Nystrom approximation of a positive semidefinite A,
% alpha>0 is a shift parameter for the preconditioner, and we are solving
% A + mu*Id.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5, alpha = min(diag(Lambda)); end
if alpha==0; warning('Smallest estimated eigenvalue is 0'); end

UTx = U'*x;

Px = x - U*UTx + sqrt(alpha + mu)* U * (((diag(Lambda)+mu).^(-1/2)).*UTx);