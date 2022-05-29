function Px = Nystrom_Pinv(x, U, Lambda, mu, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (vct) <- (vct,mtx,mtx,real, real)
% Returns vector Px obtained by multiplying x and the inverse of Nystrom 
% preconditioner P:
% P^{-1} = (alpha + mu) * U * (Lambda + mu*Id)^{-1} * U^T + (Id - U*U^T),
% where U*Lambda*U^T is Nystrom approximation of a positive semidefinite A,
% alpha>0 is a shift parameter for the preconditioner, and we are solving
% A + mu*Id.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5, alpha = min(diag(Lambda)); end
if alpha==0; warning('Smallest estimated eigenvalue is 0'); end

UTx = U'*x;

Px = x - U*UTx + (alpha + mu)* U * (((diag(Lambda)+mu).^(-1)).*UTx);