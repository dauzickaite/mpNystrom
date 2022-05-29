function [U, Lambda, Y, Go] = NystromSketch(A, n, l, mvp, rndseed,mpon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns n x l orthogonal matrix U and an l x l diagonal matrix Lambda
% of the Nystrom approximation U*Lambda*U^T of a positive semidefinite 
% n x n matrix A. The matrix-vector product Y = A*G is performed in
% precision mvp. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(rndseed)
if mpon 
    G = mp(randn(n,l),64);
else
    G = randn(n,l);
end

[Go,~] = qr(G,0); 


if mvp == 'd'
    Y = A*Go;
else
    opt.format = mvp;
    chop([],opt)
    Y = zeros(n,l);
    for i = 1:n
        Y = chop(Y + chop(A(:,i)*Go(i,:)));
    end
end

if mvp == 'd'
    shiftc = eps;
elseif mvp == 's'
    shiftc = 1.2e-7;
elseif mvp == 'h'
    shiftc = 9.8e-4;
elseif strcmp(mvp, 'q43')
    shiftc = 0.25;
elseif strcmp(mvp, 'q52')
    shiftc = 0.5;
end

shift = shiftc*norm(Y);
Ys = Y + shift*Go;
B = Go'*Ys;
C = chol((B+B')*0.5);
F = Ys/C;
[U, Sigma, ~] = svd(F,0);
Lambda = max(0, Sigma^2 - shift*eye(l));