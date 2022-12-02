%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to compute the 2-norm of A*G*pinv(G^T*A*G), where A is n x n 
% symmetric positive semidefinite and G is n x k random Gaussian matrix.
% Considered: two choices of A with different spectral properties.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% path to advanpix
addpath '/Users/id1917/Documents/MATLAB/AdvanpixMCT-4.8.5.14607' 

mp.Digits(64);

rng(123)
n = 1e2;
rdnno = 10; % runs with random data

[W,~] = qr(mp(rand(n))); 

% generate A
eigsA1 = sort(mp(rand(1,n)),'descend'); % eig slow
A1 = W*diag(eigsA1)*W';
A1sqrt = W*diag(eigsA1.^(1/2))*W';

R = 15; p=1; nnull = 10; % nnull: dimension of null space of A2 
eigsA2 = sort(mp([1e2*rand(1,R), (2:(n-R-nnull+1)).^(-p), zeros(1,nnull)]),...
    'descend'); % gap+fast+zeros
A2 = W*diag(eigsA2)*W';
A2sqrt = W*diag(eigsA2.^(1/2))*W';

figure; semilogy(eigsA1,'x','Color',[0.4660 0.6740 0.1880],'Linewidth', 10,'MarkerSize',15); hold on
semilogy(eigsA2,'o','Color',[0.6350 0.0780 0.1840],'Linewidth', 10,'MarkerSize',15); hold on
[hleg1, hleg2] = legend({' A_1',' A_2'},'FontSize',60);
objhl = findobj(hleg2, 'type', 'line'); 
set(objhl, 'Markersize', 50); 
set(hleg1,'position',[0.6 0.6 0.4 0.2])
set(gca, 'fontsize',60);
%legend boxoff  
%%
% compute norms
[nAGinvGTAG_A1G,minsvalGTAsq_A1G] = nGwinvA(A1,A1sqrt,W,1:n,rdnno); 
[nAQinvGTAG_A2G,minsvalGTAsq_A2G] = nGwinvA(A2,A2sqrt,W,1:n-nnull,rdnno); 

% compute bounds
condAk1 = (eigsA1(1)./eigsA1).^(1/2);
condAk2 = (eigsA2(1)./eigsA2(1:n-nnull)).^(1/2);

gapA1G = (eigsA1.^(1/2))./minsvalGTAsq_A1G;
gapA2G = (eigsA2(1:n-nnull).^(1/2))./minsvalGTAsq_A2G;

bound_A1G_gap = eigsA1(1)^(1/2)./minsvalGTAsq_A1G;
bound_A2G_gap = eigsA2(1)^(1/2)./minsvalGTAsq_A2G;

alpha = 1e-1;
bound_A1G_prob = (1/alpha)*condAk1.*((1:n).^(1/2));
bound_A2G_prob = (1/alpha)*condAk2.*((1:n-nnull).^(1/2));

%% plots

figure;

for rndi=1:rdnno
    semilogy(nAGinvGTAG_A1G(rndi,:),'Color',[0.4660 0.6740 0.1880],'Linewidth', 5); hold on
end
for rndi=1:rdnno
    semilogy(bound_A1G_gap(rndi,:),'m:','Linewidth', 3,'MarkerSize',5); hold on
end
for rndi=1:rdnno
    semilogy(gapA1G(rndi,:),'k--','Linewidth', 3); hold on
end
semilogy(bound_A1G_prob,'r--x','Linewidth', 5,'MarkerSize',15); hold on

xlabel('k')
% ylim([0 3000])
yticks([1e0 1e2 1e3])
set(gca, 'fontsize',60);

figure;
for rndi=1:rdnno
    semilogy(nAQinvGTAG_A2G(rndi,:),'Color',[0.6350 0.0780 0.1840],'Linewidth', 5); hold on
end
for rndi=1:rdnno
    semilogy(bound_A2G_gap(rndi,:),'m:','Linewidth', 3,'MarkerSize',5); hold on
end
for rndi=1:rdnno
    semilogy(gapA2G(rndi,:),'k--','Linewidth', 3); hold on
end

semilogy(bound_A2G_prob,'r--x','Linewidth', 5,'MarkerSize',15); hold on

xlabel('k')
xlim([0,n-nnull])
yticks([1e0 1e2 1e4])
set(gca, 'fontsize',60);
%%
function [nAGinvGTAG,minsvalGTAsq] = nGwinvA(A,Asqrt,W,ranks,rndruns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mtx <- [mtx,mtx,vct,int]
% Returns norm(A*G*pinv(G^T*A*G)) for a given n x n matrix A, and G is a 
% random Gaussian matrix. The experiment is repeated rndruns times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n = length(A);
nranks = length(ranks);

nAGinvGTAG = mp(zeros(rndruns,nranks)); 
minsvalGTAsq = mp(zeros(rndruns,nranks)); 

for k = ranks
    fprintf('k=%i \n',k)
            
  for rndrun = 1:rndruns
      rng(rndrun)
      G = mp(randn(n,k));
   
      minsvalGTAsq(rndrun,k) = min(svd(mp(G')*mp(Asqrt)));
      
      AG = mp(A)*G;
      Gwi = AG/(G'*AG);
      nAGinvGTAG(rndrun,k) = norm(Gwi);

  end
               
end

end

