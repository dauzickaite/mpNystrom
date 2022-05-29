%% choose the problem and parameters
clear
% path to chop
addpath '/Users/id1917/Documents/MATLAB/chop-master' 
% path to advanpix
addpath '/Users/id1917/Documents/MATLAB/AdvanpixMCT-4.8.5.14607' 
% path to SuiteSparse problems
addpath '/Users/id1917/Documents/MATLAB/suite_sparse_problems'

mp.Digits(64);

problem = 'Journals';  %SuiteSparse problem
k_loop = 10:20:110; % choose the approximation rank
rnd_runs = 10; % choose count of random runs

solveOn = true; % switch to solve the linear systems

% choose u_p
doubleOn = true;
singleOn = true;
halfOn = true;

q43on = false;
q52on = true;

%% load the problem and solve unpreconditioned system
isAspd = true;
load([problem,'.mat'])
A = Problem.A;
n = length(A);
normA = norm(mp(full(A)));
rng(1234)
b = rand(n,1);

eigsA = sort(eig(full(A)),'descend');
lambda_Amin = min(eigsA);

mu_list = [0,0.1,0.5,1];
mu_count = length(mu_list);
condAs = zeros(1,mu_count);
FLAG = zeros(1,mu_count);
ITER = zeros(1,mu_count); 
X = zeros(n,mu_count); 
if solveOn
   for mu_ind = 1:mu_count
   % solve unpreconditioned problem
     mu = mu_list(mu_ind);
     As = A + mu*eye(n); 
     condAs(mu_ind) = cond(mp(full(As)));
     [X(:,mu_ind),FLAG(mu_ind),~,ITER(mu_ind)] =...
        pcg_nostag(As,b,1e-6,3*n);
   end
end

%% compute the approximations and errors

loop_ind = length(k_loop);

nystr_error = zeros(rnd_runs,loop_ind);

% finite precision error
error_double = zeros(rnd_runs,loop_ind);
error_single = zeros(rnd_runs,loop_ind);
error_half = zeros(rnd_runs,loop_ind);
error_q43 = zeros(rnd_runs,loop_ind);
error_q52 = zeros(rnd_runs,loop_ind);

% total approximation error
error_tot_double = zeros(rnd_runs,loop_ind);
error_tot_single = zeros(rnd_runs,loop_ind);
error_tot_half = zeros(rnd_runs,loop_ind);
error_tot_q43 = zeros(rnd_runs,loop_ind);
error_tot_q52 = zeros(rnd_runs,loop_ind);

% smallest eigenvalues
lambdamin_exact = zeros(rnd_runs,loop_ind);
lambdamin_double = zeros(rnd_runs,loop_ind);
lambdamin_single = zeros(rnd_runs,loop_ind);
lambdamin_half = zeros(rnd_runs,loop_ind);
lambdamin_q43 = zeros(rnd_runs,loop_ind);
lambdamin_q52 = zeros(rnd_runs,loop_ind);


if solveOn
    Xexact = cell(1,loop_ind);
    Xd = cell(1,loop_ind);
    Xs = cell(1,loop_ind);
    Xh = cell(1,loop_ind);
    Xq43 = cell(1,loop_ind);
    Xq52 = cell(1,loop_ind);

    FLAGexact = zeros(rnd_runs,loop_ind,mu_count);
    FLAGd = zeros(rnd_runs,loop_ind,mu_count);
    FLAGs = zeros(rnd_runs,loop_ind,mu_count);
    FLAGh = zeros(rnd_runs,loop_ind,mu_count);
    FLAGq43 = zeros(rnd_runs,loop_ind,mu_count);
    FLAGq52 = zeros(rnd_runs,loop_ind,mu_count);

    RELRESexact = zeros(rnd_runs,loop_ind,mu_count);
    RELRESd = zeros(rnd_runs,loop_ind,mu_count);
    RELRESs = zeros(rnd_runs,loop_ind,mu_count);
    RELRESh = zeros(rnd_runs,loop_ind,mu_count);
    RELRESq43 = zeros(rnd_runs,loop_ind,mu_count);
    RELRESq52 = zeros(rnd_runs,loop_ind,mu_count);

    ITERexact = zeros(rnd_runs,loop_ind,mu_count);
    ITERd = zeros(rnd_runs,loop_ind,mu_count);
    ITERs = zeros(rnd_runs,loop_ind,mu_count);
    ITERh = zeros(rnd_runs,loop_ind,mu_count);
    ITERq43 = zeros(rnd_runs,loop_ind,mu_count);
    ITERq52 = zeros(rnd_runs,loop_ind,mu_count);
    
    condPAsd_split = zeros(rnd_runs,loop_ind,mu_count);
    condPAss_split = zeros(rnd_runs,loop_ind,mu_count);
    condPAsh_split = zeros(rnd_runs,loop_ind,mu_count);   
    condPAsexact_split = zeros(rnd_runs,loop_ind,mu_count);
    condPAsq43_split = zeros(rnd_runs,loop_ind,mu_count);
    condPAsq52_split = zeros(rnd_runs,loop_ind,mu_count);
    
end
    ind = 0;

% compute rank k Nystrom approximation using mp, double, single, half and
% quarter precisions;
% double, single, and half only affect computation of AQ, other
% computations are in double    
    
for k=k_loop 
    ind = ind+1;
    l=0;

    rnd_ind = 0;
    for rndseed = 1:rnd_runs
        rnd_ind = rnd_ind+1;

        fprintf('rank of approximation: %d \n, rnd run: %d \n',k,rndseed)    

        fprintf('hig precision: 64 digits \n')
        [U, Lambda] = NystromSketch(mp(A), n, k+l,'d',rndseed,true);
        U = U(:,1:k);
        Lambda = Lambda(1:k,1:k);  
        A_nyst = U*Lambda*U';

        nystr_error(rndseed,ind) = norm(mp(A) - mp(A_nyst));
        lambdamin_exact(rndseed,ind) = Lambda(k,k);

         if doubleOn
        % double precision A*Q
            fprintf('double precis. approx \n')
            [Ud, Lambdad] = NystromSketch(A, n, k+l,'d',rndseed, false);
            Ud = Ud(:,1:k);
            Lambdad = Lambdad(1:k,1:k);
            A_nystd = Ud*Lambdad*Ud';
            error_double(rndseed,ind) = norm(mp(A_nyst) - mp(A_nystd)); 
            error_tot_double(rndseed,ind) = norm(mp(A) - mp(A_nystd));
            lambdamin_double(rndseed,ind) = Lambdad(k,k);
          end

        if singleOn
        % single precision A*Q
            fprintf('single precis. approx \n')
            [Us, Lambdas] = NystromSketch(A, n, k+l,'s',rndseed,false);
            Us = double(Us(:,1:k));
            Lambdas = double(Lambdas(1:k,1:k));
            A_nysts = Us*Lambdas*Us';
            error_single(rndseed,ind) = norm(mp(A_nyst) - mp(A_nysts));
            error_tot_single(rndseed,ind) = norm(mp(A) - mp(A_nysts));
            lambdamin_single(rndseed,ind) = Lambdas(k,k);
         end

        if halfOn
            % half precision A*Q
            fprintf('half precis. approx \n')
            [Uh, Lambdah] = NystromSketch(A, n, k+l,'h',rndseed,false);
            Uh = double(Uh(:,1:k));
            Lambdah = double(Lambdah(1:k,1:k));
            A_nysth = Uh*Lambdah*Uh';    
            error_half(rndseed,ind) = norm(mp(A_nyst) - mp(A_nysth));            
            error_tot_half(rndseed,ind) = norm(mp(A) - mp(A_nysth));
            lambdamin_half(rndseed,ind) = Lambdah(k,k);
        end
        
        if q43on
            % q43 precision A*Q
            fprintf('q43 precis. approx \n')
            [Uq43, Lambdaq43] = NystromSketch(A, n, k+l,'q43',rndseed,false);
            Uq43 = double(Uq43(:,1:k));
            Lambdaq43 = double(Lambdaq43(1:k,1:k));
            A_nystq43 = Uq43*Lambdaq43*Uq43';    
            error_q43(rndseed,ind) = norm(mp(A_nyst) - mp(A_nystq43));            
            error_tot_q43(rndseed,ind) = norm(mp(A) - mp(A_nystq43));
            lambdamin_q43(rndseed,ind) = Lambdaq43(k,k);
        end

        if q52on
            % q52 precision A*Q
            fprintf('q52 precis. approx \n')
            [Uq52, Lambdaq52] = NystromSketch(A, n, k+l,'q52',rndseed,false);
            Uq52 = double(Uq52(:,1:k));
            Lambdaq52 = double(Lambdaq52(1:k,1:k));
            A_nystq52= Uq52*Lambdaq52*Uq52';    
            error_q52(rndseed,ind) = norm(mp(A_nyst) - mp(A_nystq52));            
            error_tot_q52(rndseed,ind) = norm(mp(A) - mp(A_nystq52));
            lambdamin_q52(rndseed,ind) = Lambdaq52(k,k);
        end
        
        if solveOn
            for mu_ind = 1:mu_count
                mu = mu_list(mu_ind);
                As = A + mu*eye(n);
                
                % extended precision case
                PAsexact = @(x) Nystrom_Pinv(x, double(U), double(Lambda), double(mu));

                % form preconditioned matrix
                PAsexact_split = Nystrom_Pinvsplit(mp(As,64)*...
                    Nystrom_Pinvsplit(mp(eye(n),64),mp(U,64), mp(Lambda,64), mp(mu,64)),...
                    mp(U,64), mp(Lambda,64), mp(mu,64));

               condPAsexact_split(rndseed,ind,mu_ind) = cond(PAsexact_split);    
    

               [Xexact{ind}(:,rndseed,mu_ind),FLAGexact(rndseed,ind,mu_ind),...
                    RELRESexact(rndseed,ind,mu_ind),ITERexact(rndseed,ind,mu_ind)] = ...
                    pcg_nostag(A+mu*eye(n),b,1e-6,3*n,PAsexact);
                
                if doubleOn

                    PAsd  = @(x) Nystrom_Pinv(x, Ud, Lambdad, mu);

                    PAsd_split = Nystrom_Pinvsplit(mp(As,64)*...
                        Nystrom_Pinvsplit(mp(eye(n),64),mp(Ud,64), mp(Lambdad,64), mp(mu,64)),...
                        mp(Ud,64), mp(Lambdad,64), mp(mu,64));

                    condPAsd_split(rndseed,ind,mu_ind) = cond(PAsd_split);

                    [Xd{ind}(:,rndseed,mu_ind),FLAGd(rndseed,ind,mu_ind),...
                        RELRESd(rndseed,ind,mu_ind),...
                        ITERd(rndseed,ind,mu_ind)] = ...
                        pcg_nostag(A+mu*eye(n),b,1e-6,3*n,PAsd);
                end


                if singleOn
                    PAss =  @(x) Nystrom_Pinv(x, Us, Lambdas, mu);

                    PAss_split = Nystrom_Pinvsplit(mp(As,64)*...
                        Nystrom_Pinvsplit(mp(eye(n),64),mp(Us,64), mp(Lambdas,64), mp(mu,64)),...
                        mp(Us,64), mp(Lambdas,64), mp(mu,64));

                    condPAss_split(rndseed,ind,mu_ind) = cond(PAss_split);

                    [Xs{ind}(:,rndseed,mu_ind),FLAGs(rndseed,ind,mu_ind),...
                        RELRESs(rndseed,ind,mu_ind),...
                        ITERs(rndseed,ind,mu_ind)] = ...
                        pcg_nostag(A+mu*eye(n),b,1e-6,3*n,PAss);

                end


                if halfOn
                    PAsh =  @(x) Nystrom_Pinv(x, Uh, Lambdah, mu);
                    
                    PAsh_split = Nystrom_Pinvsplit(mp(As,64)*...
                        Nystrom_Pinvsplit(mp(eye(n),64),mp(Uh,64), mp(Lambdah,64), mp(mu,64)),...
                        mp(Uh,64), mp(Lambdah,64), mp(mu,64));

                    condPAsh_split(rndseed,ind,mu_ind) = cond(PAsh_split);

                    [Xh{ind}(:,rndseed,mu_ind),FLAGh(rndseed,ind,mu_ind),...
                        RELRESh(rndseed,ind,mu_ind),...
                         ITERh(rndseed,ind,mu_ind)] = ...
                         pcg_nostag(A+mu*eye(n),b,1e-6,3*n,PAsh);

                end
                
                if q43on
                    PAsq43 =  @(x) Nystrom_Pinv(x, Uq43, Lambdaq43, mu);
                    
                    PAsq43_split = Nystrom_Pinvsplit(mp(As,64)*...
                        Nystrom_Pinvsplit(mp(eye(n),64),mp(Uq43,64), mp(Lambdaq43,64), mp(mu,64)),...
                        mp(Uq43,64), mp(Lambdaq43,64), mp(mu,64));

                    condPAsq43_split(rndseed,ind,mu_ind) = cond(PAsq43_split);

                    [Xq43{ind}(:,rndseed,mu_ind),FLAGq43(rndseed,ind,mu_ind),...
                        RELRESq43(rndseed,ind,mu_ind),...
                         ITERq43(rndseed,ind,mu_ind)] = ...
                         pcg_nostag(A+mu*eye(n),b,1e-6,3*n,PAsq43);

                end
                
                if q52on
                    PAsq52 =  @(x) Nystrom_Pinv(x, Uq52, Lambdaq52, mu);
                    
                    PAsq52_split = Nystrom_Pinvsplit(mp(As,64)*...
                        Nystrom_Pinvsplit(mp(eye(n),64),mp(Uq52,64), mp(Lambdaq52,64), mp(mu,64)),...
                        mp(Uq52,64), mp(Lambdaq52,64), mp(mu,64));

                    condPAsq52_split(rndseed,ind,mu_ind) = cond(PAsq52_split);

                    [Xq52{ind}(:,rndseed,mu_ind),FLAGq52(rndseed,ind,mu_ind),...
                        RELRESq52(rndseed,ind,mu_ind),...
                         ITERq52(rndseed,ind,mu_ind)] = ...
                         pcg_nostag(A+mu*eye(n),b,1e-6,3*n,PAsq52);
                end


            end

        end

    end
end

%% means and std
nystr_error_mean = mean(nystr_error,1);
nystr_error_std = std(nystr_error,1);

if doubleOn
    error_double_mean = mean(error_double,1);    
    error_double_st = std(error_double,1);
    error_tot_double_mean = mean(error_tot_double,1);
    error_tot_double_std = std(error_tot_double,1);
end

if singleOn
    error_single_mean = mean(error_single,1);
    error_single_std = std(error_single,1);
    error_tot_single_mean = mean(error_tot_single,1);
    error_tot_single_std = std(error_tot_single,1);
end

if halfOn
    error_half_mean = mean(error_half,1);
    error_half_std = std(error_half,1);
    error_tot_half_mean = mean(error_tot_half,1);
    error_tot_half_std = std(error_tot_half,1);
end

if q43on
    error_q43_mean = mean(error_q43,1);
    error_q43_std = std(error_q43,1);
    error_tot_q43_mean = mean(error_tot_q43,1);
    error_tot_q43_std = std(error_tot_q43,1);
end

if q52on
    error_q52_mean = mean(error_q52,1);
    error_q52_std = std(error_q52,1);
    error_tot_q52_mean = mean(error_tot_q52,1);
    error_tot_q52_std = std(error_tot_q52,1);
end

%% plot normalised eigenvalues and sqrt(n)*u_p
figure; 

if singleOn; semilogy(0:n-1,sqrt(n)*2^(-24)*ones(n,1),'--','Color',...
        [0.9290 0.6940 0.1250],'LineWidth',20); hold on; end


if halfOn; semilogy(0:n-1,sqrt(n)*2^(-11)*ones(n,1),'--','Color',...
        [0.4940 0.1840 0.5560],'LineWidth',20); hold on; end

semilogy(0:n-1,eigsA./eigsA(1),'x', 'Color', [0 0.4470 0.7410], ...
    'MarkerSize',25,'Linewidth', 1); hold on

xlim([0 n-1])
xticks([0 50 100])
ylim([1e-7 1])
yticks([1e-6 1e-2 1e0])
set(gca, 'fontsize',80);

text(78,2*1e-2,' k=77','FontSize',80)
%text(320,4*1e-7,' k=320','HorizontalAlignment','right','FontSize',80)

%% plot the total error
xfin = loop_ind-1;
figure; semilogy(0:xfin,nystr_error_mean,'o','Color',[0 0.4470 0.7410],'MarkerSize',25,'Linewidth', 5);hold on
if doubleOn; semilogy(0:xfin,error_tot_double_mean,'o','Color',[0.8500 0.3250 0.0980],'MarkerSize',25,'Linewidth', 5);hold on; end
if singleOn; semilogy(0:xfin,error_tot_single_mean,'o','Color',[0.9290 0.6940 0.1250],'MarkerSize',25,'Linewidth', 5);hold on; end
if halfOn; semilogy(0:xfin,error_tot_half_mean,'o','Color',[0.4940 0.1840 0.5560],'MarkerSize',25,'LineWidth',5); hold on;end

xlabel('$k$', 'interpreter','latex')
xticks(0:xfin)
xticklabels(k_loop)
xlim([0 xfin])
set(gca, 'fontsize',60);

%% plot the finite precison error
% finite precision error bound
nud = n*mp(2^(-53),64);
nus = n*mp(2^(-24),64);
nuh = n*mp(2^(-11),64);
nuq43 = n*mp(2^(-3),64);
nuq52 = n*mp(2^(-2),64);

boundfin_d = n^(3/2)*(nud/(1-nud))*normA;
boundfin_s = n^(3/2)*(nus/(1-nus))*normA;
boundfin_h = n^(3/2)*(nuh/(1-nuh))*normA;
boundfin_q43 = n^(3/2)*(nuq43/(1-nuq43))*normA;
boundfin_q52 = n^(3/2)*(nuq52/(1-nuq52))*normA;

figure; semilogy(0:xfin,nystr_error_mean,'o','Color',[0 0.4470 0.7410],'MarkerSize',35,'Linewidth', 5);hold on
if doubleOn; semilogy(0:xfin,error_double_mean,'x','Color',[0.8500 0.3250 0.0980],'MarkerSize',35,'Linewidth', 5);hold on; end
if singleOn; semilogy(0:xfin,error_single_mean,'x','Color',[0.9290 0.6940 0.1250],'MarkerSize',35,'Linewidth', 5);hold on; end 
if halfOn; semilogy(0:xfin,error_half_mean,'x','Color',[0.4940 0.1840 0.5560],'MarkerSize',35,'Linewidth', 5);hold on; end 

if doubleOn; semilogy(0:xfin,boundfin_d*ones(xfin+1,1),':','Color',[0.8500 0.3250 0.0980],'Linewidth', 10);hold on; end
if singleOn; semilogy(0:xfin,boundfin_s*ones(xfin+1,1),':','Color',[0.9290 0.6940 0.1250],'Linewidth', 10);hold on; end
if halfOn; semilogy(0:xfin,boundfin_h*ones(xfin+1,1),':','Color',[0.4940 0.1840 0.5560],'Linewidth', 10);hold on; end
    
xlabel('$k$', 'interpreter','latex')
xticks(0:xfin)
xticklabels(k_loop)
xlim([0 xfin])
set(gca, 'fontsize',60);

%% plot the quarter precision error

if q52on
    figure; semilogy(0:xfin,nystr_error_mean,'o','Color',[0 0.4470 0.7410],...
        'MarkerSize',25,'Linewidth', 5);hold on
    semilogy(0:xfin,error_tot_q52_mean,'o','Color',[0.6350 0.0780 0.1840],...
        'MarkerSize',25,'Linewidth', 5);hold on

    xlabel('$k$', 'interpreter','latex')
    xticks(0:xfin)
    xticklabels(k_loop)
    xlim([0 xfin])
    set(gca, 'fontsize',60);

end


%% plot the condition number and estimates

if ~solveOn
    return
end


% condition number means
condPAsexact_split_mean = mean(condPAsexact_split,1);
condPAsexact_split_std = std(condPAsexact_split,1);

if doubleOn
    condPAsd_split_mean = mean(condPAsd_split,1);
    condPAsd_split_std = std(condPAsd_split,1);
end

if singleOn
    condPAss_split_mean = mean(condPAss_split,1);
    condPAss_split_std = std(condPAss_split,1);
end

if halfOn
    condPAsh_split_mean = mean(condPAsh_split,1);
    condPAsh_split_std = std(condPAsh_split,1);
end

if q43on
    condPAsq43_split_mean = mean(condPAsq43_split,1);
    condPAsq43_split_std = std(condPAsq43_split,1);
end

if q52on
    condPAsq52_split_mean = mean(condPAsq52_split,1);
    condPAsq52_split_std = std(condPAsq52_split,1);
end

% Frangella, Tropp, and Udell (2021) bound for the expected exact approx. error
indk = 0;
errorbound = zeros(1,loop_ind);
for k = k_loop
    indk = indk +1;
    if k < 4
        errorbound(indk) = nan;
    else
        errorbound(indk) = Nystrom_boundFTU(eigsA,n,mp(k));
    end
end

% estimates of the condition number bound

cond_boundl_exact = zeros(rnd_runs,loop_ind,mu_count);
cond_boundumu_exact = zeros(rnd_runs,loop_ind,mu_count);
cond_bounduspd_exact = zeros(rnd_runs,loop_ind,mu_count);

if doubleOn
    cond_boundl_d = zeros(rnd_runs,loop_ind,mu_count);
    cond_boundumu_d = zeros(rnd_runs,loop_ind,mu_count);
    cond_bounduspd_d = zeros(rnd_runs,loop_ind,mu_count);
end

if singleOn
    cond_boundl_s = zeros(rnd_runs,loop_ind,mu_count);
    cond_boundumu_s = zeros(rnd_runs,loop_ind,mu_count);
    cond_bounduspd_s = zeros(rnd_runs,loop_ind,mu_count);
end

if halfOn
    cond_boundl_h = zeros(rnd_runs,loop_ind,mu_count);
    cond_boundumu_h = zeros(rnd_runs,loop_ind,mu_count);
    cond_bounduspd_h = zeros(rnd_runs,loop_ind,mu_count);
end

for mu_ind = 1:mu_count
     mu = mu_list(mu_ind);
     indk = 0;
     for k = k_loop
        indk = indk+1;
        rnd_ind = 0;
        for rndseed = 1:rnd_runs
            rnd_ind = rnd_ind+1;           
                        
            [bl_ex,bu_mu_ex,bu_spd_ex] = ...
                condNoBound(lambdamin_exact(rndseed,indk),errorbound(indk),...
               0,mu,lambda_Amin,isAspd);
            cond_boundl_exact(rndseed,indk,mu_ind) = bl_ex;
            cond_boundumu_exact(rndseed,indk,mu_ind) = bu_mu_ex;
            cond_bounduspd_exact(rndseed,indk,mu_ind) = bu_spd_ex;

            if doubleOn
                [bl_d,bu_mu_d,bu_spd_d] = ...
                    condNoBound(lambdamin_double(rndseed,indk),errorbound(indk),...
                   boundfin_d,mu,lambda_Amin,isAspd);
                cond_boundl_d(rndseed,indk,mu_ind) = bl_d;
                cond_boundumu_d(rndseed,indk,mu_ind) = bu_mu_d;
                cond_bounduspd_d(rndseed,indk,mu_ind) = bu_spd_d;
            end
        
            if singleOn
                [bl_s,bu_mu_s,bu_spd_s] = ...
                    condNoBound(lambdamin_single(rndseed,indk),errorbound(indk),...
                   boundfin_s,mu,lambda_Amin,isAspd);
                cond_boundl_s(rndseed,indk,mu_ind) = bl_s;
                cond_boundumu_s(rndseed,indk,mu_ind) = bu_mu_s;
                cond_bounduspd_s(rndseed,indk,mu_ind) = bu_spd_s;
            end
            
            if halfOn
                [bl_h,bu_mu_h,bu_spd_h] = ...
                    condNoBound(lambdamin_half(rndseed,indk),errorbound(indk),...
                   boundfin_h,mu,lambda_Amin,isAspd);
                cond_boundl_h(rndseed,indk,mu_ind) = bl_h;
                cond_boundumu_h(rndseed,indk,mu_ind) = bu_mu_h;
                cond_bounduspd_h(rndseed,indk,mu_ind) = bu_spd_h;
            end

        end
     end
end

% mean and std
cond_boundl_exact_mean = mean(cond_boundl_exact,1);
cond_boundumu_exact_mean = mean(cond_boundumu_exact,1);
cond_bounduspd_exact_mean = mean(cond_bounduspd_exact,1);
cond_boundl_exact_std = std(cond_boundl_exact,1);
cond_boundumu_exact_std = std(cond_boundumu_exact,1);
cond_bounduspd_exact_std = std(cond_bounduspd_exact,1);

if doubleOn
    cond_boundl_d_mean = mean(cond_boundl_d,1);
    cond_boundumu_d_mean = mean(cond_boundumu_d,1);
    cond_bounduspd_d_mean = mean(cond_bounduspd_d,1);
    cond_boundl_d_std = std(cond_boundl_d,1);
    cond_boundumu_d_std = std(cond_boundumu_d,1);
    cond_bounduspd_d_std = std(cond_bounduspd_d,1);
end

if singleOn
    cond_boundl_s_mean = mean(cond_boundl_s,1);
    cond_boundumu_s_mean = mean(cond_boundumu_s,1);
    cond_bounduspd_s_mean = mean(cond_bounduspd_s,1);
    cond_boundl_s_std = std(cond_boundl_s,1);
    cond_boundumu_s_std = std(cond_boundumu_s,1);
    cond_bounduspd_s_std = std(cond_bounduspd_s,1);

    cond_boundl_s_mean(cond_boundl_s_mean<1) = 1;
end

if halfOn
    cond_boundl_h_mean = mean(cond_boundl_h,1);
    cond_boundumu_h_mean = mean(cond_boundumu_h,1);
    cond_bounduspd_h_mean = mean(cond_bounduspd_h,1);
    cond_boundl_h_std = std(cond_boundl_h,1);
    cond_boundumu_h_std = std(cond_boundumu_h,1);
    cond_bounduspd_h_std = std(cond_bounduspd_h,1);
    
    cond_boundl_h_mean(cond_boundl_h_mean<1) = 1;
end

for mu_ind =  3 %1:mu_count
 
    figure; semilogy(0:xfin,condAs(mu_ind)*ones(1,xfin+1),'k',...
        'MarkerSize',25,'Linewidth', 10); hold on
    semilogy(0:xfin,condPAsexact_split_mean(:,:,mu_ind),'d','Color',...
        [0 0.4470 0.7410],'MarkerSize',25,'Linewidth', 5);hold on
    if doubleOn; semilogy(0:xfin,condPAsd_split_mean(:,:,mu_ind),'o','Color',...
            [0.8500 0.3250 0.0980],'MarkerSize',25,'Linewidth', 5);hold on; end
    if singleOn; semilogy(0:xfin,condPAss_split_mean(:,:,mu_ind),'+',...
            'Color',[0.9290 0.6940 0.1250],'MarkerSize',25,'Linewidth', 5);hold on; end 
    if halfOn; semilogy(0:xfin,condPAsh_split_mean(:,:,mu_ind),'x','Color',...
            [0.4940 0.1840 0.5560],'MarkerSize',25,'Linewidth', 5);hold on; end 

    semilogy(0:xfin,cond_bounduspd_exact_mean(:,:,mu_ind),...
        'b:','Linewidth', 10);hold on 
    if doubleOn; semilogy(0:xfin,cond_bounduspd_d_mean(:,:,mu_ind),':',...
        'Color',[0.8500 0.3250 0.0980],'Linewidth', 10);hold on; end
    if singleOn; semilogy(0:xfin,cond_bounduspd_s_mean(:,:,mu_ind),...
            ':','Color',[0.9290 0.6940 0.1250],'Linewidth', 10);hold on; end 
    if halfOn; semilogy(0:xfin,cond_bounduspd_h_mean(:,:,mu_ind),...
            ':','Color',[0.4940 0.1840 0.5560],'Linewidth', 10);hold on; end 
    
    semilogy(0:xfin,cond_boundl_exact_mean(:,:,mu_ind),'b',...
        'Linewidth', 10);hold on
    if doubleOn; semilogy(0:xfin,cond_boundl_d_mean(:,:,mu_ind),'Color',...
        [0.8500 0.3250 0.0980],'Linewidth', 10);hold on; end
    if singleOn; semilogy(0:xfin,cond_boundl_s_mean(:,:,mu_ind),'Color',...
            [0.9290 0.6940 0.1250],'Linewidth', 10);hold on; end
    if halfOn; semilogy(0:xfin,cond_boundl_h_mean(:,:,mu_ind),'Color',...
            [0.4940 0.1840 0.5560],'Linewidth', 10);hold on; end

    semilogy(0:xfin,cond_boundumu_exact_mean(:,:,mu_ind),'b-','Linewidth',...
        10);hold on
    if doubleOn; semilogy(0:xfin,cond_boundumu_d_mean(:,:,mu_ind),'-','Color',...
        [0.8500 0.3250 0.0980],'Linewidth', 10);hold on; end 
    if singleOn; semilogy(0:xfin,cond_boundumu_s_mean(:,:,mu_ind),'-',...
            'Color',[0.9290 0.6940 0.1250],'Linewidth', 10);hold on; end 
    if halfOn; semilogy(0:xfin,cond_boundumu_h_mean(:,:,mu_ind),'-',...
            'Color',[0.4940 0.1840 0.5560],'Linewidth', 10);hold on; end

    xlabel('$k$', 'interpreter','latex')
    xticks(0:xfin)
    xlim([0 xfin])
    xticklabels(k_loop)
    set(gca, 'fontsize',60);
end

%% plot the iteration count

ITERexact_mean = mean(ITERexact,1);
if doubleOn; ITERd_mean= mean(ITERd,1); end
if singleOn; ITERs_mean= mean(ITERs,1); end
if halfOn; ITERh_mean= mean(ITERh,1); end

ITERexact_std = std(ITERexact,1);
if doubleOn; ITERd_std= std(ITERd,1); end
if singleOn; ITERs_std= std(ITERs,1); end
if halfOn; ITERh_std= std(ITERh,1); end

% plot

for mu_ind = 3 %1:mu_count
    figure; plot(0:xfin,ITER(mu_ind)*ones(1,xfin+1),'k','MarkerSize',25,...
        'Linewidth', 10); hold on
    plot(0:xfin,ITERexact_mean(:,:,mu_ind),'d','Color',...
        [0 0.4470 0.7410],'MarkerSize',25,'Linewidth', 5); hold on
    if doubleOn; plot(0:xfin,ITERd_mean(:,:,mu_ind),'o','Color',...
            [0.8500 0.3250 0.0980],'MarkerSize',25,'Linewidth', 5); hold on; end
    if singleOn; plot(0:xfin,ITERs_mean(:,:,mu_ind),'+','Color',...
            [0.9290 0.6940 0.1250],'MarkerSize',25,'Linewidth', 5); hold on; end
    if halfOn; plot(0:xfin,ITERh_mean(:,:,mu_ind),'x','Color',...
            [0.4940 0.1840 0.5560],'MarkerSize',25,'Linewidth', 5); hold on; end

    xlabel('$k$', 'interpreter','latex')
    xticks(0:xfin)
    xlim([0 xfin])
    xticklabels(k_loop)
    ylabel('iteration count')
    set(gca, 'fontsize',60);
end

%% plot the quarter precision results
if q52on

    for mu_ind =  3 %1:mu_count

        figure; semilogy(0:xfin,condAs(mu_ind)*ones(1,xfin+1),'k',...
            'MarkerSize',25,'Linewidth', 10); hold on
        semilogy(0:xfin,condPAsq52_split_mean(:,:,mu_ind),'o','Color',...
                [0.6350 0.0780 0.1840],'MarkerSize',25,'Linewidth', 5);hold on 

        xlabel('$k$', 'interpreter','latex')
        xticks(0:xfin)
        xlim([0 xfin])
        xticklabels(k_loop)
        set(gca, 'fontsize',60);

    end


    % iteration count  
    ITERq52_mean= mean(ITERq52,1); 
    ITERq52_std= std(ITERq52,1); 

    for mu_ind = 3 % 1:mu_count
        figure; plot(0:xfin,ITER(mu_ind)*ones(1,xfin+1),'k','MarkerSize',25,...
            'Linewidth', 10); hold on
        plot(0:xfin,ITERq52_mean(:,:,mu_ind),'x','Color',...
                [0.6350 0.0780 0.1840],'MarkerSize',25,'Linewidth', 5); hold on

        xlabel('$k$', 'interpreter','latex')
        xticks(0:xfin)
        xlim([0 xfin])
        xticklabels(k_loop)
        ylabel('iteration count')
        set(gca, 'fontsize',60);
    end
end