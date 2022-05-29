clear
% path to chop
addpath '/Users/id1917/Documents/MATLAB/chop-master' 
% path to advanpix
addpath '/Users/id1917/Documents/MATLAB/AdvanpixMCT-4.8.5.14607' 

mp.Digits(64);

example = 'poldecay'; % choose the problem
n = 1e2; % size of A
R = 10;  % effective rank A (n-R eigenvalues decay rapidly)

k_loop = 1:R;
rnd_runs = 10;
loop_ind = length(k_loop);

singleOn = true;

normind = 0;

normAlist = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,...
    1e13,1e14,1e15,1e16];

norm_count = length(normAlist);

nystr_error_mean = zeros(loop_ind,norm_count);
nystr_error_std = zeros(loop_ind,norm_count);

error_double_mean = zeros(loop_ind,norm_count);
error_double_std = zeros(loop_ind,norm_count);
error_tot_double_mean = zeros(loop_ind,norm_count);
error_tot_double_std = zeros(loop_ind,norm_count);

error_single_mean = zeros(loop_ind,norm_count);
error_single_std = zeros(loop_ind,norm_count);
error_tot_single_mean = zeros(loop_ind,norm_count);
error_tot_single_std = zeros(loop_ind,norm_count);

error_half_mean = zeros(loop_ind,norm_count);
error_half_std = zeros(loop_ind,norm_count);
error_tot_half_mean = zeros(loop_ind,norm_count);
error_tot_half_std = zeros(loop_ind,norm_count);



for normA = normAlist
    normind = normind+1;
    fprintf('norm(A) = %d \n',normA)

    % generate matrix that has R large eigenvalues equal to normA and the
    % rest are decaying

    switch example
        case 'expdecay'
            row = 1:n;
            col = 1:n;
            val = zeros(1,n);

            val(1:R) = normA; 
            q = 0.25; % rate of exp. decay: 0.1 slow, 0.25 med, 1 fast 
            for j=R+1:n
                val(j) = 10^(-q * (j-R) );
            end
            A = sparse(row,col,val);

        case 'psdNoise'
            A = diag([normA*ones(1,R),zeros(1,n-R)]);
            rng(0)
            G = randn(n);
            ksi = 1e-4; % 1e-4; 1e-2, 1e-1
            A = A + ksi*(1/n)*(G*G');

       case 'poldecay'
             % polynomial decay
            p = 2; % rate of decay: 0.5 slow, 1 med, 2 fast 
            row = 1:n;
            col = 1:n;
            val = [normA*ones(1,R), (2:(n-R+1)).^(-p)];
            A = sparse(row,col,val);

    end


    if max(max(A)) < 6.5*1e4
        halfOn = true;
        halfind = normind;
    else
        halfOn = false;
    end


    %% compute the approximations and errors

    nystr_error = zeros(rnd_runs,loop_ind);

    error_double = zeros(rnd_runs,loop_ind);
    error_single = zeros(rnd_runs,loop_ind);
    error_half = zeros(rnd_runs,loop_ind);

    error_tot_double = zeros(rnd_runs,loop_ind);
    error_tot_single = zeros(rnd_runs,loop_ind);
    error_tot_half = zeros(rnd_runs,loop_ind);

    ind = 0;

    for k=k_loop 
        ind = ind+1;
        l=0;
        

        rnd_ind = 0;
        for rndseed = 1:rnd_runs
            rnd_ind = rnd_ind+1;

            fprintf('rank of approximation: %d \n, rnd run: %d \n',k,rndseed)    
            fprintf('hig precision: 64 digits \n')
            [U, Lambda] = NystromSketch(mp(A), n, k+l, 'd',rndseed,true);
            U = U(:,1:k);
            Lambda = Lambda(1:k,1:k);  
            A_nyst = U*Lambda*U';

            nystr_error(rndseed,ind) = norm(mp(A) - mp(A_nyst));

            % double precision AQ
            fprintf('double precis. approx \n')
            [Ud, Lambdad, AQd] = NystromSketch(A, n, k+l, 'd',rndseed, false);
            Ud = Ud(:,1:k);
            Lambdad = Lambdad(1:k,1:k);
            A_nystd = Ud*Lambdad*Ud';
            error_double(rndseed,ind) = norm(mp(A_nyst) - mp(A_nystd));
            error_tot_double(rndseed,ind) = norm(mp(A) - mp(A_nystd));


            if singleOn
            % single precision A*G
                fprintf('single precis. approx \n')
                [Us, Lambdas, AQs] = NystromSketch(A, n, k+l, 's',rndseed,false);
                Us = Us(:,1:k);
                Lambdas = Lambdas(1:k,1:k);
                A_nysts = Us*Lambdas*Us';
                error_single(rndseed,ind) = norm(mp(A_nyst) - mp(A_nysts));
                error_tot_single(rndseed,ind) = norm(mp(A) - mp(A_nysts));
            end

            if halfOn
                % half precision A*G
                fprintf('half precis. approx \n')
                [Uh, Lambdah, AQh] = NystromSketch(A, n, k+l,'h',rndseed,false);
                Uh = Uh(:,1:k);
                Lambdah = Lambdah(1:k,1:k);
                A_nysth = Uh*Lambdah*Uh';    
                error_half(rndseed,ind) = norm(mp(A_nyst) - mp(A_nysth));
                error_tot_half(rndseed,ind) = norm(mp(A) - mp(A_nysth));
            end
        end
    end

    %means
    nystr_error_mean(:,normind) = mean(nystr_error,1);
    nystr_error_std(:,normind) = std(nystr_error,1);

    error_double_mean(:,normind) = mean(error_double,1);
    error_double_std(:,normind) = std(error_double,1);

    error_tot_double_mean(:,normind) = mean(error_tot_double,1);
    error_tot_double_std(:,normind) = std(error_tot_double,1);


    if singleOn
        error_single_mean(:,normind) = mean(error_single,1);
        error_single_std(:,normind) = std(error_single,1);
        error_tot_single_mean(:,normind) = mean(error_tot_single,1);
        error_tot_single_std(:,normind) = std(error_tot_single,1);
    end

    if halfOn
        error_half_mean(:,normind) = mean(error_half,1);
        error_half_std(:,normind) = std(error_half,1);
        error_tot_half_mean(:,normind) = mean(error_tot_half,1);
        error_tot_half_std(:,normind) = std(error_tot_half,1);
    end

    
end

if halfOn
    error_half_mean = error_half_mean(:,1:halfind);
    error_half_std = error_half_std(:,1:halfind);
end
 
%% plots

nud = n*mp(2^(-53),64);
nus = n*mp(2^(-24),64);
nuh = n*mp(2^(-11),64);

boundfin_d = n^(3/2)*(nud/(1-nud))*normAlist;
boundfin_s = n^(3/2)*(nus/(1-nus))*normAlist;
boundfin_h = n^(3/2)*(nuh/(1-nuh))*normAlist(1:halfind);

% k=1:9
figure;
ind_st = 1;
ind_end = 9;

for j=ind_st:ind_end
    semilogy(0:length(nystr_error_mean(j,:))-1,nystr_error_mean(j,:),...
        '-x','Color',[0 0.4470 0.7410],'MarkerSize',25,'LineWidth',10); hold on
    semilogy(0:length(error_tot_double_mean(j,:))-1,error_tot_double_mean(j,:),...
        '-x','Color',[0.8500 0.3250 0.0980],'MarkerSize',25,'LineWidth',10); hold on
    semilogy(0:length(error_tot_single_mean(j,:))-1,error_tot_single_mean(j,:),...
        '-x','Color',[0.9290 0.6940 0.1250],'MarkerSize',25,'LineWidth',10); hold on
    semilogy(0:length(error_tot_half_mean(j,:))-1,error_tot_half_mean(j,:),...
        '-x','Color',[0.4940 0.1840 0.5560],'MarkerSize',25,'LineWidth',10); hold on
end

xlim([0,16])
xticks([0 4 9 16])
xticklabels({'10^0','10^4','10^9', '10^{16}'})
ylim([1e-20 1e20])
set(gca, 'fontsize',60);

figure;
for j=ind_st:ind_end
    semilogy(0:length(error_double_mean(j,:))-1,error_double_mean(j,:),...
        'Color',[0.8500 0.3250 0.0980],'LineWidth',10); hold on
    semilogy(0:length(error_single_mean(j,:))-1,error_single_mean(j,:),...
        'Color',[0.9290 0.6940 0.1250],'LineWidth',10); hold on
    semilogy(0:length(error_half_mean(j,:))-1,error_half_mean(j,:),...
        'Color',[0.4940 0.1840 0.5560],'LineWidth',10); hold on
end

semilogy(0:length(boundfin_d)-1,boundfin_d,':','Color',[0.8500 0.3250 0.0980],...
    'LineWidth',10); hold on
semilogy(0:length(boundfin_s)-1,boundfin_s,':','Color',[0.9290 0.6940 0.1250],...
    'LineWidth',10); hold on
semilogy(0:length(boundfin_h)-1,boundfin_h,':','Color',[0.4940 0.1840 0.5560],...
    'LineWidth',10); hold on

xlim([0,16])
xticks([0 4 9 16])
xticklabels({'10^0','10^4','10^9', '10^{16}'})
set(gca, 'fontsize',60);

% k=10
figure;
ind_st = 10;
ind_end = 10;

for j=ind_st:ind_end
    semilogy(0:length(nystr_error_mean(j,:))-1,nystr_error_mean(j,:),'-x',...
        'Color',[0 0.4470 0.7410],'MarkerSize',25,'LineWidth',10); hold on
    semilogy(0:length(error_tot_double_mean(j,:))-1,error_tot_double_mean(j,:),...
        '-x','Color',[0.8500 0.3250 0.0980],'MarkerSize',25,'LineWidth',10); hold on
    semilogy(0:length(error_tot_single_mean(j,:))-1,error_tot_single_mean(j,:),...
        '-x','Color',[0.9290 0.6940 0.1250],'MarkerSize',25,'LineWidth',10); hold on
    semilogy(0:length(error_tot_half_mean(j,:))-1,error_tot_half_mean(j,:),...
        '-x','Color',[0.4940 0.1840 0.5560],'MarkerSize',25,'LineWidth',10); hold on
end

xlim([0,16])
xticks([0 4 9 16])
xticklabels({'10^0','10^4','10^9', '10^{16}'})
ylim([1e-20 1e20])
set(gca, 'fontsize',60);

figure;
for j=ind_st:ind_end
    semilogy(0:length(error_double_mean(j,:))-1,error_double_mean(j,:),...
        'Color',[0.8500 0.3250 0.0980],'LineWidth',10); hold on
    semilogy(0:length(error_single_mean(j,:))-1,error_single_mean(j,:),...
        'Color',[0.9290 0.6940 0.1250],'LineWidth',10); hold on
    semilogy(0:length(error_half_mean(j,:))-1,error_half_mean(j,:),...
        'Color',[0.4940 0.1840 0.5560],'LineWidth',10); hold on
end

semilogy(0:length(boundfin_d)-1,boundfin_d,':','Color',[0.8500 0.3250 0.0980],...
    'LineWidth',10); hold on
semilogy(0:length(boundfin_s)-1,boundfin_s,':','Color',[0.9290 0.6940 0.1250],...
    'LineWidth',10); hold on
semilogy(0:length(boundfin_h)-1,boundfin_h,':','Color',[0.4940 0.1840 0.5560],...
    'LineWidth',10); hold on


xlim([0,16])
xticks([0 4 9 16])
xticklabels({'10^0','10^4','10^9', '10^{16}'})
set(gca, 'fontsize',60);



