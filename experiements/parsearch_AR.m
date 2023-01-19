clear all, close all, clc
addpath('../algorithms/');
addpath('../tools/');

%% load data
load('source_signal.mat');
N=size(h_star,1);

%% hyperparameters
NUM_OF_TRIALS=100; % default=100
K=40000; % default=40000
SNR=10; % default=10 dB
GAMMA=0.8; % correlation factor of the AR process

%% parameter list
METHOD='CS-APA (r=2, rho=1)';
switch(METHOD)
    case 'ZA-LMS'
        % ZA-LMS
        par.gamma=3e-2; % weight of l1-norm, default=3e-2
        par.mu=1.5e-3; % stepsize, default=1.5e-3
    case 'l0-NLMS'
        % l0-NLMS
        par.gamma=2e-5; % weight of l0-norm, default=2e-5
        par.mu=1; % stepsize, default=1
        par.beta=5; % parameter of f_beta, default=5
    case 'TS-NLMS'
        % TS-NLMS
        par.lmd=3e-4; % weight of l1-norm, default=3e-4
        par.b2=1/par.lmd/N; % parameter of minimax concave penalty
    case 'OLBI'
        % OLBI
        par.gamma=3e0; % weight of l1-norm, default=3e0
        par.delta=2e-3; % stepsize, default=2e-3
    case 'HT-LMS'
        % HT-LMS
        par.lmd=1e-2; % weight of l1-norm, default=1e-2
        par.alpha=1.5e-3; % stepsize, default=1/N=1.5e-3
        par.rho=3000; % weight of vector update, default=3000
    case 'CS-APA (r=2, rho=0)'
        % the proposed CS-APA with r=2, rho=0 and variable stepsize
        par.lmd=5e-4; % weight of l1-norm, default=5e-4
        par.eta=8e-2; % the threshold of, default=8e-2
        par.C=1e-2;
        par.mu_max=0.99/(1+par.lmd);
        par.alpha=0.97;
    case 'CS-APA (r=2, rho=1)'
        % the proposed CS-APA with r=2, rho=1 and variable stepsize
        par.lmd=5e-4; % weight of l1-norm, default=5e-4
        par.eta=4e-2; % the threshold of, default=8e-2
        par.C=1e-2;
        par.mu_max=0.99/(1+par.lmd);
        par.alpha=0.97;
end

%% numerical experiment
hk=zeros(N,K+1);
eta_hk=zeros(1,K+1); % system mismatch
xi_hk=zeros(1,K+1); % sparseness measure
for tt=1:NUM_OF_TRIALS
    % generate uk and dk
    par.uk=generate_AR(GAMMA,K+N-1);
    par.uk=signal2mat(par.uk,N);
    par.dk=awgn(h_star.'*par.uk,SNR,'measured');
    
    switch(METHOD)
        case 'ZA-LMS'
            % ZA-LMS
            hk=ZA_LMS(par);
        case 'l0-NLMS'
            % l0-NLMS
            hk=l0_NLMS(par);
        case 'TS-NLMS'
            % TS-NLMS
            hk=TS_NLMS(par);
        case 'OLBI'
            % OLBI
            hk=OLBI(par);
        case 'HT-LMS'
            % HT-LMS
            hk=HT_LMS(par);
        case 'CS-APA (r=2, rho=0)'
            % the proposed CS-APA with r=2, rho=0 and variable stepsize
            hk=CS_APA_r2_rho0(par);
        case 'CS-APA (r=2, rho=1)'
            % the proposed CS-APA with r=2, rho=1 and variable stepsize
            hk=CS_APA_r2_rho1(par);
    end

    % compute the system mismatch
    eta_hk=eta_hk+10*log10(sum((hk-h_star).^2,1)/norm(h_star)^2);
    % compute the sparseness measure
    xi_hk=xi_hk+sum(1-exp(-1000*abs(hk)),1);

    if mod(tt,10)==0
        fprintf('num of trials=%d\n',tt);
    end
end

eta_hk=eta_hk/NUM_OF_TRIALS;
xi_hk=xi_hk/NUM_OF_TRIALS;

%% plot the results
figure;
subplot(1,2,1);
plot(eta_hk); % plot the system mismatch
subplot(1,2,2);
plot(xi_hk); % plot the sparseness measure

figure;
stem(h_star,'kx'); % plot the true signal
hold on;
stem(hk(:,end),'ro'); % plot the final estimate of the last trial

fprintf('System mismatch=%.2f',eta_hk(end)); % report the performance of the current method

%% save the results
save(['results/AR_',METHOD],'par','hk','eta_hk','xi_hk');
