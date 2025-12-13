clear; close all; clc;
addpath("functions");
% rng(0) % the random number generator start from the same starting point
%% System parameters
[para] = parameter(); % System parameters stored in function

Movable_region = [9,16,25,36,49,64,81,100]; % available region for antenna move

CASE = Movable_region;
XX = length(CASE);
L  = para.path_num;
K  = para.user_num;
%% Determination
SINR = zeros(K);
Signal = zeros(K);
Noise = zeros(K);
all_rate  = zeros(para.monte_carlo,XX);
mean_rate = zeros(1,XX);
%% Main
for mon = 1:para.monte_carlo
    % disp(mon);
    [beta,phi,theta] = generate_DOA(para); % virtual DoAs of l path of the k user channel
    for cse = 1:XX
        % [H] = dictionary_channel(para,beta,phi,theta,CASE(cse));
        H_transpose=zeros(CASE(cse),K);
        for k = 1:K
            H_transpose(:,k) = dictionary_channel(para,beta(:,k),phi(:,k),theta(:,k),CASE(cse));
        end
        H = H_transpose';
        % ---------- PSO Algorithm ----------
        F = zeros(CASE(cse),K);
        objective_func = (norm(eye(K) - H*F,"fro"))^2 + para.alpha * (norm(F,"fro"))^2;
        nvars = CASE(cse)*K;
        F = particleswarm(objective_func, nvars);
        % ---------- SINR ----------
        for k = 1:K
            Signal(k) = (H_transpose(:,k)*F(:,k))^2;
            for j = 1:K
                if(j ~= k)
                    Noise(k) = Noise(j) + (H_transpose(:,k)*F(:,j))^2;
                end
            end
            Noise(k) = Noise(k) + para.sigma_2;
            SINR(k) = Signal/Noise;
            all_rate(mon,cse) = all_rate(mon,cse) + log2(1+SINR(k));
        end
    end
end
%% Average rate
for i = 1:XX
    mean_rate(i) = mean(all_rate(:,i));   
end
%% Plot