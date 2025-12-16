clear; close all; clc;
addpath("functions");
% rng(0) % the random number generator start from the same starting point
%% System parameters
[para] = parameter(); % System parameters stored in function

Movable_region = [9,16,25,36,49,64,81,100]; % available region for antenna move

CASE = Movable_region;
XX = length(CASE);
L = para.path_num;
K = para.user_num;
N = para.ant_num;
alpha = para.alpha;
P_max = para.power;
sigma_2 = para.sigma_2; 
%% Determination
SINR = zeros(K,1);
Signal_Power = zeros(K,1);
Noise_Power = zeros(K,1);
all_rate  = zeros(para.monte_carlo,XX);
mean_rate = zeros(1,XX);
%% Main
for mon = 1:para.monte_carlo
    disp(mon);
    [beta,phi,theta] = generate_DOA(para); % virtual DoAs of l path of the k user channel
    for cse = 1:XX
        H=zeros(CASE(cse),K);
        for k = 1:K
            H(:,k) = dictionary_channel(para,beta(:,k),phi(:,k),theta(:,k),CASE(cse));
        end
        % ---------- GA Algorithm ----------
        H_sel = GA_sel(H,N,CASE(cse),K,P_max,sigma_2);
        F = H_sel/(H_sel'*H_sel+alpha*eye(K));
        H_conj_trans = H_sel';
        for k=1:K
            F(:,k) = F(:,k) / norm(F(:,k)); % Normalization
        end
        % ---------- SINR ----------
        for k = 1:K
            Signal_Power(k) = abs(H(:,k)'*F(:,k))^2;
            for j = 1:K
                if(j ~= k)
                    Noise_Power(k,1) = Noise_Power(k,1) + abs(H(:,k)'*F(:,j))^2;
                end
            end
            Noise_Power(k,1) = Noise_Power(k,1) + para.sigma_2;
            SINR(k,1) = Signal_Power(k,1)/Noise_Power(k,1);
            all_rate(mon,cse) = all_rate(mon,cse) + log2(1 + SINR(k));
        end
    end
end
%% Average rate
for i = 1:XX
    mean_rate(i) = mean(all_rate(:,i));   
end
%% Plot
plot(Movable_region, mean_rate, '-s', 'LineWidth', 2, 'MarkerSize', 6,'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'Color', 'b');
xlabel('Size of Movable Region G');
ylabel('Sum Rate [bit/s/Hz]');
grid on;

function [H_sel] = GA_SEL(H,N,CASE(cse),K,P_max,sigma_2)
    
end

function GA_output = func(x, H, N, K, alpha)

end