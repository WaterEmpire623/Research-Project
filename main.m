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
%% Determination
SINR = zeros(K,1);
Signal = zeros(K,1);
Noise = zeros(K,1);
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
        % % ---------- PSO Algorithm ----------
        % objective_func = @(x) func(x, H, CASE(cse), K, alpha);
        % nvars = CASE(cse) * K;
        % options = optimoptions('particleswarm', 'MaxIterations', 20, 'Display', 'none');
        % [F_flattened] = particleswarm(objective_func, nvars, [], [], options);
        % F = reshape(F_flattened, CASE(cse), K);
        % ---------- GA Algorithm ----------
        objective_func = @(x) func(x, H, CASE(cse), K, alpha);
        nvars = CASE(cse) * K;
        options = optimoptions('ga', 'MaxGenerations', 200, 'Display', 'none');
        [F_flattened] = ga(objective_func, nvars, [], [], [], [], [], [], [], options);
        F = reshape(F_flattened, CASE(cse), K);
        % ---------- SINR ----------
        for k = 1:K
            Signal(k) = abs(H(:,k)'*F(:,k))^2;
            for j = 1:K
                if(j ~= k)
                    Noise(k,1) = Noise(k,1) + abs(H(:,k)'*F(:,j))^2;
                end
            end
            Noise(k,1) = Noise(k,1) + para.sigma_2;
            SINR(k,1) = Signal(k,1)/Noise(k,1);
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

function Result = func(x, H, N, K, alpha)
    F = reshape(x, N, K);
    Result = norm(eye(K) - H'*F,"fro")^2 + alpha*norm(F,"fro")^2;
end