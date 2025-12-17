clear; close all; clc;
addpath("functions");
rng(0) % the random number generator start from the same starting point
%% System parameters
[para] = parameter(); % System parameters stored in function
Movable_region = [9,16,25,36,49,64,81,100]; % available region for antenna move
CASE = Movable_region;
XX = length(CASE);
L  = para.path_num;
K  = para.user_num;
N  = para.ant_num;
alpha = para.alpha;
%% Determination
all_rate = zeros(para.monte_carlo, XX);
mean_rate = zeros(1, XX);
%% Main
for mon = 1:para.monte_carlo
    disp(mon);
    [beta, phi, theta] = generate_DOA(para); % virtual DoAs of l path of the k user channel
    for cse = 1:XX
        G = CASE(cse);  % Size of movable region
        H_conj_trans=zeros(CASE(cse),K);
        for k = 1:K
            H_conj_trans(:,k) = dictionary_channel(para,beta(:,k),phi(:,k),theta(:,k),CASE(cse));
        end
        H = H_conj_trans';
        % GA Algorithm
        objective_func = @(positions) -1 * ga_func(positions, H, G, K, alpha, para); % find maximum so take negative of output
        nvars = N; % N movable antennas
        options = optimoptions('ga', 'Display', 'off', 'MaxGenerations', 200);
        [pos_opt, neg_sum_rate] = ga(objective_func, nvars, [], [], [], [], [], [], [], [], options);
        % Calculate F
        pos_indices = round(pos_opt); % round indices
        pos_indices = max(1, min(G, pos_indices)); % change illegal indices
        H_sel = H(:, pos_indices)';
        F_sel = H_sel/(H_sel'*H_sel+alpha*eye(K));
        % Normalize Power
        total_power = sum(abs(F_sel(:)).^2);
        F_sel = F_sel * sqrt(para.power/total_power);
        % Calculate Sum Rate
        sum_rate = calculate_sum_rate(H_sel, F_sel, K, para.sigma_2);
        all_rate(mon, cse) = sum_rate;
    end
end
%% Average rates
for i = 1:XX
    mean_rate(i) = mean(all_rate(:,i));   
end
%% Plot
plot(Movable_region, mean_rate, '-s', 'LineWidth', 2, 'MarkerSize', 6,'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'Color', 'b');
xlabel('Size of Movable Region G');
ylabel('Sum Rate [bit/s/Hz]');
grid on;
%% Functions
function sum_rate = calculate_sum_rate(H, F, K, sigma_2)
    sum_rate = 0;
    signal_power = zeros(K:1);
    noise_power = zeros(K:1);
    SINR = zeros(K:1);
    for k = 1:K
        signal_power(k) = abs(H(:,k)' * F(:,k))^2;
        noise_power(k) = 0;
        for j = 1:K
            if j ~= k
                noise_power(k) = noise_power(k) + abs(H(:,k)' * F(:,j))^2;
            end
        end
        SINR(k) = signal_power(k) / (noise_power(k) + sigma_2);
        sum_rate = sum_rate + log2(1 + SINR(k));
    end
end

function sum_rate = ga_func(positions, H, G, K, alpha, para)
    pos_indices = round(positions);
    pos_indices = max(1, min(G, pos_indices));
    H = H(:, pos_indices)';
    F = H/(H'*H+alpha*eye(K)); 
    total_power = sum(abs(F(:)).^2);
    F = F * sqrt(para.power/total_power);
    sum_rate = calculate_sum_rate(H, F, K, para.sigma_2);
    if(~unique(pos_indices)) % Remove results that has duplicates
        sum_rate = 1000;
    end
end