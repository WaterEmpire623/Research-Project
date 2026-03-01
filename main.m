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
N  = para.ant_num;
alpha = para.alpha;
%% Determination
all_rate_ga = zeros(para.monte_carlo, XX);
all_rate_pso = zeros(para.monte_carlo, XX);
all_rate_abc = zeros(para.monte_carlo, XX);
mean_rate_ga = zeros(1, XX);
mean_rate_pso = zeros(1, XX);
mean_rate_abc = zeros(1, XX);
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
        objective_func = @(positions) -1 * ga_func(positions, H, G, K, alpha, para); % find maximum so take negative of output
        nvars = N; % N movable antennas
        lb = ones(1, nvars);      % Lower bound: index 1
        ub = G * ones(1, nvars);  % Upper bound: grid size G
        % GA Algorithm
        options_ga = optimoptions('ga', 'Display', 'off');
        [pos_opt_ga, neg_sum_rate_ga] = ga(objective_func, nvars, [], [], [], [], [], [], [], [], options_ga);
        % PSO Algorithm
        options_pso = optimoptions('particleswarm', 'Display', 'off');
        [pos_opt_pso, neg_sum_rate_pso] = particleswarm(objective_func, nvars, lb, ub, options_pso);
        % ABC Algorithm
        [pos_opt_abc, neg_sum_rate_abc] = ABC(objective_func, lb, ub, 50, 100, 20);
        % Calculate F
        pos_indices_ga = round(pos_opt_ga); % round indices
        pos_indices_pso = round(pos_opt_pso);
        pos_indices_abc = round(pos_opt_abc);
        pos_indices_ga = max(1, min(G, pos_indices_ga)); % change illegal indices
        pos_indices_pso = max(1, min(G, pos_indices_pso));
        pos_indices_abc = max(1, min(G, pos_indices_abc));
        H_sel_ga = H(:, pos_indices_ga)'; % H_sel is N x K matrix
        H_sel_pso = H(:, pos_indices_pso)';
        H_sel_abc = H(:, pos_indices_abc)';
        F_sel_ga = H_sel_ga/(H_sel_ga'*H_sel_ga+alpha*eye(K)); % F_sel is N x K matrix
        F_sel_pso = H_sel_pso/(H_sel_pso'*H_sel_pso+alpha*eye(K));
        F_sel_abc = H_sel_abc/(H_sel_abc'*H_sel_abc+alpha*eye(K));
        % Normalize Power
        total_power_ga = sum(abs(F_sel_ga(:)).^2);
        total_power_pso = sum(abs(F_sel_pso(:)).^2);
        total_power_abc = sum(abs(F_sel_abc(:)).^2);
        F_sel_ga = F_sel_ga * sqrt(para.power/total_power_ga);
        F_sel_pso = F_sel_pso * sqrt(para.power/total_power_pso);
        F_sel_abc = F_sel_abc * sqrt(para.power/total_power_abc);
        % Calculate Sum Rate
        sum_rate_ga = calculate_sum_rate(H_sel_ga, F_sel_ga, K, para.sigma_2);
        sum_rate_pso = calculate_sum_rate(H_sel_pso, F_sel_pso, K, para.sigma_2);
        sum_rate_abc = calculate_sum_rate(H_sel_abc, F_sel_abc, K, para.sigma_2);
        all_rate_ga(mon, cse) = sum_rate_ga;
        all_rate_pso(mon, cse) = sum_rate_pso;
        all_rate_abc(mon, cse) = sum_rate_abc;
    end
end
%% Average rates
for i = 1:XX
    mean_rate_ga(i) = mean(all_rate_ga(:,i));   
    mean_rate_pso(i) = mean(all_rate_pso(:,i));   
    mean_rate_abc(i) = mean(all_rate_abc(:,i));     
end
%% Plot
plot(Movable_region, mean_rate_ga, '-s', 'LineWidth', 2, 'MarkerSize', 6,'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'Color', 'b');
hold on;
plot(Movable_region, mean_rate_pso, '-s', 'LineWidth', 2, 'MarkerSize', 6,'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'Color', 'r');
hold on;
plot(Movable_region, mean_rate_abc, '-s', 'LineWidth', 2, 'MarkerSize', 6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'Color', 'g');
xlabel('Size of Movable Region G');
ylabel('Sum Rate [bit/s/Hz]');
legend('GA','PSO','ABC');
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
    if(numel(unique(pos_indices)) < numel(pos_indices)) % Remove results that has duplicates
        sum_rate = -1000;
    end
end