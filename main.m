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
mean_rate_ga = zeros(1, XX);
all_rate_pso = zeros(para.monte_carlo, XX);
mean_rate_pso = zeros(1, XX);
all_rate_abc = zeros(para.monte_carlo, XX);
mean_rate_abc = zeros(1, XX);
all_rate_mfo = zeros(para.monte_carlo, XX);
mean_rate_mfo = zeros(1, XX);
all_rate_mpa = zeros(para.monte_carlo, XX);
mean_rate_mpa = zeros(1, XX);
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
        % GA Algorithm (Genetic Algorithm)
        options_ga = optimoptions('ga', 'Display', 'off');
        [pos_opt_ga, neg_sum_rate_ga] = ga(objective_func, nvars, [], [], [], [], [], [], [], [], options_ga);
        pos_indices_ga = round(pos_opt_ga); % round indices
        pos_indices_ga = max(1, min(G, pos_indices_ga)); % change illegal indices
        H_sel_ga = H(:, pos_indices_ga)'; % H_sel is N x K matrix
        F_sel_ga = H_sel_ga/(H_sel_ga'*H_sel_ga+alpha*eye(K)); % F_sel is N x K matrix
        total_power_ga = sum(abs(F_sel_ga(:)).^2); 
        F_sel_ga = F_sel_ga * sqrt(para.power/total_power_ga); % Normalize Power
        sum_rate_ga = calculate_sum_rate(H_sel_ga, F_sel_ga, K, para.sigma_2); % Calculate Sum Rate
        all_rate_ga(mon, cse) = sum_rate_ga;
        % PSO Algorithm (Particle Swarm Optimization)
        options_pso = optimoptions('particleswarm', 'Display', 'off');
        [pos_opt_pso, neg_sum_rate_pso] = particleswarm(objective_func, nvars, lb, ub, options_pso);
        pos_indices_pso = round(pos_opt_pso);
        pos_indices_pso = max(1, min(G, pos_indices_pso));
        H_sel_pso = H(:, pos_indices_pso)';
        F_sel_pso = H_sel_pso/(H_sel_pso'*H_sel_pso+alpha*eye(K));
        total_power_pso = sum(abs(F_sel_pso(:)).^2);
        F_sel_pso = F_sel_pso * sqrt(para.power/total_power_pso);
        sum_rate_pso = calculate_sum_rate(H_sel_pso, F_sel_pso, K, para.sigma_2);
        all_rate_pso(mon, cse) = sum_rate_pso;
        % ABC Algorithm (Artificial Bee Colony Algorithm)
        [pos_opt_abc, neg_sum_rate_abc] = ABC(objective_func, lb, ub, 50, 100, 20);
        pos_indices_abc = round(pos_opt_abc);
        pos_indices_abc = max(1, min(G, pos_indices_abc));
        H_sel_abc = H(:, pos_indices_abc)';
        F_sel_abc = H_sel_abc/(H_sel_abc'*H_sel_abc+alpha*eye(K));
        total_power_abc = sum(abs(F_sel_abc(:)).^2);
        F_sel_abc = F_sel_abc * sqrt(para.power/total_power_abc);
        sum_rate_abc = calculate_sum_rate(H_sel_abc, F_sel_abc, K, para.sigma_2);
        all_rate_abc(mon, cse) = sum_rate_abc;
        % MFO Algorithm (Moth Flame Optimization)
        [mfo_score, mfo_pos, ~] = MFO(50, 100, lb, ub, nvars, objective_func);
        pos_indices_mfo = round(mfo_pos);
        pos_indices_mfo = max(1, min(G, pos_indices_mfo));
        H_sel_mfo = H(:, pos_indices_mfo)';
        F_sel_mfo = H_sel_mfo/(H_sel_mfo'*H_sel_mfo + alpha*eye(K));
        total_power_mfo = sum(abs(F_sel_mfo(:)).^2);
        F_sel_mfo = F_sel_mfo * sqrt(para.power/total_power_mfo);
        sum_rate_mfo = calculate_sum_rate(H_sel_mfo, F_sel_mfo, K, para.sigma_2);
        all_rate_mfo(mon, cse) = sum_rate_mfo;
        % MPA Algorithm (Marine Predators Algorithm)
        [mpa_score, pos_opt_mpa, ~] = MPA(50, 100, lb, ub, nvars, objective_func);
        pos_indices_mpa = round(pos_opt_mpa);
        pos_indices_mpa = max(1, min(G, pos_indices_mpa));
        H_sel_mpa = H(:, pos_indices_mpa)';
        F_sel_mpa = H_sel_mpa/(H_sel_mpa'*H_sel_mpa + alpha*eye(K));
        total_power_mpa = sum(abs(F_sel_mpa(:)).^2);
        F_sel_mpa = F_sel_mpa * sqrt(para.power/total_power_mpa);
        sum_rate_mpa = calculate_sum_rate(H_sel_mpa, F_sel_mpa, K, para.sigma_2);
        all_rate_mpa(mon, cse) = sum_rate_mpa;
    end
end
%% Average rates
for i = 1:XX
    mean_rate_ga(i) = mean(all_rate_ga(:,i));
    mean_rate_pso(i) = mean(all_rate_pso(:,i));
    mean_rate_abc(i) = mean(all_rate_abc(:,i));
    mean_rate_mfo(i) = mean(all_rate_mfo(:,i));
    mean_rate_mpa(i) = mean(all_rate_mpa(:,i));
end
%% Plot
plot(Movable_region, mean_rate_ga, '-s', 'LineWidth', 2, 'MarkerSize', 6,'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'Color', 'b');
hold on;
plot(Movable_region, mean_rate_pso, '-s', 'LineWidth', 2, 'MarkerSize', 6,'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'Color', 'r');
hold on;
plot(Movable_region, mean_rate_abc, '-s', 'LineWidth', 2, 'MarkerSize', 6,'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'Color', 'g');
hold on;
plot(Movable_region, mean_rate_mfo, '-s', 'LineWidth', 2, 'MarkerSize', 6,'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y', 'Color', 'y');
hold on;
plot(Movable_region, mean_rate_mpa, '-s', 'LineWidth', 2, 'MarkerSize', 6,'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'c', 'Color', 'c');
xlabel('Size of Movable Region G');
ylabel('Sum Rate [bit/s/Hz]');
legend('GA','PSO','ABC','MFO','MPA');
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