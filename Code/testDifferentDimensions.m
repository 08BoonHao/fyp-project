clear all 
close all
clc

N = 40;             % Number of search agents
F_name = 'F6';      % Nafme of the test function
T = 500;            % Maximum number of iterations
Run_no = 100;       % Number of independent runs

[lb, ub, dim, fobj] = Get_F(F_name); % Get details of the benchmark functions

% Initialize accumulators
best_funs_13 = zeros(Run_no, 1);
global_Cov_13_all = zeros(Run_no, T);
best_funs_2 = zeros(Run_no, 1);
global_Cov_2_all = zeros(Run_no, T);
best_funs_3 = zeros(Run_no, 1);
global_Cov_3_all = zeros(Run_no, T);
best_funs_4 = zeros(Run_no, 1);
global_Cov_4_all = zeros(Run_no, T);
best_funs_5 = zeros(Run_no, 1);
global_Cov_5_all = zeros(Run_no, T);

for run = 1:Run_no
    [best_fun_13, ~, global_Cov_13] = HybridMGO_COA_after(N, T, lb, ub, dim, fobj);
    best_funs_13(run) = best_fun_13;
    global_Cov_13_all(run, :) = global_Cov_13;

    [best_fun_2, ~, global_Cov_2] = HybridMGO_COA_after_dim_2to1(N, T, lb, ub, dim, fobj);
    best_funs_2(run) = best_fun_2;
    global_Cov_2_all(run, :) = global_Cov_2;

    [best_fun_3, ~, global_Cov_3] = HybridMGO_COA_after_dim_1to2(N, T, lb, ub, dim, fobj);
    best_funs_3(run) = best_fun_3;
    global_Cov_3_all(run, :) = global_Cov_3;
    
    [best_fun_4, ~, global_Cov_4] = HybridMGO_COA_after_dim_5to1(N, T, lb, ub, dim, fobj);
    best_funs_4(run) = best_fun_4;
    global_Cov_4_all(run, :) = global_Cov_4;
    
    [best_fun_5, ~, global_Cov_5] = HybridMGO_COA_after_dim_1to5(N, T, lb, ub, dim, fobj);
    best_funs_5(run) = best_fun_5;
    global_Cov_5_all(run, :) = global_Cov_5;
end

% Compute average results
avg_global_Cov_13 = mean(global_Cov_13_all, 1);
avg_global_Cov_2 = mean(global_Cov_2_all, 1);
avg_global_Cov_3 = mean(global_Cov_3_all, 1);
avg_global_Cov_4 = mean(global_Cov_4_all, 1);
avg_global_Cov_5 = mean(global_Cov_5_all, 1);

avg_best_fun_13 = mean(best_funs_13);
avg_best_fun_2 = mean(best_funs_2);
avg_best_fun_3 = mean(best_funs_3);
avg_best_fun_4 = mean(best_funs_4);
avg_best_fun_5 = mean(best_funs_5);

std_1 = std(best_funs_13);
std_2 = std(best_funs_2);
std_3 = std(best_funs_3);
std_4 = std(best_funs_4);
std_5 = std(best_funs_5);

% Plot results
figure;
hold on;
semilogy(avg_global_Cov_13, 'LineWidth', 1.5, 'Color', [1 0 0]);
semilogy(avg_global_Cov_2, 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560]);
semilogy(avg_global_Cov_3, 'LineWidth', 1.5, 'Color', [0 0 1]);
semilogy(avg_global_Cov_4, 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880]);
semilogy(avg_global_Cov_5, 'LineWidth', 1.5, 'Color', [0 0 0]);

xlabel('Iteration#');
ylabel('Average Best Fitness');
legend('HCCMGA(15:15)', 'HCCMGA(10:20)', 'HCCMGA(20:10)', 'HCCMGA(5:25)', 'HCCMGA(25:5)');
grid on;
title(F_name);
hold off;

% Display average best fitness values
disp(['The average optimal value of the objective function found by HCCMGA(1:1) is : ', num2str(avg_best_fun_13)]);
disp(['The average optimal value of the objective function found by HCCMGA(1:2) is : ', num2str(avg_best_fun_2)]);
disp(['The average optimal value of the objective function found by HCCMGA(2:1) is : ', num2str(avg_best_fun_3)]);
disp(['The average optimal value of the objective function found by HCCMGA(1:5) is : ', num2str(avg_best_fun_4)]);
disp(['The average optimal value of the objective function found by HCCMGA(5:1) is : ', num2str(avg_best_fun_5)]);

disp(['The standard deviation found by HCCMGA(1:1) is : ', num2str(std_1)]);
disp(['The standard deviation found by HCCMGA(1:2) is : ', num2str(std_2)]);
disp(['The standard deviation found by HCCMGA(2:1) is : ', num2str(std_3)]);
disp(['The standard deviation found by HCCMGA(1:5) is : ', num2str(std_4)]);
disp(['The standard deviation found by HCCMGA(5:1) is : ', num2str(std_5)]);