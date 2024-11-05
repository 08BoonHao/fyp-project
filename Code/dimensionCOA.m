clear all 
close all
clc

N=40;  %Number of search agents
F_name='F5';     %Name of the test function
global temp_rand;
T=500;           %Maximum number of iterations
Run_no=30;         % Number of independent runs
S_m = 1;
S_m2 = 2;
S_m5 = 5;
S_n = 1;
S_n2 = 2;
S_n5 = 5;
[lb,ub,dim,fobj]=Get_F(F_name); %Get details of the benchmark functions
for run = 1:Run_no
    display(['run no. ', num2str(run)])
    temp_rand = rand;
    [best_fun(run),best_position,global_Cov_1(run,:)]=HybridMGO_COA_after_new(N,T,lb,ub,dim,fobj, S_m, S_n);
    [best_fun_2(run),best_position_2,global_Cov_2(run,:)]=HybridMGO_COA_after_new(N,T,lb,ub,dim,fobj, S_m2, S_n);
    [best_fun_3(run),best_position_3,global_Cov_3(run,:)]=HybridMGO_COA_after_new(N,T,lb,ub,dim,fobj, S_m, S_n2);
    [best_fun_4(run),best_position_4,global_Cov_4(run,:)]=HybridMGO_COA_after_new(N,T,lb,ub,dim,fobj, S_m5, S_n);
    [best_fun_5(run),best_position_5,global_Cov_5(run,:)]=HybridMGO_COA_after_new(N,T,lb,ub,dim,fobj, S_m, S_n5);
end
hold on; % Add this line to hold the current plot

semilogy(mean(global_Cov_1), 'LineWidth',1.5, 'Color', [1 0 0])
semilogy(mean(global_Cov_2), 'LineWidth',1.5, 'Color', [0.4940 0.1840 0.5560])
semilogy(mean(global_Cov_3), 'LineWidth',1.5, 'Color', [0 0 1])
semilogy(mean(global_Cov_4), 'LineWidth',1.5, 'Color', [0.4660 0.6740 0.1880])
semilogy(mean(global_Cov_5), 'LineWidth',1.5, 'Color', [0 0 0])

xlabel('Iteration#');
ylabel('Best fitness');
legend('HCCMGA(15:15)', 'HCCMGA(10:20)', 'HCCMGA(20:10)', 'HCCMGA(5:25)', 'HCCMGA(25:5)');
grid on;

title(F_name);

hold off; % Release the hold on the plot

display(['The mean value found by HCCMGA(15:15) is : ', num2str(mean(best_fun)), ' and std is ', num2str(std(best_fun))]);
display(['The mean value found by HCCMGA(10:20) is : ', num2str(mean(best_fun_2)), ' and std is ', num2str(std(best_fun_2))]);
display(['The mean value found by HCCMGA(20:10) is : ', num2str(mean(best_fun_3)), ' and std is ', num2str(std(best_fun_3))]);
display(['The mean value found by HCCMGA(5:25) is : ', num2str(mean(best_fun_4)), ' and std is ', num2str(std(best_fun_4))]);
display(['The mean value found by HCCMGA(25:5) is : ', num2str(mean(best_fun_5)), ' and std is ', num2str(std(best_fun_5))]);
