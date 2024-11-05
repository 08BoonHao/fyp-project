clear all 
close all
clc

N=40;  %Number of search agents
F_name='F1';     %Name of the test function
global temp_rand;
T=500;           %Maximum number of iterations
Run_no=30;         % Number of independent runs
S_m = 1;
S_n = 1;
[lb,ub,dim,fobj]=Get_F(F_name); %Get details of the benchmark functions
for run = 1:Run_no
    display(['run no. ', num2str(run)])
    temp_rand = rand;
    [best_fun(run),best_position,global_Cov_1(run,:)]=COA(N,T,lb,ub,dim,fobj);
    % [best_fun_1,best_position_1,global_Cov_1]=SSA(N,T,lb,ub,dim,fobj);
    [best_fun_2(run),best_position_2,global_Cov_2(run,:)]=GWO(N,T,lb,ub,dim,fobj);
    [best_fun_3(run),best_position_3,global_Cov_3(run,:)]=WOA(N,T,lb,ub,dim,fobj);
    [best_fun_4(run),best_position_4,global_Cov_4(run,:)]=SCA(N,T,lb,ub,dim,fobj);
    % [best_fun_5,best_position_5,global_Cov_5]=ALO(N,T,lb,ub,dim,fobj);
    % [best_fun_6,best_position_6,global_Cov_6]=MFO(N,T,lb,ub,dim,fobj);
    % [best_fun_7,best_position_7,global_Cov_7]=DA(N,T,lb,ub,dim,fobj);
    % [best_fun_8,best_position_8,global_Cov_8]=MVO(N,T,lb,ub,dim,fobj);
    % [Convergence_curve,Ave,Sd]=EO(N,T,lb,ub,dim,fobj,Run_no);
    % [best_fun_9,best_position_9,global_Cov_9]=MPA(N,T,lb,ub,dim,fobj);
    % [best_fun_10,best_position_10,global_Cov_10]=AOA(N,T,lb,ub,dim,fobj);
    [best_fun_11(run),best_position_11,global_Cov_11(run,:)]=MGO(N,T,lb,ub,dim,fobj);
    % [best_fun_12,best_position_12,global_Cov_12]=HybridMGO_COA(N,T,lb,ub,dim,fobj);
    [best_fun_13(run),best_position_13,global_Cov_13(run,:)]=HybridMGO_COA_after_new(N,T,lb,ub,dim,fobj, S_m, S_n);

end
hold on; % Add this line to hold the current plot

semilogy(mean(global_Cov_1), 'LineWidth',1.5, 'Color', [1 0 0])
semilogy(mean(global_Cov_2), 'LineWidth',1.5, 'Color', [0 1 0])
semilogy(mean(global_Cov_3), 'LineWidth',1.5, 'Color', [1 0 1])
semilogy(mean(global_Cov_4), 'LineWidth',1.5, 'Color', [0 1 1])
% semilogy(global_Cov_5, 'LineWidth',3, 'Color', [1 0.5 0])
% semilogy(global_Cov_6, 'LineWidth',3, 'Color', [0.5 0.5 0.5])
% semilogy(global_Cov_7, 'LineWidth',3, 'Color', [0.5 0 0.5])
% semilogy(global_Cov_8, 'LineWidth',3, 'Color', [0.5 0.5 0])
% semilogy(Convergence_curve, 'LineWidth',3, 'Color', [0 0.5 0.5])
% semilogy(global_Cov_9, 'LineWidth',3, 'Color', [0.3 0.5 0.2])
% semilogy(global_Cov_10, 'LineWidth',3, 'Color', [0.5 0.1 0.2])
semilogy(mean(global_Cov_11), 'LineWidth',1.5, 'Color', [0 0 0])
% semilogy(global_Cov_12, 'LineWidth',3,'LineStyle', '--','Color',[0 0 1])
semilogy(mean(global_Cov_13), 'LineWidth',1.5,'LineStyle', '--','Color',[0.5 0.5 0.5])

xlabel('Iteration#');
ylabel('Best fitness');
legend('COA', 'GWO', 'WOA', 'SCA', 'MGO', 'HCCMGA');
grid on;

title(F_name);

hold off; % Release the hold on the plot

display(['The best optimal value of the objective function found by COA is : ', num2str(mean(best_fun)), ' and std is ', num2str(std(best_fun))]);
display(['The best optimal value of the objective function found by GWO is : ', num2str(mean(best_fun_2)), ' and std is ', num2str(std(best_fun_2))]);
display(['The best optimal value of the objective function found by WOA is : ', num2str(mean(best_fun_3)), ' and std is ', num2str(std(best_fun_3))]);
display(['The best optimal value of the objective function found by SCA is : ', num2str(mean(best_fun_4)), ' and std is ', num2str(std(best_fun_4))]);
display(['The best optimal value of the objective function found by MGO is : ', num2str(mean(best_fun_11)), ' and std is ', num2str(std(best_fun_11))]);
display(['The best optimal value of the objective function found by Hybrid is : ', num2str(mean(best_fun_13)), ' and std is ', num2str(std(best_fun_13))]);

