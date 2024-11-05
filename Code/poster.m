clear all 
close all
clc

N=40;  %Number of search agents
F_name='F1';     %Name of the test function
T=500;           %Maximum number of iterations
Run_no=30;         % Number of independent runs

[lb,ub,dim,fobj]=Get_F(F_name); %Get details of the benchmark functions

[best_fun_1,best_position_1,global_Cov_1]=COA(N,T,lb,ub,dim,fobj);
[best_fun_2,best_position_2,global_Cov_2]=GWO(N,T,lb,ub,dim,fobj);
[best_fun_3,best_position_3,global_Cov_3]=WOA(N,T,lb,ub,dim,fobj);
[best_fun_4,best_position_4,global_Cov_4]=MGO(N,T,lb,ub,dim,fobj);
[best_fun_13,best_position_13,global_Cov_13]=HybridMGO_COA(N,T,lb,ub,dim,fobj);
hold on; % Add this line to hold the current plot

semilogy(global_Cov_13, 'LineWidth',1.5, 'Color', [1 0 0])
semilogy(global_Cov_1, 'LineWidth',1.5, 'Color', [1 0 0])
semilogy(global_Cov_4, 'LineWidth',1.5, 'Color', [1 0 0])
semilogy(global_Cov_2, 'LineWidth',1.5, 'Color', [0 0 0])
semilogy(global_Cov_3, 'LineWidth',1.5,'Color',[0 0 1])

xlabel('Iteration#');
ylabel('Best fitness');
legend('HCCMGA', 'COA', 'MGO', 'GWO', 'WOA');
grid on;

title(F_name);

hold off; % Release the hold on the plot

display(['The best optimal value of the objective function found by COA is : ', num2str(best_fun_2)]);
display(['The best optimal value of the objective function found by MGO is : ', num2str(best_fun_3)]);
display(['The best optimal value of the objective function found by Hybrid is : ', num2str(best_fun_13)]);