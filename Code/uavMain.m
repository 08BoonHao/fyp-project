clc;
clear all;
close all;

%% Parameters
XMAX = 1500; XMIN = 0;
YMAX = 1500; YMIN = 0;
ZMAX = 450; ZMIN = 100;
UAV_num = 10;  % Number of UAVs
UE_num = 100; % Number of Users
Dim = 3 * UAV_num; % Dim to represent 3 coordinate for each UAV
N = 20; % Number of search agents
MaxIter = 100; % Max iterations
numRuns = 100; % Number of runs for averaging
S_m = 1;
S_n = 1;

%% Initialize the average convergence matrices
avgCnvgHCCMGA = zeros(MaxIter, 1);
avgCnvgCOA = zeros(MaxIter, 1);
avgCnvgMGO = zeros(MaxIter, 1);

%% Run the simulation multiple times to get the average convergence
for runIdx = 1:numRuns
    display(['run no. ', num2str(runIdx)]);
    % ---- Initialize position of Users, their position is fixed -----
    Pos_U = [(XMAX-XMIN)*rand(UE_num,1)+XMIN, (YMAX-YMIN)*rand(UE_num,1)+YMIN, zeros(UE_num,1)];

    %% Initial Search Agent Position
    PopPos = zeros(N,Dim);
    for i=1:N
        for j = 1:3:Dim
            PopPos(i,j) = (XMAX-XMIN)*rand+XMIN;
        end
        for j = 2:3:Dim
            PopPos(i,j) = (YMAX-YMIN)*rand+YMIN;
        end
        for j = 3:3:Dim
            PopPos(i,j) = (ZMAX-ZMIN)*rand+ZMIN;
        end
    end

    %% Boundary of Service Area
    xMax = zeros(1, Dim);
    xMin = zeros(1, Dim);
    for i = 1:3:Dim
        xMax(i) = XMAX; xMin(i) = XMIN;
    end
    for i = 2:3:Dim
        xMax(i) = YMAX; xMin(i) = YMIN;
    end
    for i = 3:3:Dim
        xMax(i) = ZMAX; xMin(i) = ZMIN;
    end

    %% Implementing HCCMGA
    [~, ~, cnvg] = HybridMGO_COA_after_new2(N, MaxIter, xMin, xMax, Dim, S_m, S_n, UAV_num, UE_num, Pos_U, PopPos);

    %% Implementing COA
    [~, ~, cnvg_coa] = COA2(N, MaxIter, xMin, xMax, Dim, UAV_num, UE_num, Pos_U, PopPos);

    %% Implementing MGO
    [~, ~, cnvg_mgo] = MGO2(N, MaxIter, xMin, xMax, Dim, UAV_num, UE_num, Pos_U, PopPos);

    %% Accumulate the convergence data
    avgCnvgHCCMGA = avgCnvgHCCMGA + cnvg';
    avgCnvgCOA = avgCnvgCOA + cnvg_coa';
    avgCnvgMGO = avgCnvgMGO + cnvg_mgo';
end

%% Average the convergence data
avgCnvgHCCMGA = avgCnvgHCCMGA / numRuns;
avgCnvgCOA = avgCnvgCOA / numRuns;
avgCnvgMGO = avgCnvgMGO / numRuns;

%% Plot the average convergence results
figure;

plot(avgCnvgHCCMGA, 'LineWidth', 1.5);
hold on;
plot(avgCnvgCOA, 'LineWidth', 1.5);
plot(avgCnvgMGO, 'LineWidth', 1.5);

xlabel('Iteration#');
ylabel('Average Capacity (bits/s/Hz)');
grid on;

legend('HCCMGA', 'COA', 'MGO');
hold off;
