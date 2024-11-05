clc;
clear all;
close all;

%% Parameters
XMAX = 1500; XMIN = 0;
YMAX = 1500; YMIN = 0;
ZMAX = 450; ZMIN = 100;
N = 20; % Number of search agents
MaxIter = 100; % max iteration
numRuns = 100; % Number of runs for averaging
S_m = 1;
S_n = 1;

UAV_nums = [2, 4, 6, 8, 10]; % Different numbers of UAVs
UE_num = 100; % Number of Users

%% Storage for averaged results
AvgBestF = zeros(length(UAV_nums), 3); % To store average BestF for each UAV scenario and each algorithm

for uavIdx = 1:length(UAV_nums)
    UAV_num = UAV_nums(uavIdx);
    Dim = 3 * UAV_num; % Dimension to represent 3 coordinates for each UAV

    % Temporary storage for 100 runs
    tempBestF = zeros(numRuns, 3);

    for runIdx = 1:numRuns
        display(['run no. ', num2str(runIdx)]);
        display(['UAV Num. ', num2str(UAV_num)]);
        %% Initialize position of Users, their position is fixed 100 UEs*3
        Pos_U = [(XMAX-XMIN)*rand(UE_num,1)+XMIN, (YMAX-YMIN)*rand(UE_num,1)+YMIN, zeros(UE_num,1)];

        %% Initialize Population Positions
        PopPos = zeros(N, Dim);
        for i = 1:N
            for j = 1:3:Dim
                PopPos(i,j) = (XMAX-XMIN)*rand + XMIN;
            end
            for j = 2:3:Dim
                PopPos(i,j) = (YMAX-YMIN)*rand + YMIN;
            end
            for j = 3:3:Dim
                PopPos(i,j) = (ZMAX-ZMIN)*rand + ZMIN;
            end
        end

        % Define the search space bounds for the optimization
        xMax = [repmat(XMAX, 1, Dim/3), repmat(YMAX, 1, Dim/3), repmat(ZMAX, 1, Dim/3)];
        xMin = [repmat(XMIN, 1, Dim/3), repmat(YMIN, 1, Dim/3), repmat(ZMIN, 1, Dim/3)];

        %% Run the optimization algorithms
        [BestF, ~, ~] = HybridMGO_COA_after_new2(N, MaxIter, xMin, xMax, Dim, S_m, S_n, UAV_num, UE_num, Pos_U, PopPos);
        [BestF_coa, ~, ~] = COA2(N, MaxIter, xMin, xMax, Dim, UAV_num, UE_num, Pos_U, PopPos);
        [BestF_mgo, ~, ~] = MGO2(N, MaxIter, xMin, xMax, Dim, UAV_num, UE_num, Pos_U, PopPos);

        %% Store results for this run
        tempBestF(runIdx, :) = [BestF, BestF_coa, BestF_mgo];
    end

    %% Average the results for this UAV scenario
    AvgBestF(uavIdx, :) = mean(tempBestF, 1);
end

%% Plot the averaged results
x = UAV_nums';
y = AvgBestF;

bar(x, y);
title('Averaged Performance for Different DBS Numbers');
xlabel('Number of DBSs');
ylabel('Average Capacity (bits/s/Hz)');
legend('HCCMGA', 'COA', 'MGO');
