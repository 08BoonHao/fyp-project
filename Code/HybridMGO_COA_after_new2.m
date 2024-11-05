function [BestF, BestX, cnvg] = HybridMGO_COA_after_new2(N, MaxIter, LB, UB, dim, S_m, S_n, UAV_num, UE_num, Pos_U, PopPos)
    % Split population and dimensions into halves
    half_pop1 =  ceil(S_m * N / (S_m + S_n));
    half_pop2 = N - half_pop1;
    half_dim1 = ceil(S_m * dim / (S_m+S_n));
    half_dim2 = dim - half_dim1;
    % Check if LB and UB are scalars, and expand them to arrays if needed
    if numel(LB) == 1
        LB = LB * ones(1, dim);
    end
    if numel(UB) == 1
        UB = UB * ones(1, dim);
    end

    % Initialize MGO parameters (for half of the population and half of the dimensions)
    X_mgo = PopPos(1:half_pop1, 1:half_dim1);
    BestX_mgo = [];
    BestFitness_mgo = -inf;

    % Initialize COA parameters (for the other half of the population and other half of the dimensions)
    X_coa = PopPos(half_pop1+1:end, half_dim1+1:end);
    cuve_f = zeros(1, MaxIter);
    global_Cov = zeros(1, MaxIter);
    BestFitness_coa = -inf;
    best_position_coa = zeros(1, half_dim2);
    fitness_f = zeros(1, half_pop2);
    best_fitness_ff = -inf;
    BestX = [];
    % Initial fitness evaluation for MGO and COA
    Sol_Cost_mgo = zeros(half_pop1, 1); % Initialize as column vector
    for i = 1:half_pop2
        for j = 1:half_pop1
            % COA initial fitness evaluation
            fitness_ff(i,j) = ObjFunction(UAV_num, UE_num, Pos_U, [X_mgo(j, :), X_coa(i, :)]);
            if best_fitness_ff <= fitness_ff(i,j) 
                best_coa_index = i;
                best_mgo_index = j;
                best_fitness_ff = fitness_ff(i,j);
            end
        end
    end
    best_position_coa1 = X_coa(best_coa_index, :);
    BestX_mgo1 = X_mgo(best_mgo_index, :);
    for i = 1:half_pop2
        fitness_f(i) = fitness_ff(i,best_mgo_index);
    end
    for i = 1:half_pop1
        Sol_Cost_mgo(i) = fitness_ff(best_coa_index,i);
    end
    BestX =[BestX_mgo, best_position_coa];
    BestX_mgo = BestX_mgo1;
    best_position_coa = best_position_coa1;
    
    global_position_coa = best_position_coa;
    global_fitness_coa = BestFitness_coa;
    cuve_f(1) = BestFitness_coa;
    cnvg = zeros(1, MaxIter);

    % Main hybrid optimization loop
    for Iter = 1:MaxIter
        %% COA Update
        C = 2 - (Iter / MaxIter);
        temp = rand * 15 + 20;
        xf = (best_position_coa + global_position_coa) / 2;
        Xfood = best_position_coa;
        Xnew_coa = zeros(half_pop2, half_dim2); % Initialize new positions
        for i = 1:half_pop2
            if temp > 30
                if rand < 0.5
                    Xnew_coa(i, :) = X_coa(i, :) + C * rand(1, half_dim2) .* (xf - X_coa(i, :));
                else
                    for j = 1:half_dim2
                        z = round(rand * (half_pop2 - 1)) + 1;
                        Xnew_coa(i, j) = X_coa(i, j) - X_coa(z, j) + xf(j);
                    end
                end
            else
                P = 3 * rand * fitness_f(i) / ObjFunction(UAV_num, UE_num, Pos_U, [zeros(1, half_dim2), Xfood]);
                if P > 2
                    Xfood = exp(-1 / P) .* Xfood;
                    for j = 1:half_dim2
                        Xnew_coa(i, j) = X_coa(i, j) + cos(2 * pi * rand) * Xfood(j) * p_obj(temp) - sin(2 * pi * rand) * Xfood(j) * p_obj(temp);
                    end
                else
                    Xnew_coa(i, :) = (X_coa(i, :) - Xfood) .* p_obj(temp) + p_obj(temp) .* rand(1, half_dim2) .* X_coa(i, :);

                end
            end
        end

        %% Boundary conditions for COA
        for i = 1:half_pop2
            for j = 1:half_dim2
                if length(UB) == 1
                    Xnew_coa(i, j) = min(UB(half_dim1+j), Xnew_coa(i, j));
                    Xnew_coa(i, j) = max(LB(half_dim1+j), Xnew_coa(i, j));
                else
                    Xnew_coa(i, j) = min(UB(half_dim1+j), Xnew_coa(i, j));
                    Xnew_coa(i, j) = max(LB(half_dim1+j), Xnew_coa(i, j));
                end
            end
        end

        global_position_coa = Xnew_coa(1, :);
        global_fitness_coa = ObjFunction(UAV_num, UE_num, Pos_U, [BestX_mgo, global_position_coa]);

        for i = 1:half_pop2
            new_fitness = ObjFunction(UAV_num, UE_num, Pos_U, [BestX_mgo, Xnew_coa(i, :)]);
            if new_fitness > global_fitness_coa
                global_fitness_coa = new_fitness;
                global_position_coa = Xnew_coa(i, :);
            end
            if new_fitness > fitness_f(i)
                fitness_f(i) = new_fitness;
                X_coa(i, :) = Xnew_coa(i, :);
                if fitness_f(i) > BestFitness_coa
                    BestFitness_coa = fitness_f(i);
                    best_position_coa = X_coa(i, :);
                end
            end
        end
        global_Cov(Iter) = global_fitness_coa;

        %% MGO Update

        for i = 1:half_pop1
            % MGO initial fitness evaluation
            Sol_Cost_mgo(i) = ObjFunction(UAV_num, UE_num, Pos_U, [X_mgo(i, :), best_position_coa]);
            if Sol_Cost_mgo(i) >= BestFitness_mgo
                BestFitness_mgo = Sol_Cost_mgo(i);
                BestX_mgo = X_mgo(i, :);
            end
        end
        for i = 1:half_pop1
            RandomSolution = randperm(half_pop1, ceil(half_pop1 / 3));
            M = X_mgo(randi([(ceil(half_pop1 / 3)), half_pop1]), :) * floor(rand) + mean(X_mgo(RandomSolution, :)) .* ceil(rand);
            cofi = Coefficient_Vector(half_dim1, Iter, MaxIter);
            A = randn(1, half_dim1) .* exp(2 - Iter * (2 / MaxIter));
            D = (abs(X_mgo(i, :)) + abs(BestX_mgo)) * (2 * rand - 1);
            NewX_mgo = Solution_Imp(X_mgo, BestX_mgo, LB(1:half_dim1), UB(1:half_dim1), half_pop1, cofi, M, A, D, i);
            [NewX_mgo, Sol_CostNew_mgo] = Boundary_Check1(NewX_mgo, @(x) ObjFunction(UAV_num, UE_num, Pos_U, x), LB(1:half_dim1), UB(1:half_dim1), best_position_coa);
            X_mgo = [X_mgo; NewX_mgo]; % Append new solution
            Sol_Cost_mgo = [Sol_Cost_mgo; Sol_CostNew_mgo]; % Append new cost
            [~, idbest] = max(Sol_Cost_mgo);
            BestX_mgo = X_mgo(idbest, :);
        end
        [Sol_Cost_mgo, SortOrder] = sort(Sol_Cost_mgo, 'descend');
        X_mgo = X_mgo(SortOrder, :);
        [BestFitness_mgo, idbest] = max(Sol_Cost_mgo);
        BestX_mgo = X_mgo(idbest, :);
        X_mgo = X_mgo(1:half_pop1, :);
        Sol_Cost_mgo = Sol_Cost_mgo(1:half_pop1, :);

        %% Hybrid Solution Exchange
        
        % Select bottom 50% solutions from each algorithm to exchange
        num_exchange_1 = ceil(S_m/(S_m+S_n) * half_pop2); %modified by yllee
        num_exchange_2 = ceil(S_n/(S_m+S_n) * half_pop1); %modified by yllee
        % Sort MGO population
        [~, idx_mgo] = sort(Sol_Cost_mgo);
        
%         % Sort COA population
        [~, idx_coa] = sort(fitness_f);
        
        % Exchange solutions
        if half_dim1 < half_dim2
            for count = 1:num_exchange_2
                X_mgo(idx_mgo(count), 1:half_dim2) = best_position_coa(1:half_dim2); %modified by yllee %best_coa_solutions(:, 1:half_dim);
            end
            for count = 1:num_exchange_1
                X_coa(idx_coa(count), 1:half_dim2) = BestX_mgo(1:half_dim2); %modified by yllee %best_mgo_solutions(:, 1:half_dim);
            end
        else
            for count = 1:num_exchange_2
                X_mgo(idx_mgo(count), 1:half_dim1) = best_position_coa(1:half_dim1); %modified by yllee %best_coa_solutions(:, 1:half_dim);
            end
            for count = 1:num_exchange_1
                X_coa(idx_coa(count), 1:half_dim1) = BestX_mgo(1:half_dim1); %modified by yllee %best_mgo_solutions(:, 1:half_dim);
            end
        end
        %amended by yllee
        for i = 1:half_pop2
            for j = 1:half_pop1
                fitness_ff(i,j) = ObjFunction(UAV_num, UE_num, Pos_U, [X_mgo(j, :), X_coa(i, :)]);
                if best_fitness_ff > fitness_ff(i,j)
                    best_position_coa1 = X_coa(i, :);
                    BestX_mgo1 = X_mgo(j, :);
                end
            end
        end
        
        best_position_coa = best_position_coa1;
        BestX_mgo = BestX_mgo1;
        
        for i = 1:half_pop2
            fitness_f(i) = ObjFunction(UAV_num, UE_num, Pos_U, [BestX_mgo, X_coa(i, :)]);
        end
        for i = 1:half_pop1
            Sol_Cost_mgo(i) = ObjFunction(UAV_num, UE_num, Pos_U, [X_mgo(i, :), best_position_coa]);
        end

        BestX = [BestX_mgo, best_position_coa];
        
        % Find the highest value between the best fitness values of COA and MGO
        BestF_COA = max(fitness_f);  % Best fitness from the COA part
        BestF_MGO = max(Sol_Cost_mgo);  % Best fitness from the MGO part

        % Take the maximum of both fitness values
        BestF = max(BestF_COA, BestF_MGO);  % This ensures you get the highest fitness value
        
        cnvg(Iter) = BestF;
    end
end

function y = p_obj(x)
    y = 0.2 * (1 / (sqrt(2 * pi) * 3)) * exp(-(x - 25).^2 / (2 * 3.^2));
end
