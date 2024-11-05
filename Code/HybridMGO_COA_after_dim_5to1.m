function [BestF, BestX, cnvg] = HybridMGO_COA_after_dim_5to1(N, MaxIter, LB, UB, dim, fobj)
    % Split dimensions
    one_six_dim = floor(dim / 6);

    % Check if LB and UB are scalars, and expand them to arrays if needed
    if numel(LB) == 1
        LB = LB * ones(1, dim);
    end
    if numel(UB) == 1
        UB = UB * ones(1, dim);
    end

    % Initialize MGO parameters (first half of dimensions)
    lb_mgo = LB(1:one_six_dim);
    ub_mgo = UB(1:one_six_dim);
    X_mgo = initialization(N, one_six_dim, ub_mgo, lb_mgo);
    BestX_mgo = [];
    BestFitness_mgo = inf;

    % Initialize COA parameters (second half of dimensions)
    lb_coa = LB(one_six_dim+1:end);
    ub_coa = UB(one_six_dim+1:end);
    X_coa = initialization(N, dim - one_six_dim, ub_coa, lb_coa);
    cuve_f = zeros(1, MaxIter);
    global_Cov = zeros(1, MaxIter);
    BestFitness_coa = inf;
    best_position_coa = zeros(1, dim - one_six_dim);
    fitness_f = zeros(1, N);

    % Initial fitness evaluation for MGO and COA
    Sol_Cost_mgo = zeros(N, 1); % Initialize as column vector
    for i = 1:N
        % MGO initial fitness evaluation
        Sol_Cost_mgo(i) = fobj([X_mgo(i, :), zeros(1, dim - one_six_dim)]);
        if Sol_Cost_mgo(i) <= BestFitness_mgo
            BestFitness_mgo = Sol_Cost_mgo(i);
            BestX_mgo = X_mgo(i, :);
        end

        % COA initial fitness evaluation
        fitness_f(i) = fobj([zeros(1, one_six_dim), X_coa(i, :)]);
        if fitness_f(i) < BestFitness_coa
            BestFitness_coa = fitness_f(i);
            best_position_coa = X_coa(i, :);
        end
    end

    global_position_coa = best_position_coa;
    global_fitness_coa = BestFitness_coa;
    cuve_f(1) = BestFitness_coa;
    cnvg = zeros(1, MaxIter);

    % Main hybrid optimization loop
    for Iter = 1:MaxIter
        %% MGO Update
        for i = 1:N
            RandomSolution = randperm(N, ceil(N / 3));
            M = X_mgo(randi([(ceil(N / 3)), N]), :) * floor(rand) + mean(X_mgo(RandomSolution, :)) .* ceil(rand);
            cofi = Coefficient_Vector(one_six_dim, Iter, MaxIter);
            A = randn(1, one_six_dim) .* exp(2 - Iter * (2 / MaxIter));
            D = (abs(X_mgo(i, :)) + abs(BestX_mgo)) * (2 * rand - 1);
            NewX_mgo = Solution_Imp(X_mgo, BestX_mgo, lb_mgo, ub_mgo, N, cofi, M, A, D, i);
            [NewX_mgo, Sol_CostNew_mgo] = Boundary_Check(NewX_mgo, fobj, LB(1:one_six_dim), UB(1:one_six_dim));
            X_mgo = [X_mgo; NewX_mgo]; % Append new solution
            Sol_Cost_mgo = [Sol_Cost_mgo; Sol_CostNew_mgo]; % Append new cost
            [~, idbest] = min(Sol_Cost_mgo);
            BestX_mgo = X_mgo(idbest, :);
        end
        [Sol_Cost_mgo, SortOrder] = sort(Sol_Cost_mgo);
        X_mgo = X_mgo(SortOrder, :);
        [BestFitness_mgo, idbest] = min(Sol_Cost_mgo);
        BestX_mgo = X_mgo(idbest, :);
        X_mgo = X_mgo(1:N, :);
        Sol_Cost_mgo = Sol_Cost_mgo(1:N, :);

        %% COA Update
        C = 2 - (Iter / MaxIter);
        temp = rand * 15 + 20;
        xf = (best_position_coa + global_position_coa) / 2;
        Xfood = best_position_coa;
        Xnew_coa = zeros(N, dim - one_six_dim); % Initialize new positions
        for i = 1:N
            if temp > 30
                if rand < 0.5
                    Xnew_coa(i, :) = X_coa(i, :) + C * rand(1, dim - one_six_dim) .* (xf - X_coa(i, :));
                else
                    for j = 1:dim - one_six_dim
                        z = round(rand * (N - 1)) + 1;
                        Xnew_coa(i, j) = X_coa(i, j) - X_coa(z, j) + xf(j);
                    end
                end
            else
                P = 3 * rand * fitness_f(i) / fobj([zeros(1, one_six_dim), Xfood]);
                if P > 2
                    Xfood = exp(-1 / P) .* Xfood;
                    for j = 1:dim - one_six_dim
                        Xnew_coa(i, j) = X_coa(i, j) + cos(2 * pi * rand) * Xfood(j) * p_obj(temp) - sin(2 * pi * rand) * Xfood(j) * p_obj(temp);
                    end
                else
                    Xnew_coa(i, :) = (X_coa(i, :) - Xfood) * p_obj(temp) + p_obj(temp) .* rand(1, dim - one_six_dim) .* X_coa(i, :);
                end
            end
        end

        %% Boundary conditions for COA
        for i = 1:N
            for j = 1:dim - one_six_dim
                if length(ub_coa) == 1
                    Xnew_coa(i, j) = min(ub_coa, Xnew_coa(i, j));
                    Xnew_coa(i, j) = max(lb_coa, Xnew_coa(i, j));
                else
                    Xnew_coa(i, j) = min(ub_coa(j), Xnew_coa(i, j));
                    Xnew_coa(i, j) = max(lb_coa(j), Xnew_coa(i, j));
                end
            end
        end

        global_position_coa = Xnew_coa(1, :);
        global_fitness_coa = fobj([zeros(1, one_six_dim), global_position_coa]);

        for i = 1:N
            new_fitness = fobj([zeros(1, one_six_dim), Xnew_coa(i, :)]);
            if new_fitness < global_fitness_coa
                global_fitness_coa = new_fitness;
                global_position_coa = Xnew_coa(i, :);
            end
            if new_fitness < fitness_f(i)
                fitness_f(i) = new_fitness;
                X_coa(i, :) = Xnew_coa(i, :);
                if fitness_f(i) < BestFitness_coa
                    BestFitness_coa = fitness_f(i);
                    best_position_coa = X_coa(i, :);
                end
            end
        end
        global_Cov(Iter) = global_fitness_coa;

        % Update convergence curve
        if BestFitness_mgo < BestFitness_coa
            cnvg(Iter) = BestFitness_mgo;
        else
            cnvg(Iter) = BestFitness_coa;
        end

        %% Hybrid Solution Exchange
        
        % Select bottom 50% solutions from each algorithm to exchange
        num_exchange_1 = ceil(0.9 * N);
        num_exchange_2 = ceil(0.1 * N);
        % Sort MGO population
        [~, idx_mgo] = sort(Sol_Cost_mgo, 'descend');
        best_mgo_solutions = X_mgo(idx_mgo(1:num_exchange_1), :);
        % Sort COA population
        [~, idx_coa] = sort(fitness_f, 'descend');
        best_coa_solutions = X_coa(idx_coa(1:num_exchange_2), :);
        
        
        % Exchange solutions
        X_mgo(idx_mgo(1:num_exchange_2), :) = best_coa_solutions(:, 1:one_six_dim);
        X_coa(idx_coa(1:num_exchange_1), 1:one_six_dim) = best_mgo_solutions(:, 1:one_six_dim);
    end

    % Compare the best solutions from MGO and COA and display the overall best
    if BestFitness_mgo < BestFitness_coa
        BestF = BestFitness_mgo;
        BestX = [BestX_mgo, zeros(1, dim - one_six_dim)];
    else
        BestF = BestFitness_coa;
        BestX = [zeros(1, one_six_dim), best_position_coa];
    end
end

function y = p_obj(x)
    y = 0.2 * (1 / (sqrt(2 * pi) * 3)) * exp(-(x - 25).^2 / (2 * 3.^2));
end
