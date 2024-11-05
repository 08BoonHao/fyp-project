function [fitness]=ObjFunction(UAV_num,UE_num,Pos_U,PopPos)
    %% Coefficient Settings
    beta_zero = 0.001;
    B = 1;
    P = 10;
    gamma_square = 1e-12;
    max_users_per_uav = 100;
    
    %% PROBLEM start
    UAV_UE_distance = zeros(UAV_num, UE_num);
    
    QOS = P * beta_zero / gamma_square;
    SNR = QOS / (135.^2 + 140);
    
    for d = 1:UAV_num
        n = 3 * d;
        % --- distance of UAV with UE
        for u = 1:UE_num
            manhatten_distance = abs(PopPos(n-2) - Pos_U(u,1)) + abs(PopPos(n-1) - Pos_U(u,2));
            UAV_UE_distance(d,u) = manhatten_distance.^2 + PopPos(n).^2;
        end
    end
    
    channel_gain = beta_zero ./ UAV_UE_distance;
    % disp(['Channel gain: ' num2str(channel_gain)]);
    capacity = B * log2(1 + P * channel_gain / gamma_square);
    maximum_coverage_radius = sqrt(QOS / SNR - PopPos(3:3:end)');
    
    %% Greedy algorithm
    A = zeros(UAV_num,UE_num);
    users_associated_with_uav = zeros(UAV_num,1);
    
    for u = 1:UE_num
        % Sort UAVs by capacity in descending order for user 'u'
        [sorted_capacities, sorted_idx] = sort(capacity(:,u), 'descend');
    
        % Try to assign user to the best UAVs
        assigned = false;
        for i = 1:length(sorted_idx)
            idx = sorted_idx(i);
            if users_associated_with_uav(idx) < max_users_per_uav
                A(idx, u) = 1;
                users_associated_with_uav(idx) = users_associated_with_uav(idx) + 1;
                assigned = true;
                break;
            end
        end
    end
        
    %% Calculation for fitness
    fitness = sum(sum((A .* capacity))) / UAV_num;
    
end