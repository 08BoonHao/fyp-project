function [NewX , Sol_CostNew] = Boundary_Check1(NewX,fobj,LB,UB,best_position_coa)

    for j=1:4
        NewX(j,:) = boundaryCheck(NewX(j,:), LB, UB);
        Sol_CostNew(j,:)=fobj([NewX(j,:),best_position_coa]);%#ok
    end
            
end