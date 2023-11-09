function [STD1, STD2, STD3] = F_data_std1(CNT, num_ave, A1, A2, A3)
    
    % STD1(CNT,1) = std(A1(CNT,2:num_ave+1));
    % STD2(CNT,1) = std(A2(CNT,2:num_ave+1));
    % STD3(CNT,1) = std(A3(CNT,2:num_ave+1));
    STD1 = std(A1(:,2:num_ave+1),0,2);
    STD2 = std(A2(:,2:num_ave+1),0,2);
    STD3 = std(A3(:,2:num_ave+1),0,2);
    
end
