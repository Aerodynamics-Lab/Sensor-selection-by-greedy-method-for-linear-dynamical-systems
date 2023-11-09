function [A1, A2, A3, A4, A5, A6] = F_data_ave2(CNT_num, num_ave, A1, A2, A3, A4, A5, A6)

    for CNT=1:1:CNT_num
        A1(CNT,1) = mean(A1(CNT,2:num_ave+1));
        A2(CNT,1) = mean(A2(CNT,2:num_ave+1));
        A3(CNT,1) = mean(A3(CNT,2:num_ave+1));
        A4(CNT,1) = mean(A4(CNT,2:num_ave+1));
        A5(CNT,1) = mean(A5(CNT,2:num_ave+1));
        A6(CNT,1) = mean(A6(CNT,2:num_ave+1));
    end
    % A1(:,1) = mean(A1(:,2:num_ave+1));
    % A2(:,1) = mean(A2(:,2:num_ave+1));
    % A3(:,1) = mean(A3(:,2:num_ave+1));

end
