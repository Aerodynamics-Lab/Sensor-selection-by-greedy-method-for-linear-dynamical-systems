% function [logdet]=F_calc_det(p, H, U)

%     [~,r]=size(U);
%     C = H*U;
%     if p <= r
%         logdet = log(det(C*C'));
%     else
%         logdet = log(det(C'*C));
%     end
% end
function [detCC]=F_calc_det(p, H, U)

    [~,r]=size(U);
    C = H*U;
    if p <= r
        detCC = det(C*C');
    else
        detCC = det(C'*C);
    end
end
