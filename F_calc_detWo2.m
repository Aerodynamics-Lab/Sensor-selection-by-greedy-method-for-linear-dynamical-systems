function [detWo,S]=F_calc_detWo2(p, H,Aorg,Borg,Corg,Dorg)

    [~,r]=size(Corg);
    C = H*Corg;
    D = H*Dorg;
    CTC = C'*C;
    vecCTC = reshape(CTC,[],1);
    % vecWo = - ( kron(Aorg',Aorg')-eye(r^2) ) \ vecCTC;
    vecWo = - pinv( kron(Aorg',Aorg')-eye(r^2) ) * vecCTC;%(2021/4/22)
    Gram = reshape(vecWo,[r,r]);
    [U,S,V] = svd(Gram);
    if p <= r
        maxrank = rank(Gram);
        % Ur=U(:        ,1:maxrank);
        Sr = S(1:maxrank,1:maxrank);
        % Vr=V(:        ,1:maxrank);
        detWo = det(Sr);
    else
        detWo = det(S);
    end
    
end