function [detWo]=F_calc_detWo3(p, H,Aorg,Borg,Corg,Dorg)

    [~,r]=size(Corg);
    C = H*Corg;
    D = H*Dorg;
    CTC = C'*C;
    vecCTC = reshape(CTC,[],1);
    % vecWo = - ( kron(Aorg',Aorg')-eye(r^2) ) \ vecCTC;
    vecWo = - pinv( kron(Aorg',Aorg')-eye(r^2) ) * vecCTC;%(2021/4/22)
    Gram = reshape(vecWo,[r,r]);
    maxrank = rank(Gram);
    
    % sys = ss(Aorg,Borg,C,D,[]); % sys=ss(Aorg,Borg,C,D);
    % G(:,:,i) = gram(sys,'o'); % 可観測グラム行列が求まるコマンド
    % maxrank = G(:,:,i);

    text = [ '  rank(Wo) = ', num2str(maxrank) ]; 
    disp(text);
    eps=1e-8;
    detWo = det(Gram+eps*eye(r,r));

end