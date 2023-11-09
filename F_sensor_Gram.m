function [time, H, isensor]=F_sensor_Gram(Aorg,Borg,Corg,Dorg,p,sigma_s2,sigma_o2)

    tic
    [n,r]= size(Corg);
    H = zeros(p, n);
    isensor = zeros(p, 1);
    ranks = zeros(n, 1);
    obj = zeros(n, 1);

    Ctmp=Corg;
    
    flag_full_rank = 0;
    for pp=1:1:p
        if pp ~= 1
            Cprev=H*Corg;
        else 
            Cprev=[];
            H=zeros(1,n);
        end
        
        for i=1:1:n
            C=[Cprev;Ctmp(i,:)];
            CTC = C'*C;
            G = dlyap(Aorg', CTC);
            if ~flag_full_rank
                s = svd(G);
                rankG = sum(s > max(size(Aorg)) * eps(norm(s,inf))); % same as rank in MATLAB
                ranks(i) = rankG;
                obj(i) = prod(s(1:rankG)); % determinant of G in observable subspace
            else
                obj(i) = det(G);
            end
        end

        maxrank = max(ranks);
        if ~flag_full_rank
            obj(ranks < maxrank) = 0; % exclude candidates without max rank
        end
        [~,isensor(pp)] = max(obj);

        H(pp,isensor(pp)) = 1;
        text = [ 'Result(Gramian-based): isensor=' num2str(isensor(pp))];
        disp ( text ) ;
        Ctmp(isensor(pp),:) = zeros(1,r); %zero reset
        if maxrank == r
            flag_full_rank = 1;
        end
    end
    time=toc;
    isensor=isensor';
end
