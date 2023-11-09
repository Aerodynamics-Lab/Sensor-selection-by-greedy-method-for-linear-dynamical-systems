function [time, H, isensor]=F_sensor_KF(Aorg,Borg,Corg,Dorg,p,sigma_s2,sigma_o2)

    tic
    [n,r]= size(Corg);
    Ctmp=Corg;
    obj = zeros(n, 1);
    isensor = zeros(p, 1);
    H = zeros(p, n);
    for pp=1:1:p
        if pp ~= 1 %not equal
            Cprev=H*Corg;
        else %k=1
            Cprev=[];
            H=zeros(1,n);
        end

        obj(:) = inf; % objmin=10^100;
        %% Noise setting
        Q = sigma_s2 * eye( r ); 
        R = sigma_o2 * eye( pp ); 
        for i=1:1:n
            %% Update C and system for next sensor selection
            C = [ Cprev ; Ctmp(i,:)];
            
            %% Steady covariance matrix (Steady Kalman filter)
            [Pst,~,~] = idare(Aorg', C', Q, R, [], []); %B*Q*B' (R2019a-)
            obj(i) = det(Pst); % obj2(i) = trace(Pst) ;
        end
        [~,isensor(pp)] = min(obj); % [~,Im] = min(obj3); isensor(pp) = Im;
        H(pp,isensor(pp)) = 1;
        text = [ 'Result(KF-based): isensor=' num2str(isensor(pp))];
        disp ( text ) ;
        Ctmp(isensor(pp),:) = zeros(1,r);

    end
    time=toc;
    isensor=isensor';

end
