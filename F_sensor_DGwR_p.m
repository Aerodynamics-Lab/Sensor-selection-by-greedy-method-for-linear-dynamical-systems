 function [sensors,Rinv,det_test]= F_sensor_DGwR_p(U,Unoi,Snoi,p,sensors,Rinv)
    [n,r1]=size(U);
    Unoi=zeros(n,r1); %!!!!!
%     [~,rn]=size(Un);
    [ps,~]=size(sensors);
    Snoi_sq=Snoi; %Snoi_sq=Snoi*Snoi;!!!!!
    Cpp=zeros(ps,r1);
    initial=true;
    det_test=zeros(p,1);
    if p>r1
    for pp=(ps+1):p
        if initial==true   % initialize W&Cpp
%             if pp>1
            for l=1:ps
                Cpp(l,:)=U(sensors(l,1),:); 
            end
            RC=Rinv*Cpp;
            W=Cpp'*RC;
            Winv=inv(W);
            initial=false;
        end
        det_vec=zeros(1,n);
        %% searching
        for nn=1:n
            u_i=U(nn,:);
            s=zeros(1,pp-1);
            t=0;
            for l=1:(pp-1)
%                         for i=1:rn   % if truncate higher modes from noise
%                             s(1,l)=s(1,l)+Unoi(nn,i)*Unoi(sensors(l,1),i)*S(i,i)*S(i,i);
%                         end
                % s(1,l)=Unoi(nn,:)*Snoi_sq*Unoi(sensors(l,1),:)';
                s(1,l)=0; %!!!!!
            end
                        
%                     for i=1:rn   % if truncate higher modes from noise
%                         t=t+Unoi(sensors(pp,1),i)*Unoi(sensors(pp,1),i)*S(i,i)*S(i,i);
%                     end
            % t=Unoi(nn,:)*Snoi_sq*Unoi(nn,:)';   % same as rank of Unoi
            t=Snoi; %!!!!!
                %% objective value
            diff=s*RC-u_i;
            dnm=t-s*Rinv*s';
            NUME=diff*Winv*diff';
%                 det_vec(1,nn)=det(eye(s)+inv(dnm)*NUME); 
%                 det_v=det(eye(s)+dnm\NUME);
            det_vec(1,nn)=1+dnm\NUME;   
        end
        
        for l=1:(pp-1)
            det_vec(1,sensors(l,1))=0;
        end
        [det_test(pp,1),sensors(pp,1)]=max(det_vec);   % argmaxdet
        
%%   Update Rinv&C after we get pp-th sensor  
        s=zeros(1,pp-1);
        t=0;
%         for dd=1:d
        u_i=U(sensors(pp,1),:);
        
        for l=1:(pp-1)
%             for i=1:rn   % if truncate higher modes from noise
%                 s(1,l)=s(1,l)+Unoi(nn,i)*Unoi(sensors(l,1),i)*S(i,i)*S(i,i);
%             end
            % s(1,l)=Unoi(sensors(pp,1),:)*Snoi_sq*Unoi(sensors(l,1),:)';
            s(1,l)=0; %!!!!!
        end
        % t=Unoi(sensors(pp,1),:)*Snoi_sq*Unoi(sensors(pp,1),:)';
        t=Snoi; %!!!!!
        
        diff=s*RC-u_i;
        dnm=t-s*Rinv*s';
        W=W+diff'*(dnm\diff);
        Winv=inv(W);

        Cpp=[Cpp;u_i];
        
        % sR=s*Rinv; %!!out
        % Rinv_new=zeros(pp,pp); %!!out
        % Rinv_new(1:pp-1,1:pp-1)=Rinv; %!!out
        % Rinv=Rinv_new+[sR';-1]*(dnm\[sR -1]); %!!out
        Rtemp = Snoi * eye( pp ); %!!!!!white noise
        Rinv = inv(Rtemp); %!!!!!
        RC=Rinv*Cpp;
    end
    else
        disp('something wrong with p')
    end
end