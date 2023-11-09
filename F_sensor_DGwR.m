function [time_DGwR, H, sensors] = F_sensor_DGwR(U, p,sigma_o2)

    sensors=[]; %ok?
    [n,r]=size(U);
    if p <= r
        tic;
        [sensors,Rinv,det_test] = F_sensor_DGwR_r(U,[],sigma_o2,p,sensors,[]);
        % time_DGwR=toc;
        [H] = F_calc_sensormatrix(p, n, sensors);
        time_DGwR=toc;
    else
        tic;
        [isensors,Rinv,det_test] = F_sensor_DGwR_r(U,[],sigma_o2,r,sensors,[]);
        % [H] = F_calc_sensormatrix(r, n, isensors);
        [sensors] = F_sensor_DGwR_p(U,[],sigma_o2,p,isensors,Rinv); %H?
        % time_DGwR=toc;
        [H] = F_calc_sensormatrix(p, n, sensors);
        time_DGwR=toc;
    end
    text = ['S2SwR-based : sensor=',mat2str(sensors)];
    disp(text)

end