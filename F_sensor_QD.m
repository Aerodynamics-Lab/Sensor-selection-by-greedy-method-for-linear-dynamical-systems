function [time_QD, H, sensors]=F_sensor_QD(U, p)

    [n,r]=size(U);
    if p <= r
        tic;
        [sensors] = F_sensor_QR_pivot(p, U);
        [H] = F_calc_sensormatrix(p, n, sensors);
        time_QD = toc;
    else
        tic;
        [isensors] = F_sensor_QR_pivot(r, U);
        [H] = F_calc_sensormatrix(r, n, isensors);
        [sensors] = F_sensor_DG_p(U, p, H, isensors');
        [H] = F_calc_sensormatrix(p, n, sensors);
        time_QD = toc;
    end
    text = ['S2S-based : sensor=',mat2str(sensors)];
    disp(text)
    
end