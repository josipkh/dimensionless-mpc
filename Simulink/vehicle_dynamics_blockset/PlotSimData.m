clc;close all
data_file = 'data/sine_steer_pi_45';

if ~isempty(data_file)
    clearvars -except data_file VEHICLE SIM
    load(data_file)
else
    data = logsout;
    tsample = SIM.T;
end

if length(data) > 1
    data = data{1};
end

t =  data{1}.Values.Time;
vx = data{1}.Values.Data;
yaw_rate = data{2}.Values.Data;
wheel_speed = data{3}.Values.Data;
reference = data{4}.Values.Data;
inputs = squeeze(data{5}.Values.Data);
steer = data{6}.Values.Data;

if ~isempty(data_file) && length(t) > length(inputs)
    tsample = 0 : round(t(end)/length(inputs),2) : t(end);
    assert(length(tsample)==length(inputs))
else
    tsample = t;
end

figure
subplot(4,1,1); plot(tsample,steer); title('Steering angle [deg]'); ylim(get(gca,'YLim')+[-1,1])
%     hold on; yline(180/pi*VEHICLE.MAX_STEERING_ANGLE,'--r'); yline(-180/pi*VEHICLE.MAX_STEERING_ANGLE,'--r')
subplot(4,1,2); plot(t,vx); title('Longitudinal speed [km/h]'); ylim(get(gca,'YLim')+[-1,1])
    hold on; plot(tsample,reference(:,1),'r--')
subplot(4,1,3); plot(t,yaw_rate); title('Yaw rate [deg/s]'); ylim(get(gca,'YLim')+[-10,10])
    hold on; plot(tsample,reference(:,3),'r--')
subplot(4,1,4); stairs(tsample,inputs'); title('Motor torques [Nm]'); ylim(get(gca,'YLim')+[-1,1])
    legend('T_{fl}','T_{fr}','T_{rl}','T_{rr}')
xlabel('Time [s]')

