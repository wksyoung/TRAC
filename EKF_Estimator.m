function obj = EKF_Estimator(InitState)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
obj = extendedKalmanFilter(@state_trans,@measureFunc,InitState,'HasAdditiveProcessNoise',false);%unscentedKalmanFilter
obj.MeasurementNoise = 0.0001;
obj.ProcessNoise = 0.8*eye(3);
end

function x = state_trans(x,w,u)
    global StM;
    S = [eye(3) zeros(3)];
    D = [zeros(3) eye(3)];
    R = diag([-0.4,-0.2,-0.2])*S*x;
    M = StM;
    new_s = S*x + (R + 250*M*w)*0.02;%new_s = S*x + (R + 25*M*w)*0.02;
%     new_f = WA*F*x + WB*w;
%     wy = WC*new_f;
    new_d = D*x + (-0.5*M'*S*x + u + 150*w)*0.02;%new_d = D*x + (-0.5*M'*S*x + u + 15*w)*0.02;
    x = [new_s;new_d];%new_f
end

function y = measureFunc(x,u)
%    C = [eye(4),zeros(4,20)];
    %C = [eye(4),zeros(4)];
    y = x;
end

% function output=satu(S)
%     if S>1
%         output=1;
%     elseif S<-1
%         output=-1;
%     else
%         output=S;
%     end
% end
