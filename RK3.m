function [y, w_ieP, RPY] = RK3(IMUdata, y, Cpe, Omega_iee, w_ieP, h)

    % Messungen der Beschleunigungen
    aP = IMUdata(:,1:3)';
    w_ipP = IMUdata(:,4:6)';

    % 对IMU数据进行积分
    for i = 1:length(aP)
        Omegaiep = inv(Cpe)*Omega_iee*Cpe;
        wie1p = Omegaiep(3,2);
        wie2p = Omegaiep(1,3);
        wie3p = Omegaiep(2,1);
        %Messungen:
        wip1p = IMUdata(i,4);   %Rotationen
        wip2p = IMUdata(i,5);
        wip3p = IMUdata(i,6);
        ap = [IMUdata(i,1);
          IMUdata(i,2);
          IMUdata(i,3)];  %Beschleunigungen
        k1x = y(4:6,i);
        k1v = Cpe * aP(:,i) - 2*Omega_iee*y(4:6,i)-Omega_iee*Omega_iee*y(1:3,i);
        k1q = q_dot(y(7:10,i),w_ipP(:,i),w_ieP);


        k2x = y(4:6,i) + h/2*k1v;
        k2v = Cpe * aP(:,i) - 2*Omega_iee*(y(4:6,i)+h/2*k1v)-Omega_iee*Omega_iee*(y(1:3,i)+h/2*k1x);
        k2q = q_dot(y(7:10,i) + h/2*k1q, w_ipP(:,i),w_ieP);

        k3x = y(4:6,i) - h*k1v + 2*h*k2v;
        k3v = Cpe * aP(:,i) - 2*Omega_iee*(y(4:6,i)-h*k1v+2*h*k2v)-Omega_iee*Omega_iee*(y(1:3,i)-h*k1x+2*h*k2x);
        k3q = q_dot(y(7:10,i) - h*k1q + 2*h*k2q, w_ipP(:,i),w_ieP);

        y(:,i+1) = y(:,i) + h*(1/6*[k1x; k1v; k1q] + 2/3*[k2x; k2v; k2q] + 1/6*[k3x; k3v; k3q]); % RK3 update
      %Cpe aus neuen Quaternionen berechnen
        q0 = y(7,i+1);
        q1 = y(8,i+1);
        q2 = y(9,i+1);
        q3 = y(10,i+1);
        Cpe = [q0^2+q1^2-q2^2-q3^2 2*(q1*q2+q3*q0) 2*(q1*q3-q2*q0);
           2*(q1*q2-q3*q0) q0^2-q1^2+q2^2-q3^2 2*(q2*q3+q1*q0);
           2*(q1*q3+q2*q0) 2*(q2*q3-q1*q0) q0^2-q1^2-q2^2+q3^2];  %InS 3.3
        RPY(1,i) = atan2(Cpe(3,2), Cpe(3,3));
        RPY(2,i) = asin(-Cpe(3,1));
        RPY(3,i) = atan2(Cpe(2,1), Cpe(1,1));

        % Update rotation matrix and Euler angles
        % (Code for updating Cpe and calculating RPY remains the same)
    end
end
function [qD] = q_dot(q, w_ipP, w_ieP)
    % Berechnung der Ableitung der Quaternionen
    w_ipP1 = w_ipP(1);
    w_ipP2 = w_ipP(2);
    w_ipP3 = w_ipP(3);
    w_ieP1 = w_ieP(1);
    w_ieP2 = w_ieP(2);
    w_ieP3 = w_ieP(3);
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    qD = 1/2 * [0, w_ipP1-w_ieP1, w_ipP2-w_ieP2, w_ipP3-w_ieP3;
               -w_ipP1+w_ieP1, 0, w_ipP3-w_ieP3, -w_ipP2+w_ieP2;
               -w_ipP2+w_ieP2, -w_ipP3+w_ieP3, 0, w_ipP1-w_ieP1;
               -w_ipP3+w_ieP3, w_ipP2-w_ieP2, -w_ipP1+w_ieP1, 0] * [q0; q1; q2; q3];
end