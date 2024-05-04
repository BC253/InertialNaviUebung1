function [y, w_ieP,RPY] = RK2(IMUdata, y, Cpe, Omega_iee, w_ieP, h)

   
    aP = IMUdata(:,1:3)';
   %  % 初始化输出数组
    w_ipP = IMUdata(:,4:6)';

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
    %DGL, bzw. f(yn,tn), bzw. k1
    k1x = y(4:6,i);
    k1v = Cpe*ap-2*Omega_iee*y(4:6,i)-Omega_iee*Omega_iee*y(1:3,i);
    k1q = 1/2*[0 wip1p-wie1p wip2p-wie2p wip3p-wie3p;
              -wip1p+wie1p 0 wip3p-wie3p -wip2p+wie2p;
              -wip2p+wie2p -wip3p+wie3p 0 wip1p-wie1p;
              -wip3p+wie3p wip2p-wie2p -wip1p+wie1p 0]*y(7:10,i);
    % y(7:10,i+1) = y(7:10,i)+h*(k1q);
    % %Normierung der Quaternionen
    % y(7:10,i+1) = y(7:10,i+1)/sqrt(y(7:10,i+1)'*y(7:10,i+1));
    %Cpe aus neuen Quaternionen berechnen
  
    k2x = y(4:6,i)+h*k1v;
    k2v = Cpe*ap-2*Omega_iee*(y(4:6,i)+h*k1v)-Omega_iee*Omega_iee*(y(1:3,i)+h*k1x);
    k2q = 1/2*[0 wip1p-wie1p wip2p-wie2p wip3p-wie3p;
              -wip1p+wie1p 0 wip3p-wie3p -wip2p+wie2p;
              -wip2p+wie2p -wip3p+wie3p 0 wip1p-wie1p;
              -wip3p+wie3p wip2p-wie2p -wip1p+wie1p 0]*(y(7:10,i)+k1q*h);
          
    %numerisch Integrierte Werte RK 2.Ordnung
    y(1:3,i+1) = y(1:3,i)+h/2*(k1x+k2x);
    y(4:6,i+1) = y(4:6,i)+h/2*(k1v+k2v);
    y(7:10,i+1) = y(7:10,i)+h/2*(k1q+k2q);
    
    %Normierung der Quaternionen
    y(7:10,i+1) = y(7:10,i+1)/sqrt(y(7:10,i+1)'*y(7:10,i+1));
    %Cpe aus neuen Quaternionen berechnen
    q0 = y(7,i+1);
    q1 = y(8,i+1);
    q2 = y(9,i+1);
    q3 = y(10,i+1);
    Cpe = [q0^2+q1^2-q2^2-q3^2 2*(q1*q2+q3*q0) 2*(q1*q3-q2*q0);
           2*(q1*q2-q3*q0) q0^2-q1^2+q2^2-q3^2 2*(q2*q3+q1*q0);
           2*(q1*q3+q2*q0) 2*(q2*q3-q1*q0) q0^2-q1^2-q2^2+q3^2];  %InS 3.3
    
    
    % 假设C_pe是按照Z-Y-X顺序构造的旋转矩阵
    RPY(1,i) = atan2(Cpe(3,2), Cpe(3,3));
    RPY(2,i) = asin(-Cpe(3,1));
    RPY(3,i) = atan2(Cpe(2,1), Cpe(1,1));


  end

end

% function [qD] = q_dot(q, w_ipP, w_ieP)
%     % Berechnung der Ableitung der Quaternionen
%     w_ipP1 = w_ipP(1);
%     w_ipP2 = w_ipP(2);
%     w_ipP3 = w_ipP(3);
%     w_ieP1 = w_ieP(1);
%     w_ieP2 = w_ieP(2);
%     w_ieP3 = w_ieP(3);
%     q0 = q(1);
%     q1 = q(2);
%     q2 = q(3);
%     q3 = q(4);
% 
%     qD = 1/2 * [0, w_ipP1-w_ieP1, w_ipP2-w_ieP2, w_ipP3-w_ieP3;
%                -w_ipP1+w_ieP1, 0, w_ipP3-w_ieP3, -w_ipP2+w_ieP2;
%                -w_ipP2+w_ieP2, -w_ipP3+w_ieP3, 0, w_ipP1-w_ieP1;
%                -w_ipP3+w_ieP3, w_ipP2-w_ieP2, -w_ipP1+w_ieP1, 0] * [q0; q1; q2; q3];
% end
