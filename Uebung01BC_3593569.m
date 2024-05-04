clear all
close all
clc
% Aufgabe1 BC

% Erdrotation
w_E = 2*pi/24/60/60;
% w_E = 7.292115*1e-5;
Omega_iee = [0 -w_E 0
             w_E 0 0
             0 0 0];


% 读取IMU数据 ax ay az(m/s^2) gx gy gz(rad/s)
IMU.raw = csvread('uebung01-IMU.csv', 1, 0);
IMUdata(:,1:3) = IMU.raw(:,2:4);
IMUdata(:,4:6) = IMU.raw(:,5:7);

REF.raw = csvread('uebung01-REF.csv', 1, 0);
refp = REF.raw(:,2:4);
refv = REF.raw(:,5:7)';

%角度转弧度 deg2rad
refRPY = deg2rad(REF.raw(:,8:10));
refp(:,1:2) = deg2rad(refp(:,1:2)); 

%Ref Position in e-system
wgs84 = wgs84Ellipsoid('meter');
[xRef(:),yRef(:),zRef(:)] = geodetic2ecef(refp(:,1),refp(:,2),refp(:,3),wgs84);



%DCM berechnen
Cne = C(3,-refp(1,2))*C(2,refp(1,1)+pi/2);
Cbn = C(3,-refRPY(1,3))*C(2,-refRPY(1,2))*C(1,-refRPY(1,1));
Cpe = Cne*Cbn;

Omega_iep = inv(Cpe)*Omega_iee*Cpe;

w_ieP(1,1) = Omega_iep(3,2);
w_ieP(2,1) = Omega_iep(1,3);
w_ieP(3,1) = Omega_iep(2,1);




refv = Cne*refv;                                                  %%%%%

% start pos
x0(1,1) = xRef(1);
x0(1,2) = yRef(1);
x0(1,3) = zRef(1);
x0(1,4:6) = refv(1:3,1);                                        %%%%%%%                       


%Quaternion im e-system berechnen
q0 = 0.5*sqrt(Cpe(1,1)+Cpe(2,2)+Cpe(3,3)+1);
q1 = (Cpe(2,3)-Cpe(3,2))/(4*q0);
q2 = (Cpe(3,1)-Cpe(1,3))/(4*q0);
q3 = (Cpe(1,2)-Cpe(2,1))/(4*q0);    %这里之前+号应该改成*号

y0 = [x0(1,1:3)';x0(1,4:6)';q0;q1;q2;q3];
h = 0.1;
% % RK2 Heun [y, w_ieP] = RK2(imudata, y, Cpe, Omega_ieE, w_ieP, h)
[y_2, w_ieP_2,RPY_2] = RK2(IMUdata,y0,Cpe,Omega_iee,w_ieP,0.1);



%% Aufgabe b


figure()
plot3(xRef, yRef,zRef)
hold on
plot3(y_2(1,:),y_2(2,:),y_2(3,:))
legend('Referenz','RK2')
title('Aufgabe b')



figure()
plot(yRef,zRef)
hold on
plot(y_2(2,:),y_2(3,:))
legend('Referenz','RK2')
title('Aufgabe b')

%% Aufgabe c
t = IMU.raw(:,1);
figure
plot(t,y_2(3,1:end-1))
hold on
plot(t,zRef(:))
xlabel('Zeit [s]')
ylabel('Höhe [m]')
legend('RK2','Referenz')
title('Höhe')

figure
plot(t,y_2(4,1:end-1))
hold on
plot(t,refv(1,:))
xlabel('Zeit [s]')
ylabel('v_N [m/s]')
legend('RK2','Referenz')
title('Geschwindigkeit in Nord')

figure
plot(t,y_2(5,1:end-1))
hold on
plot(t,refv(2,:))
xlabel('Zeit [s]')
ylabel('v_E [m/s]')
legend('RK2','Referenz')
title('Geschwindigkeit in Ost')

%% Aufgabe d

figure()
subplot(3,1,1)
plot(t,RPY_2(1,:))
hold on 
plot(t,refRPY(:,1))
xlabel('Zeit [s]')
ylabel('Roll [deg]')
title('Roll')
legend('RK2','Referenz')

subplot(3,1,2)
plot(t,RPY_2(2,:))
hold on 
plot(t,refRPY(:,2))
xlabel('Zeit [s]')
ylabel('Pitch [deg]')
title('Pitch')
legend('RK2','Referenz')

subplot(3,1,3)
plot(t,RPY_2(3,:))
hold on 
plot(t,refRPY(:,3))
xlabel('Zeit [s]')
ylabel('Yaw [deg]')
title('Yaw')
legend('RK2','Referenz')

%% Aufgabe 2 RK1
[y_1, w_ieP_1,RPY_1] = RK1(IMUdata,y0,Cpe,Omega_iee,w_ieP,0.1);
[y_3, w_ieP_3,RPY_3] = RK3(IMUdata,y0,Cpe,Omega_iee,w_ieP,0.1);


figure
plot3(xRef, yRef, zRef)
hold on
plot3(y_1(1,:), y_1(2,:), y_1(3,:),'b-','LineWidth', 2)
plot3(y_2(1,:), y_2(2,:), y_2(3,:),'r--')
% plot3(y_3(1,:), y_3(2,:), y_3(3,:))
legend('Referenz','RK1','RK2','RK3')


figure
plot( yRef, zRef)
hold on
plot( y_1(2,:), y_1(3,:),'b-','LineWidth', 2)
plot( y_2(2,:), y_2(3,:),'r--')
% plot( y_3(2,:), y_3(3,:))
legend('Referenz','RK1','RK2','RK3')



figure
plot(t, y_1(3,1:end-1),'b-','LineWidth', 2 ,'DisplayName', 'RK1')
hold on
plot(t, y_2(3,1:end-1),'r--', 'DisplayName', 'RK2')
plot(t, y_3(3,1:end-1), 'DisplayName', 'RK3')
plot(t, zRef(:), 'DisplayName', 'Referenz')
xlabel('Zeit [s]')
ylabel('Höhe [m]')
legend show
title('Höhe Comparison')
hold off


figure
plot(t, y_1(4,1:end-1), 'DisplayName', 'RK1')
hold on
plot(t, y_2(4,1:end-1), 'DisplayName', 'RK2')
plot(t, y_3(4,1:end-1), 'DisplayName', 'RK3')
plot(t, refv(1,:), 'DisplayName', 'Referenz')
xlabel('Zeit [s]')
ylabel('v_N [m/s]')
legend show
title('Geschwindigkeit in Nord Comparison')
hold off


figure
plot(t, y_1(5,1:end-1), 'DisplayName', 'RK1')
hold on
plot(t, y_2(5,1:end-1), 'DisplayName', 'RK2')
plot(t, y_3(5,1:end-1), 'DisplayName', 'RK3')
plot(t, refv(2,:), 'DisplayName', 'Referenz')
xlabel('Zeit [s]')
ylabel('v_E [m/s]')
legend show
title('Geschwindigkeit in Ost Comparison')
hold off


figure
subplot(3,1,1)
plot(t, RPY_1(1,:), 'DisplayName', 'RK1')
hold on 
plot(t, RPY_2(1,:), 'DisplayName', 'RK2')
plot(t, RPY_3(1,:), 'DisplayName', 'RK3')
plot(t, refRPY(:,1), 'DisplayName', 'Referenz')
xlabel('Zeit [s]')
ylabel('Roll [deg]')
title('Roll Comparison')
legend show
hold off

subplot(3,1,2)
plot(t, RPY_1(2,:), 'DisplayName', 'RK1')
hold on 
plot(t, RPY_2(2,:), 'DisplayName', 'RK2')
plot(t, RPY_3(2,:), 'DisplayName', 'RK3')
plot(t, refRPY(:,2), 'DisplayName', 'Referenz')
xlabel('Zeit [s]')
ylabel('Pitch [deg]')
title('Pitch Comparison')
legend show
hold off

subplot(3,1,3)
plot(t, RPY_1(3,:), 'DisplayName', 'RK1')
hold on 
plot(t, RPY_2(3,:), 'DisplayName', 'RK2')
plot(t, RPY_3(3,:), 'DisplayName', 'RK3')
plot(t, refRPY(:,3), 'DisplayName', 'Referenz')
xlabel('Zeit [s]')
ylabel('Yaw [deg]')
title('Yaw Comparison')
legend show
hold off

% RK2 VS RK1 
y_2_1 = y_2-y_1;
RPY_2_1 = RPY_2-RPY_1;

figure
subplot(3,1,1)
plot(t, y_2_1(3,1:end-1))
xlabel('Zeit [s]')
ylabel('Höhe [m]')
title('RK2 VS RK1 in Höhe')

subplot(3,1,2)
plot(t, y_2_1(3,1:end-1))
xlabel('Zeit [s]')
ylabel('v_N [m/s]')
title('RK2 VS RK1 Geschwindigkeit in Nord')

subplot(3,1,3)
plot(t, y_2_1(3,1:end-1))
xlabel('Zeit [s]')
ylabel('v_E [m/s]')
title('RK2 VS RK1 Geschwindigkeit in Ost')


figure
subplot(3,1,1)
plot(t, RPY_2_1(1,:))
xlabel('Zeit [s]')
ylabel('Roll [deg]')
title('RK2 VS RK1 in Roll')

subplot(3,1,2)
plot(t, RPY_2_1(2,:))
xlabel('Zeit [s]')
ylabel('Pitch [deg]')
title('RK2 VS RK1 in Pitch')

subplot(3,1,3)
plot(t, RPY_2_1(3,:))
xlabel('Zeit [s]')
ylabel('Yaw [deg]')
title('RK2 VS RK1 in Yaw')


