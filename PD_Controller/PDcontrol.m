
function []=PDcontrol()
    % Notations: For a given variable, x, dx is its time derivative, ddx is
    % 2nd-order derivative.
    clc;
    clear all;
    close all;
    x0= [-0.5,0.2,0.1,0.1];
    x1 =[-0.8,0.6,0.1,0.1];
    x_i=[x0;x1];
    tf = 10;
    % Implement the PD control for set point tracking.
    I1=10;  I2 = 10; m1=5; r1=.5; m2=5; r2=.5; l1=1; l2=1;
    % we compute the parameters in the dynamic model
    a = I1+I2+m1*r1^2+ m2*(l1^2+ r2^2);
    b = m2*l1*r2;
    d = I2+ m2*r2^2;
    
    global torque
    torque=[];
    for i=1:1:2
        options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4]);
        [T,X] = ode45(@(t,x)PDControl(t,x),[0 tf],x_i(i,:), options);

        % Plot Data
        figure('Name','Theta_1 under PD SetPoint Control');
        plot(T, X(:,1),'r-');
        hold on

        figure('Name','Theta_2 under PD SetPoint Control');
        plot(T, X(:,2),'r--');
        hold on
    end

    function dx=PDControl(t,x)
        theta_d = [0;0]; % Desired setpoint position
        dtheta_d = [0;0]; % Desired velocity 
        ddtheta_d= [0;0];
        theta = x(1:2,1);
        dtheta= x(3:4,1);
       
        global M C
        symx= sym('symx',[4,1]);
        M = [a+2*b*cos(symx(2)), d+b*cos(symx(2));d+b*cos(symx(2)), d];
        C = [-b*sin(symx(2))*symx(4), -b*sin(symx(2))*(symx(3)+symx(4)); b*sin(symx(2))*symx(3),0];
        invM = inv(M);
        invMC= invM*C;
        
        Kp = 80*eye(2);
        Kv = 50*eye(2);
        e = theta_d - theta; % position error
        de = dtheta_d - dtheta; % velocity error
        tau = Kp*e + Kv*de;
        torque =[torque, tau];
        dx=zeros(4,1);
        dx(1) = x(3); %dtheta1
        dx(2) = x(4); %dtheta2
        dx(3:4) = -subs(invMC,symx,x)* x(3:4) +subs(invM,symx,x)*tau; % because ddot theta = -M^{-1}(C \dot Theta) + M^{-1} tau
       
    end
end