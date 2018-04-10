function []=InverseDynamics()
    clc;
    clear all;
    close all;
    % the following parameters for the arm
    % initial condition
    x0= [-0.5,0.2,0.1,0.1];
    x1 =[-0.8,0.5,0.1,0.1];
    x_i=[x0;x1];
    w = 0.2;
    tf = 10; 
    I1=10;  I2 = 10; m1=5; r1=.5; m2=5; r2=.5; l1=1; l2=1;
    % we compute the parameters in the dynamic model
    a = I1+I2+m1*r1^2+ m2*(l1^2+ r2^2);
    b = m2*l1*r2;
    d = I2+ m2*r2^2;
    
    global torque
    torque=[];
    for i=1:1:2
       options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4]);
       [T,X] = ode45(@(t,x)InverseDynamics(t,x),[0 tf],x_i(i,:), options);


       % Plot Data
       figure('Name','Theta1 under Computed Torque Control');
       plot(T, X(:,1),'r-');
       hold on
       plot(T, w*ones(size(T,1),1),'b-');
       figure('Name','Theta2 under Computed Torque Control');
       plot(T, X(:,2),'r--');
       hold on
       plot(T, sin(2*T),'b-');
    end





    function dx=InverseDynamics(t,x)
        theta_d = [0.2;sin(2*t)]; % Desired setpoint position
        dtheta_d = [0;2*cos(2*t)]; % Desired velocity 
        ddtheta_d= [0;-4*sin(2*t)];
        theta = x(1:2,1);
        dtheta= x(3:4,1);
       
        % create symbolic variable for x.
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
        tau = subs(M,symx,x) * (Kp*e + Kv*de)+ subs(C,symx,x)*dtheta + subs(M,symx,x) *ddtheta_d;
        torque =[torque, tau];
        dx=zeros(4,1);
        dx(1) = x(3); %dtheta1
        dx(2) = x(4); %dtheta2
        dx(3:4) = -subs(invMC,symx,x)* x(3:4) +subs(invM,symx,x)*tau; % because ddot theta = -M^{-1}(C \dot Theta) + M^{-1} tau
       
    end
end