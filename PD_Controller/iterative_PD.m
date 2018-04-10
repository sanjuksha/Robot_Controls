function []=iterative()
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
    I1=10;  I2 = 10; m1=5; r1=.5; m2=5; r2=.5; l1=1; l2=1;g=9.8;
    % we compute the parameters in the dynamic model
    a = I1+I2+m1*r1^2+ m2*(l1^2+ r2^2);
    b = m2*l1*r2;
    d = I2+ m2*r2^2;
    
    global torque
    torque=[];
       options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4]);
       [T,X] = ode45(@(t,x)iterative(t,x),[0 tf],x0, options);


       % Plot Data
       figure('Name','Theta1 under Iterative Control Scheme');
       plot(T, X(:,1),'r-');
       hold on
       plot(T, w*ones(size(T,1),1),'b-');
       figure('Name','Theta2 under Iterative Control Scheme');
       plot(T, X(:,2),'r--');
       hold on
       plot(T,w*ones(size(T,1),1) ,'b-');
    





    function dx=iterative(t,x)
        theta_d = [0.2;0.2]; % Desired setpoint position
        dtheta_d = [0;0]; % Desired velocity 
        ddtheta_d= [0;0];
        theta = x(1:2,1);
        dtheta= x(3:4,1);
       
        % create symbolic variable for x.
        global M C u Gq d_Gq
        symx= sym('symx',[4,1]);
        M = [a+2*b*cos(symx(2)), d+b*cos(symx(2));d+b*cos(symx(2)), d];
        C = [-b*sin(symx(2))*symx(4), -b*sin(symx(2))*(symx(3)+symx(4)); b*sin(symx(2))*symx(3),0];
        invM = inv(M);
        invMC= invM*C;
        % Gravity Matrix 
        g1=-(m1+m2)*g*l1*sin(symx(2))-m2*g*l2*sin(symx(1)+symx(2)); 
        g2=-m2*g*l2*sin(symx(1)+symx(2)); 
        Gq=[g1;g2];
        % To Calculate alpha we compute derivative of Gq
%         d_Gq1 = diff(Gq,symx(1));
%         d_Gq2 = diff(Gq,symx(2));
%         d_Gq  =[d_Gq1,d_Gq2];
%        
%        alpha = double(subs(d_Gq,[symx(1);symx(2)],[x(1);x(2)]));
        Kp = 200*eye(2);
        Kv = 30*eye(2);
        e = theta_d - theta; % position error
        de = dtheta_d - dtheta; % velocity error
        beta = 0.1;
        
        if (t == 0)
            u =  zeros(2,1) ;
        end
        tau = (1/beta)*(Kp*e + Kv*de)+ u;  
     if (norm(dtheta)< 0.0001)
        u = tau
        display(t)
      end
        
        torque =[torque, tau];
        dx=zeros(4,1);
        dx(1) = x(3); %dtheta1
        dx(2) = x(4); %dtheta2
        dx(3:4) = -subs(invMC,symx,x)* x(3:4) +subs(invM,symx,x)*tau - subs(invM,symx,x)*subs(Gq,symx,x); % because ddot theta = -M^{-1}(C \dot Theta) + M^{-1} tau
       
    end
end