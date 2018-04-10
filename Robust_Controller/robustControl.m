function []= robustControl(theta10,theta20,dtheta10, dtheta20,theta1f, theta2f,dtheta1f,dtheta2f,tf)

% Robust control design for 2-D planar arm.
% input: initial and final state.
% output: Demostrate the performance of robust controller with parameter
% uncertainty.
% the nominal model parameter:
 % parameters in the paper.
% the nominal parameter vector b0 is

%% Observation of the perfromance
% The trajectories converges when the time is increased and the gamma
% values are set and the gain is increased. 

theta10=-0.5;
dtheta10 =0; 
theta1f = 0.8;
dtheta1f=0;
tf=80;

% plan a trajectory to reach target postion given by theta1f, dot theta1f,
% theta2f, dot theta2f.
theta20=-1;
dtheta20= 0.1; 
theta2f = 0.5;
dtheta2f=0;

%% Trajectory planning block
% Initial condition (TODO: CHANGE DIFFERENT INITIAL AND FINAL STATES)
x0=[-0.5,-1,0,0.1];
x0e = [-0.7,0.5,-0.2,0]; % an error in the initial state.
xf=[0.8,0.5, 0, 0];
% The parameter for planned joint trajectory 1 and 2.
global a1 a2 % two polynomial trajectory for the robot joint
nofigure=1;
% Traj generation.
a1 = planarArmTraj(theta10,dtheta10, theta1f, dtheta1f,tf, nofigure);
a2 = planarArmTraj(theta20,dtheta20, theta2f, dtheta2f,tf, nofigure);


torque=[];
torque1=[];
options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4, 1e-4, 1e-4]);
%% IMPLEMENT THE CONTROLLER
[T,X] = ode45(@(t,x)planarArmODERobust(t,x),[0 tf],x0e,options);

%% IMPLEMENT THE CONTROLLER TO AVOID CHATTERING.
[T,X] = ode45(@(t,x)planarArmODERobustApprx(t,x),[0 tf],x0e,options);

figure('Name','theta1');
plot(T, X(:,1),'r-');
hold on
plot(T, a1(1)+a1(2)*T+ a1(3)*T.^2+a1(4)*T.^3,'b-');
title('Theta_1 under Robust Control');


figure('Name','theta2');
plot(T, X(:,2),'r-');
hold on
plot(T, a2(1)+a2(2)*T+ a2(3)*T.^2+a2(4)*T.^3, 'b-');
title('Theta_2 under Robust Control');


figure('Name','torque with chatter');
plot(T,torque(1,1:size(T,1)), 'r-')
hold on
plot(T,torque(2,1:size(T,1)), 'b-')
title('Torque with chattering')

figure('Name','torque w/o chatter');
plot(T,torque1(1,1:size(T,1)), 'r-')
hold on
plot(T,torque1(2,1:size(T,1)), 'b-')
title('Torque w/o chattering ')
% PLOT THE INPUT TRAJECTORY

 
    

    function [dx ] = planarArmODERobust(t,x)
        %Todo: Select your feedback gain matrix Kp and Kd.
  
        % Compute the desired state and their time derivatives from planned
        % trajectory.
        
        m1 =10; m2=5; l1=1; l2=1; r1=0.5; r2 =.5; I1=10/12; I2=5/12;
        b0 = [ m1* r1^2 + m2*l1^2 + I1; m2*r2^2 + I2; m2*l1*r2];
        vec_t = [1; t; t^2; t^3]; % cubic polynomials
        theta_d= [a1'*vec_t; a2'*vec_t];
        %ref = [ref,theta_d];
        % compute the velocity and acceleration in both theta 1 and theta2.
        a1_vel = [a1(2), 2*a1(3), 3*a1(4), 0];
        a1_acc = [2*a1(3), 6*a1(4),0,0 ];
        a2_vel = [a2(2), 2*a2(3), 3*a2(4), 0];
        a2_acc = [2*a2(3), 6*a2(4),0,0 ];
        dtheta_d =[a1_vel*vec_t; a2_vel* vec_t];
        ddtheta_d =[a1_acc*vec_t; a2_acc* vec_t];
        theta= x(1:2,1);
        dtheta= x(3:4,1);
        
        
        %the true model
        m2t = m2+ 10*rand(1);% m1 true value is in [m1, m1+epsilon_m1] and epsilon_m1 a random number in [0,10];
        r2t = r2 + 0.5*rand(1);
        I2t = I2 + (15/12)*rand(1);
        
        a = I1+I2+m1*r1^2+ m2t*(l1^2+ r2t^2);
        b = m2t*l1*r2t;
        d = I2t+ m2t*r2t^2;
        % the actual dynamic model of the system is characterized by M and
        % C
        Mmat = [a+2*b*cos(x(2)), d+b*cos(x(2));  d+b*cos(x(2)), d];
        Cmat = [-b*sin(x(2))*x(4), -b*sin(x(2))*(x(3)+x(4)); b*sin(x(2))*x(3),0];
        invM = inv(Mmat);
        invMC = invM*Cmat;
        % For Ml lower bound
        m2t_l = m2; % m1 true value is in [m1, m1+epsilon_m1] and epsilon_m1 a random number in [0,10];
        r2t_l = r2 ;
        I2t_l = I2;
        
        a_l = I1+I2+m1*r1^2+ m2t_l*(l1^2+ r2t_l^2);
        b_l = m2t_l*l1*r2t_l;
        d_l = I2t_l+ m2t_l*r2t_l^2;
        
        M_l = [a_l+2*b_l*cos(x(2)), d_l+b_l*cos(x(2));  d_l+b_l*cos(x(2)), d_l];
        
        % For Mu upper bound
        m2t_u = m2 + 10;% m1 true value is in [m1, m1+epsilon_m1] and epsilon_m1 a random number in [0,10];
        r2t_u = r2 + 0.5;
        I2t_u = I2 + (15/12);
        
        a_u = I1+I2+m1*r1^2+ m2t_u*(l1^2+ r2t_u^2);
        b_u = m2t_u*l1*r2t_u;
        d_u = I2t_u+ m2t_u*r2t_u^2;
        
        M_u = [a_u+2*b_u*cos(x(2)), d_u+b_u*cos(x(2));  d_u+b_u*cos(x(2)), d_u];
        %invMLU = inv(M_l+M_u);
        M_bar = inv((M_l+M_u)/2);
        d_bar = M_bar(2,2);
        b_bar = M_bar(2,1) - d_bar /cos(x(2));
        C_bar = [-b_bar*sin(x(2))*x(4), -b_bar*sin(x(2))*(x(3)+x(4)); b_bar*sin(x(2))*x(3),0];
        invM_bar = inv(M_bar);
        invMC_bar = invM_bar*C_bar;
        
        Kp = 150*eye(2);
        Kd = 150*eye(2);
        gamma1 = 0.6;
        gamma2 = 0.1;
        gamma3 = 0.3;
        
        % LINK 1
        e1 = theta_d(1,1) - theta(1,1); % position error
        de1 = dtheta_d(1,1) - dtheta(1,1); % velocity error
        H1  = [e1;de1];
        
        P = eye(2);
        B = [0 ; 1];
        w1 = H1'*P*B;

        rho1 = gamma1*norm(H1) + gamma2*(norm(H1)*norm(H1)) + gamma3;
        if (w1 ~= 0)
            v1 = -w1*rho1/norm(w1);
        end
        if (w1 ==0)
            v1 = 0;
        end
        
        % LINK 2
        e2 = theta_d(2,1) - theta(2,1); % position error
        de2 = dtheta_d(2,1) - dtheta(2,1); % velocity error
        H2  = [e2;de2];
        
        P = eye(2);
        B = [0 ; 1];
        w2 = H2'*P*B;
                
        rho2 = gamma1*norm(H2) + gamma2*(norm(H2)*norm(H2))+ gamma3;
        if (w2 ~= 0)
            v2 = -w2*rho2/norm(w2);
        end
        if (w2 ==0)
            v2 = 0;
        end
        e =[e1;e2];
        de =[de1;de2];
        v =[v1;v2];
        
        aq = ddtheta_d + Kp*e + Kd*de + v ;
        tau = M_bar*aq + C_bar*theta_d ;
        torque=[torque,tau];
        %compute the robust controller
        %update the system state, compute dx
        dx=zeros(4,1);
        dx(1) = x(3);
        dx(2) = x(4);
        dx(3:4) = -invMC* x(3:4) +invM*tau; % because ddot theta = -M^{-1}(C \dot Theta) + M^{-1} tau
    end



    function [dx ] = planarArmODERobustApprx(t,x)
        %Todo: Select your feedback gain matrix Kp and Kd.
        m1 =10; m2=5; l1=1; l2=1; r1=0.5; r2 =.5; I1=10/12; I2=5/12;
        b0 = [ m1* r1^2 + m2*l1^2 + I1; m2*r2^2 + I2; m2*l1*r2];
        
        % Compute the desired state and their time derivatives from planned
        % trajectory.
        vec_t = [1; t; t^2; t^3]; % cubic polynomials
        theta_d= [a1'*vec_t; a2'*vec_t];
        %ref = [ref,theta_d];
        % compute the velocity and acceleration in both theta 1 and theta2.
        a1_vel = [a1(2), 2*a1(3), 3*a1(4), 0];
        a1_acc = [2*a1(3), 6*a1(4),0,0 ];
        a2_vel = [a2(2), 2*a2(3), 3*a2(4), 0];
        a2_acc = [2*a2(3), 6*a2(4),0,0 ];
        dtheta_d =[a1_vel*vec_t; a2_vel* vec_t];
        ddtheta_d =[a1_acc*vec_t; a2_acc* vec_t];
        theta= x(1:2,1);
        dtheta= x(3:4,1);
        
        %the true model
        m2t = m2+ 10*rand(1);% m1 true value is in [m1, m1+epsilon_m1] and epsilon_m1 a random number in [0,10];
        r2t = r2 + 0.5*rand(1);
        I2t = I2 + (15/12)*rand(1);
        
        a = I1+I2+m1*r1^2+ m2t*(l1^2+ r2t^2);
        b = m2t*l1*r2t;
        d = I2t+ m2t*r2t^2;
        % the actual dynamic model of the system is characterized by M and
        % C
        Mmat = [a+2*b*cos(x(2)), d+b*cos(x(2));  d+b*cos(x(2)), d];
        Cmat = [-b*sin(x(2))*x(4), -b*sin(x(2))*(x(3)+x(4)); b*sin(x(2))*x(3),0];
        invM = inv(Mmat);
        invMC = invM*Cmat;
        
        % For Ml lower bound
        m2t_l = m2; % m1 true value is in [m1, m1+epsilon_m1] and epsilon_m1 a random number in [0,10];
        r2t_l = r2 ;
        I2t_l = I2;
        
        a_l = I1+I2+m1*r1^2+ m2t_l*(l1^2+ r2t_l^2);
        b_l = m2t_l*l1*r2t_l;
        d_l = I2t_l+ m2t_l*r2t_l^2;
        
        M_l = [a_l+2*b_l*cos(x(2)), d_l+b_l*cos(x(2));  d_l+b_l*cos(x(2)), d_l];
        
        % For Mu upper bound
        m2t_u = m2 + 10;% m1 true value is in [m1, m1+epsilon_m1] and epsilon_m1 a random number in [0,10];
        r2t_u = r2 + 0.5;
        I2t_u = I2 + (15/12);
        
        a_u = I1+I2+m1*r1^2+ m2t_u*(l1^2+ r2t_u^2);
        b_u = m2t_u*l1*r2t_u;
        d_u = I2t_u+ m2t_u*r2t_u^2;
        
        M_u = [a_u+2*b_u*cos(x(2)), d_u+b_u*cos(x(2));  d_u+b_u*cos(x(2)), d_u];
        %invMLU = inv(M_l+M_u);
        M_bar = inv((M_l+M_u)/2);
        d_bar = M_bar(2,2);
        b_bar = M_bar(2,1) - d_bar /cos(x(2));
        C_bar = [-b_bar*sin(x(2))*x(4), -b_bar*sin(x(2))*(x(3)+x(4)); b_bar*sin(x(2))*x(3),0];
        invM_bar = inv(M_bar);
        invMC_bar = invM_bar*C_bar;
        
        Kp = 150*eye(2);
        Kd = 150*eye(2);
        gamma1 = 0.5;
        gamma2 = 0.1;
        gamma3 = 0.3;
        
        epsilon = 0.2;
        % LINK 1
        e1 = theta_d(1,1) - theta(1,1); % position error
        de1 = dtheta_d(1,1) - dtheta(1,1); % velocity error
        H1  = [e1;de1];
        
        P = eye(2);
        B = [0 ; 1];
        w1 = B'*P*H1;

        rho1 = gamma1*norm(H1) + gamma2*(norm(H1)*norm(H1)) + gamma3;
        if (norm(w1) < epsilon)
            v1 = -w1*rho1/norm(w1);
        end
        if (norm(w1) >= epsilon)
            v1 = -w1*rho1/epsilon;
        end
        
        % LINK 2
        e2 = theta_d(2,1) - theta(2,1); % position error
        de2 = dtheta_d(2,1) - dtheta(2,1); % velocity error
        H2  = [e2;de2];
        
        P = eye(2);
        B = [0 ; 1];
        w2 = B'*P*H2;

        rho2 = gamma1*norm(H2) + gamma2*(norm(H2)*norm(H2)) + gamma3;
        if (norm(w2) < epsilon)
            v2 = -w2*rho2/norm(w2);
        end
        if (norm(w2) >= epsilon)
            v2 = -w2*rho2/epsilon;
        end
        e =[e1;e2];
        de =[de1;de2];
        v =[v1;v2];
        
        aq = ddtheta_d + Kp*e + Kd*de + v ;
        tau1 = M_bar*aq + C_bar*theta_d ;
        
        torque1=[torque1, tau1];
        
        
        %compute the robust controller to avoid chattering.
        %update the system state, compute dx
        dx=zeros(4,1);
        dx(1) = x(3);
        dx(2) = x(4);
        dx(3:4) = -invMC* x(3:4) +invM*tau1; % because ddot theta = -M^{-1}(C \dot Theta) + M^{-1} tau
    end     
    %% Trajectory planning using polynomial functions.
    function [a] = planarArmTraj(theta10,dtheta10, theta1f, dtheta1f,tf, nofigure)
        % Input: Initial and final position and velocities, planning horizon [0,tf]
        % nofigure=1 then do not output the planned trajectory.
        % Cubic polynomial trajectory.

        % formulate the linear equation and solve.
        M= [1 0 0 0;
            0 1 0 0;
            1 tf tf^2 tf^3;
            0 1 2*tf 3*tf^2];
        b=[theta10; dtheta10;theta1f; dtheta1f];
        a=M\b;
        t=0:0.01:tf;

        if nofigure==1
            return
        else

        figure('Name','Position (degree)');
        plot(t,a(1)+a(2)*t+ a(3)*t.^2+a(4)*t.^3,'LineWidth',3);
        title('Position (degree)')
        grid

        figure('Name','Velocity (degree/s)');
        plot(t,a(2)*t+ 2*a(3)*t +3*a(4)*t.^2,'LineWidth',3);
        title('Velocity (degree/s)')
        grid

        figure('Name','Acceleration (degree/s^2)');
        plot(t, 2*a(3) +6*a(4)*t,'LineWidth',3);
        title('Acceleration (degree/s^2)')
        grid
    end
end

end