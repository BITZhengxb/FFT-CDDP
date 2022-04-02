% Copyright (C) 2021 Xiaobo Zheng - All Rights Reserved
% developed by Xiaobo at BIT Jan 2021

% Class for dynamics of inverted pendulum system that inherits from 
% abstract Dynamics class.
classdef UAV3D_Dynamics < Dynamics
    properties
        m; % mass of vehicle
        g; % 
        s_p;
    end
    
    methods
        % Constructor
        function obj = UAV3D_Dynamics(m, g, scale_para)
            obj.m = m;
            obj.g = g;
            obj.s_p = scale_para;
        end
        
        % Equations of motion in state space representation
        function dxdt = F(obj, x, u)
            rho = 1.15579 - 1.058*10^(-4)*x(3)*obj.s_p + 3.725*10^(-9)*(x(3)*obj.s_p)^2 - 6*10^(-14)*(x(3)*obj.s_p)^3; % air density model
            S = 0.0324;
            Cd = 0.0169;
            D = 0.5 .* rho * x(4)^2 * S * Cd;
            
            x_dot = x(4) .* cos(x(5)).* cos(x(6)) / obj.s_p;
            y_dot = x(4) .* cos(x(5)).* sin(x(6)) / obj.s_p;
            z_dot = x(4) .* sin(x(5)) / obj.s_p;
            V_dot = -D ./ obj.m - obj.g .* sin(x(5));
            gama_dot = ( -u(2) - obj.g .* cos(x(5)) ) ./ x(4);
            psi_dot = u(1) ./ ( x(4) .* cos(x(5)) );
                                       
            dxdt = [x_dot; y_dot; z_dot; V_dot; gama_dot; psi_dot];
        end
        
        % Linearized equations of motion (determined using symbolic math)
        function Phi = Phi(obj, x, u, dt)
            rho = 1.15579 - 1.058*10^(-4)*x(3)*obj.s_p + 3.725*10^(-9)*(x(3)*obj.s_p)^2 - 6*10^(-14)*(x(3)*obj.s_p)^3; % air density model
            S = 0.0324;
            Cd = 0.0169;
            
            F_x = zeros(numel(x));
            
            F_x(1,1) = 0.0;
            F_x(1,2) = 0.0;
            F_x(1,3) = 0.0;
            F_x(1,4) = cos(x(5)) .* cos(x(6)) / obj.s_p;
            F_x(1,5) = -x(4) .* sin(x(5)) .* cos(x(6)) / obj.s_p;
            F_x(1,6) = -x(4) .* cos(x(5)) .* sin(x(6)) / obj.s_p;
            
            F_x(2,1) = 0.0;
            F_x(2,2) = 0.0;
            F_x(2,3) = 0.0;
            F_x(2,4) = cos(x(5)) .* sin(x(6)) / obj.s_p;
            F_x(2,5) = -x(4) .* sin(x(5)) .* sin(x(6)) / obj.s_p;
            F_x(2,6) = x(4) .* cos(x(5)) .* cos(x(6)) / obj.s_p;
            
            F_x(3,1) = 0.0;
            F_x(3,2) = 0.0;
            F_x(3,3) = 0.0;
            F_x(3,4) = sin(x(5)) / obj.s_p;
            F_x(3,5) = x(4) .* cos(x(5)) / obj.s_p;
            F_x(3,6) = 0.0;
            
            F_x(4,1) = 0.0;
            F_x(4,2) = 0.0;
            F_x(4,3) = 0.0;
            F_x(4,4) = -rho * x(4) * S * Cd / obj.m;
            F_x(4,5) = -obj.g .* cos(x(5));
            F_x(4,6) = 0.0;
            
            F_x(5,1) = 0.0;
            F_x(5,2) = 0.0;
            F_x(5,3) = 0.0;
            F_x(5,4) = ( u(2) + obj.g .* cos(x(5)) ) ./ (x(4)^2);
            F_x(5,5) = ( obj.g .* sin(x(5)) ) ./ (x(4));
            F_x(5,6) = 0.0;
            
            F_x(6,1) = 0.0;
            F_x(6,2) = 0.0;
            F_x(6,3) = 0.0;
            F_x(6,4) = ( -u(1) ) ./ ( cos(x(5)) .* x(4)^2 );
            F_x(6,5) = ( u(1) ) .* sin(x(5)) ./ ( cos(x(5))^2 .* x(4) );
            F_x(6,6) = 0.0;
            
            Phi = eye(numel(x)) + F_x .* dt;
        end
        function beta = beta(obj, x, u, dt)
            F_u = zeros(numel(x), numel(u));
            
            F_u(1,1) = 0.0;
            F_u(1,2) = 0.0;
            
            F_u(2,1) = 0.0;
            F_u(2,2) = 0.0;
            
            F_u(3,1) = 0.0;
            F_u(3,2) = 0.0;
            
            F_u(4,1) = 0.0;
            F_u(4,2) = 0.0;
            
            F_u(5,1) = 0.0;
            F_u(5,2) = -1 ./ x(4);
            
            F_u(6,1) = 1 ./ ( x(4) .* cos(x(5)) );
            F_u(6,2) = 0.0;
                 
            beta = F_u .* dt;
        end
    end
end

