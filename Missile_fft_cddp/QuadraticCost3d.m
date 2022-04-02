% Class for quadratic terminal and running costs that inherits from
% abstract Cost class
classdef QuadraticCost3d < Cost
    properties
        Q_f; % terminal cost weights
        Q_s; % terminal cost weights
        R; % control cost weights
    end
    
    methods
        % Constructor
        function obj = QuadraticCost3d(Q_f, Q_s, R)
            obj.Q_f = Q_f;
            obj.Q_s = Q_s;
            obj.R = R;
        end
        
        % Terminal state cost
        function phi = phi(obj, x_f, x_bar)
            phi = 0.5 .* (x_f - x_bar).' * obj.Q_f * (x_f - x_bar);
        end
        
        % Terminal state cost derivatives
        function phi_x = phi_x(obj, x_f, x_star)
            phi_x = obj.Q_f * (x_f - x_star);
        end
        function phi_xx = phi_xx(obj, ~, ~)
            phi_xx = obj.Q_f;
        end
        
        % Running cost
        function L = L(obj, x, x_star, u, dt)
            L = (0.5 .* (x - x_star).' * obj.Q_s * (x - x_star) + 0.5 .* u.' * obj.R * u) .* dt;
        end
        
        % Running cost derivatives
        function L_x = L_x(obj, x, x_star, ~, dt)
            L_x = (obj.Q_s * (x - x_star)) .* dt;
        end
        function L_u = L_u(obj, ~, u, dt)
            L_u = (obj.R * u) .* dt;
        end
        function L_xx = L_xx(obj, ~, ~, ~, dt)
            L_xx = obj.Q_s .* dt;
        end
        function L_uu = L_uu(obj, ~, ~, dt)
            L_uu = obj.R .* dt;
        end
        function L_xu = L_xu(~, x, u, ~)
            L_xu = zeros(numel(x), numel(u));
        end
        function L_ux = L_ux(~, x, u, ~)
            L_ux = zeros(numel(u), numel(x));
        end
    end
end

