% Abstract class for terminal and running costs
classdef (Abstract) Cost
    methods
        % Terminal state cost
        phi = phi(obj, x_f, x_bar);
        
        % Terminal state cost derivatives
        phi_x = phi_x(obj, x_f, x_star);
        phi_xx = phi_xx(obj, x_f, x_star);
        
        % Running cost
        L = L(obj, x, u, dt);
        
        % Running cost derivatives
        L_x = L_x(obj, x, u, dt);
        L_u = L_u(obj, x, u, dt);
        L_xx = L_xx(obj, x, u, dt);
        L_uu = L_uu(obj, x, u, dt);
        L_xu = L_xu(obj, x, u, dt);
        L_ux = L_ux(obj, x, u, dt);
    end
end

