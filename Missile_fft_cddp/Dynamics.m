% Abstract class for system dynamics
classdef (Abstract) Dynamics
    methods
        % Equations of motion in state space representation
        dxdt = F(obj, x, u)
        
        % Linearized equations of motion
        Phi = Phi(obj, x, u, dt);
        beta = beta(obj, x, u, dt);
    end
end

