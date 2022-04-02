% Copyright (C) 2021 Xiaobo Zheng - All Rights Reserved
% developed by Xiaobo at BIT Jan 2021

% Class for obstacle
classdef CircleConstraints_3D < Constraints
    properties
        Center; % center of the obstacle
        r; % radius of the obstacle
    end
    
    methods
        % Constructor
        function obj = CircleConstraints_3D(Center, r)
            obj.Center = Center;
            obj.r = r;
        end
        
        % Circle constraints
        function G = G(obj, x, ~, dt)
            G = ( (x(1:2, 1)' - obj.Center) * (x(1:2, 1)' - obj.Center)' - obj.r^2 ) .*dt ;
        end
        
        % Circle constraints derivatives
        function G_x = G_x(obj, x, ~, dt)
            G_x = [2 .* (x(1:2, 1)' - obj.Center), 0, 0, 0, 0] .* dt;
        end
        function G_u = G_u(~, ~, u, ~)
            G_u = zeros(1, numel(u));
        end
    end
end

