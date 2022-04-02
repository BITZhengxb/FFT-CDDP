% Copyright (C) 2021 Xiaobo Zheng - All Rights Reserved
% developed by Xiaobo at BIT Jan 2021

% Abstract class for terminal and running costs
classdef (Abstract) Constraints
    methods
        % Constraints function
        G = G(obj, x, u, dt);
        
        % Constraints function derivatives
        G_x = G_x(obj, x, u, dt);
        G_u = G_u(obj, x, u, dt);
    end
end

