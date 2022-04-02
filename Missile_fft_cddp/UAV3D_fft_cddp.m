% Copyright (C) 2021 Xiaobo Zheng - All Rights Reserved
% developed by Xiaobo at BIT Jan 2021

%% Flexible Final time Constrained Differential Dynamic Programming Algorithm Function

% Performs Flexible Final time Constrained Differential Dynamic Programming
% Algorithm to produce a locally optimal control sequence.
%
% Inputs
% 
% x_0           : initial state vector
% x_star        : target state vector
% t_f           : time horizon
% dt            : time interval
% dyn           : instance of Dynamics class
% cost          : instance of Cost class
% constraints	: instance of Constraints class
% u_max         : maximum magnitude of control, used for clamping
% num_iter      : number of DDP iterations
% alpha         : DDP learning parameter for line searching
%
% Outputs
% 
% sol : structure with solution components
%
function sol = UAV3D_fft_cddp(x_0, x_star, t_f, dt, dyn, cost, constraints, u_max, num_iter, alpha)
    %% Allocate arrays for DDP

    % Number of states along trajectory
    N = floor(t_f ./ dt);
    
    % Time stamp array and timestep
    t = linspace(0.0, t_f - dt, N);
    t_max = t_f + 20;
    t_min = max(t_f-20, 5);
    
    % Cost history
    J = zeros(num_iter, 1);
    
    % Control energy history
    E = zeros(num_iter, 1);
    
    % Terminal time history
    Tf = zeros(num_iter, 1);
    
    % State trajectory
    x = cell(N, 1);
    x_new = cell(N, 1);
    x{1} = x_0;
    x_new{1} = x_0;
    
    % Control Input trajectory
    u = cell(N, 1);
    
    % Value function and derivatives
    V = zeros(N, 1);
    V_x = cell(N, 1);
    V_xx = cell(N, 1);
    V_t = cell(1, 1);
    V_tx = cell(1, 1);
    V_tt = cell(1, 1);
    V_p = cell(1, 1);
    V_pp = cell(1, 1);
    
    % State action value derivatives
    Q_x = cell(N, 1);
    Q_u = cell(N, 1);
    Q_xx = cell(N, 1);
    Q_uu = cell(N, 1);
    Q_xu = cell(N, 1);
    Q_ux = cell(N, 1);
    
    % Constraints derivatives
    Dx = cell(N, length(constraints));
    lambda = cell(N, length(constraints));

    %% Initialize DDP with a random input sequence

    fprintf("initializing with a random control sequence...\n")
    
    % Generate random control sequence
    for k = 1:N-1
        u{k} = zeros(numel(u_max), 1);
    end
    u{N} = zeros(numel(u_max), 1);
    
%     [path, u_PN] = Proportional_guidance3D(x_0, x_star, dyn, t_f, dt);
%     for k = 1:N
%         u{k}(1) = u_PN(1, k);
%         u{k}(2) = u_PN(2, k);
%         u{k} = u{k}';
%     end
    
    fprintf("generating initial trajectory...\n")
    
    % Generate initial trajectory using random control sequence
    for k = 1:N-1
        x_new{k+1} = x_new{k} + dyn.F(x_new{k}, u{k}) .* dt;
    end
    
    %% Perform main DDP iterations on the trajectory and input sequence
    
    fprintf("beginning DDP...\n")
    
    sigma = 1;
    for k = 1:N-1
        for constr = 1:length(constraints)
            lambda{k,constr} = 0.01;
        end
    end
    
    for i = 1:num_iter
%         fprintf("DDP iteration %d out of %d...\n", i, num_iter);
        
        % Update control sequence from previous iteration
        if i > 1
            for k = 1:N-1
                % Compute control update feed-forward and feed-back
                du_ff = -inv(Q_uu{k}) * Q_u{k};
                du_fb = -inv(Q_uu{k}) * Q_ux{k} * (x_new{k} - x{k});
                
                % Limit feed forward control modification with clamping
                du = alpha .* (du_ff + du_fb);
                for m = 1:numel(u_max)
                    du(m) = min(u_max(m), max(-u_max(m), du(m) + u{k}(m))) - u{k}(m);
                end
                
                % Update control
                u{k} = u{k} + du;
                
                % Compute next state in trajectory with new control
                x_new{k+1} = x_new{k} + dyn.F(x_new{k}, u{k}) .* dt;
                
                % Return error if problem with trajectory
                if isnan(x_new{k+1})
                    sol = assemble_solution(x, u, t, J, E, Tf, Q_u, Q_uu, Q_ux, 1, i);
                    fprintf("Problem with trajectory ...\n")
                    return
                end
            end
            
            % Update terminal time
            % Compute terminal time update feed-forward and feed-back
            dtf_ff = -inv(V_tt) * V_t;
            dtf_fb = -inv(V_tt) * V_tx * (x_new{N} - x{N});

            % Update terminal time
            dtf = 0.4 * ( dtf_ff + dtf_fb );
            dtf = min(t_max, max(t_min, dtf + t_f)) - t_f;
            
            % Update time stamp array and timestep
            dN = fix( dtf / dt );
            N_old = N;
            N = N + dN;
            
            if dN >= 0
                % Update state trajectory
                x = [x; cell(dN, 1)];
                x_new = [x_new; cell(dN, 1)];

                % Update control Input trajectory
                u = [u; cell(dN, 1)];

                % Update value function and derivatives
                V = [V; zeros(dN, 1)];
                V_x = [V_x; cell(dN, 1)];
                V_xx = [V_xx; cell(dN, 1)];

                % Update state action value derivatives
                Q_x = [Q_x; cell(dN, 1)];
                Q_u = [Q_u; cell(dN, 1)];
                Q_xx = [Q_xx; cell(dN, 1)];
                Q_uu = [Q_uu; cell(dN, 1)];
                Q_xu = [Q_xu; cell(dN, 1)];
                Q_ux = [Q_ux; cell(dN, 1)];
                
                % Update constraints derivatives
                Dx = [Dx; cell(dN, length(constraints))];
                lambda = [lambda; cell(dN, length(constraints))];
                
                % Update control sequence
                for k = 1:dN-1
                    u{N_old + k} = zeros(numel(u_max), 1);
                end

                % Update trajectory using zero control sequence
                for k = 0:dN-1
                    x_new{N_old + k+1} = x_new{N_old + k} + dyn.F(x_new{N_old + k}, u{N_old + k}) .* dt;
                end
            else
                % Update state trajectory
                x = x(1:N);
                x_new = x_new(1:N);

                % Update control Input trajectory
                u = u(1:N);

                % Update value function and derivatives
                V = V(1:N);
                V_x = V_x(1:N);
                V_xx = V_xx(1:N);

                % Update state action value derivatives
                Q_x = Q_x(1:N);
                Q_u = Q_u(1:N);
                Q_xx = Q_xx(1:N);
                Q_uu = Q_uu(1:N);
                Q_xu = Q_xu(1:N);
                Q_ux = Q_ux(1:N);
                
                % Update constraints derivatives
                Dx = Dx(1:N, :);
                lambda = lambda(1:N, :);
                
            end
            
            % Return error if problem with trajectory
            if N <= 0
                sol = assemble_solution(x, u, t, J, E, Tf, Q_u, Q_uu, Q_ux, 1, i);
                fprintf("Problem with terminal time ...\n")
                return
            end
            
            t_f = N * dt;
            t = linspace(0.0, t_f - dt, N);
            
            Num_current = min(N_old, N);
            theta_old = zeros(Num_current, 1);
            theta_new = zeros(Num_current, 1);
            for k = 1:Num_current
                theta_old(k) = x{k}(5);
                theta_new(k) = x_new{k}(5);
            end
            if i > 1 && max(abs(theta_new - theta_old)) < 1e-4 % Q、R-0.1 1e-6; Q、R-1.0 1e-5
                fprintf("Convergence theta first !!! ...\n")
                u{N} = u{N-1};
                break;
            end
            
        end
        
        % Update the current trajectory
        x = x_new;
        
        % Compute total cost
%         J(i) = cost.phi(x{N}, x_star);
        J(i) = cost.phi(x{N}, x_star) + t_f;
        for k = 1:N-1
            J(i) = J(i) + cost.L(x{k}, x_star, u{k}, dt);
        end
        
        % Compute control energy usage
        for k = 1:N-1
            E(i) = E(i) + 0.5 .* u{k}.' * u{k} .* dt;
        end
        
        % Log terminal time
        Tf(i) = t_f;
        
        if i > 1 && abs(J(i) - J(i-1)) < 1e-6 % Q、R-0.1 1e-6; Q、R-1.0 1e-5
            fprintf("Convergence !!! ...\n")
            u{N} = u{N-1};
            break;
        end
        
        % Compute terminal value function and derivatives
%         V(N) = cost.phi(x{N}, x_star);
        V(N) = cost.phi(x{N}, x_star) + t_f;
        V_x{N} = cost.phi_x(x{N}, x_star);
        V_xx{N} = cost.phi_xx(x{N}, x_star);
        u{N} = u{N-1};
%         V_t = V_x{N}' * dyn.F(x{N}, u{N});
        V_t = V_x{N}' * dyn.F(x{N}, u{N}) + 1 * dt;
        V_tx = ( V_xx{N} * dyn.F(x{N}, u{N}) )' + V_x{N}' * ( dyn.Phi(x{N}, u{N}, dt) - eye(numel(x{N})) ) ./ dt;
        V_tt = V_tx * dyn.F(x{N}, u{N});
        
        % Perform backwards pass
        for k = N-1:-1:1
            % Compute state-action value function derivatives
            if k == N-1
                V_p = V_x{N}' - V_t' * (V_tt \ V_tx);
                V_pp = V_xx{N} - V_tx' * (V_tt \ V_tx); % V_tx{N}' * inv(V_tt{N}) * V_tx{N}
                Q_x{N-1} = cost.L_x(x{N-1}, x_star, u{N-1}, dt) + dyn.Phi(x{N-1}, u{N-1}, dt)' * V_p';
                Q_u{N-1} = cost.L_u(x{N-1}, u{N-1}, dt) + dyn.beta(x{N-1}, u{N-1}, dt)' * V_p';
                Q_xx{N-1} = cost.L_xx(x{N-1}, x_star, u{N-1}, dt) + dyn.Phi(x{N-1}, u{N-1}, dt)' * V_pp * dyn.Phi(x{N-1}, u{N-1}, dt);
                Q_uu{N-1} = cost.L_uu(x{N-1}, u{N-1}, dt) + dyn.beta(x{N-1}, u{N-1}, dt)' * V_pp * dyn.beta(x{N-1}, u{N-1}, dt);
                Q_xu{N-1} = cost.L_xu(x{N-1}, u{N-1}, dt) + dyn.Phi(x{N-1}, u{N-1}, dt)' * V_pp * dyn.beta(x{N-1}, u{N-1}, dt);
                Q_ux{N-1} = cost.L_ux(x{N-1}, u{N-1}, dt) + dyn.beta(x{N-1}, u{N-1}, dt)' * V_pp * dyn.Phi(x{N-1}, u{N-1}, dt);
            else
                Q_x{k} = cost.L_x(x{k}, x_star, u{k}, dt) + dyn.Phi(x{k}, u{k}, dt)' * V_x{k+1};
                Q_u{k} = cost.L_u(x{k}, u{k}, dt) + dyn.beta(x{k}, u{k}, dt)' * V_x{k+1};
                Q_xx{k} = cost.L_xx(x{k}, x_star, u{k}, dt) + dyn.Phi(x{k}, u{k}, dt)' * V_xx{k+1} * dyn.Phi(x{k}, u{k}, dt);
                Q_uu{k} = cost.L_uu(x{k}, u{k}, dt) + dyn.beta(x{k}, u{k}, dt)' * V_xx{k+1} * dyn.beta(x{k}, u{k}, dt);
                Q_xu{k} = cost.L_xu(x{k}, u{k}, dt) + dyn.Phi(x{k}, u{k}, dt)' * V_xx{k+1} * dyn.beta(x{k}, u{k}, dt);
                Q_ux{k} = cost.L_ux(x{k}, u{k}, dt) + dyn.beta(x{k}, u{k}, dt)' * V_xx{k+1} * dyn.Phi(x{k}, u{k}, dt);
            end
            
            % construct equality quadratic programming problem
            for constr = 1:length(constraints)
                G_value =  constraints(constr).G(x{k}, u{k}, dt) ;
                if G_value < 0.00001
                    Dx{k,constr} = constraints(constr).G_x(x{k}, u{k}, dt);
%                     lambda{k,constr} = lambda{k,constr} + sigma * G_value;
                    Q_x{k} = Q_x{k} + 0.004 * Dx{k,constr}' + G_value * sigma * Dx{k,constr}';
                    Q_xx{k} = Q_xx{k} + Dx{k,constr}'* sigma *Dx{k,constr};
                end
            end
            
            % Compute the value function derivatives
            V_x{k} = Q_x{k} - Q_xu{k} * (Q_uu{k} \ Q_u{k});
            V_xx{k} = Q_xx{k} - Q_xu{k} * (Q_uu{k} \ Q_ux{k});
            
        end
        
        sigma = 1.10 * sigma;
    end
    
    x{end+1} = x{end} + dyn.F(x{end}, u{end}) .* dt;
    t = [t t_f];
    
    %% Assemble and return solution structure
    fprintf("DDP iteration %d .\n", i);
    fprintf("finished DDP, assembling results for post-processing...\n");
    
    % Assemble solution
    sol = assemble_solution(x, u, t, J, E, Tf, Q_u, Q_uu, Q_ux, 0, i);
end

%% Helper Functions for DDP Algorithm

% Assembles and returns solution structure
%
% Inputs
%
% x               : locally optimal state trajectory
% u               : locally optimal control sequence
% t               : discretized time stamps
% J               : iteration history of cost function
% E               : iteration history of control energy
% Q_u, Q_uu, Q_ux : derivatives of state action value function
% error           : zero if no error, nonzero else
% iter            : iterations
%
% Outputs
%
% sol : solution structure
function sol = assemble_solution(x, u, t, J, E, Tf, Q_u, Q_uu, Q_ux, error, iter)

    % Solution structure
    sol = struct;
    sol.error = error;
    sol.x = x;
    sol.u = u;
    sol.t = t;
    sol.dt = t(2) - t(1);
    sol.J = J;
    sol.E = E;
    sol.Tf = Tf;
    sol.Q_u = Q_u;
    sol.Q_uu = Q_uu;
    sol.Q_ux = Q_ux;
    sol.iter = iter;
    
    return
end