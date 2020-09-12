clear all
rng(2)



%% Global System Setup

% How many links?
n = 3;
% NB: tested up to 20 links. Took around 2 minutes to solve for 10
% seconds. Y'all apprecitate real-time physics simulations more now?

% system parameters
mval  = 0.5*ones(1, n);  % mass of the link, kg
mcval = 1;  % mass of the cart, kg
lval  = 0.3*ones(1, n);  % length of the links, m
gval  = 9.81;  % gravity, m/s^2

% animation speed 
animspeed = 0.5; % 1: time-true

% write output gif?
gifoutput = 0;  % 1: yes, 0: no
gifname = "TripleInverted_Sidestep.gif";




%% Simulation-Case set up

% NB on the angles: (0 --> upright, pi --> down (CW rotation))

control_on = 1;
simcase = "sidestep";

if simcase=="sidestep"
    % initial condition
    s_init = 0; ds_init = 0;  % cart position and speed
    th_init = zeros(1, n); dth_init = zeros(1, n); % link angles and rot speed

    % target condition
    s_final = 2; ds_final = 0; % cart position and speed
    th_final = zeros(1, n); dth_final = zeros(1, n); % link angles and rot speed
    
elseif simcase=="disturbance"
    % initial condition
    s_init = 0; ds_init = 0;  % cart position and speed
    th_init = 0.075*rand(1, n); dth_init = zeros(1, n); % link angles and rot speed

    % target condition
    s_final = 0; ds_final = 0; % cart position and speed
    th_final = zeros(1, n); dth_final = zeros(1, n); % link angles and rot speed
    
elseif simcase=="custom"
    % initial condition
    s_init = 0; ds_init = 0;  % cart position and speed
    th_init = 0.5*pi*ones(1, n); dth_init = 0*ones(1, n); % link angles and rot speed

    % target condition
    s_final = 0; ds_final = 0; % cart position and speed
    th_final = pi*ones(1, n); dth_final = zeros(1, n); % link angles and rot speed
end

% time range to solve for
t_range = [0 7];  % in seconds






%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% start of automatic calculations %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% declare symbolic variables

% original states
syms s(t) % cart position
syms(sprintfc('th_%d(t)', 1:n))  % angles
eval("th = [" + join(sprintfc('th_%d', 1:n), ",") + "];")  % angle vector

% save "original" state variables in one vector
origvars = [s th];

% input placeholder
syms Fc

% masses, lengths of the links
m = sym('m_%d', [1 n]); l = sym('l_%d', [1 n]); 

% mass of the cart
mc = sym('mc');

% gravity
g = sym('g');

% inertias
I = 1/12 * m .* l.^2;




%% System Description

%%%%% Speeds and Positions %%%%%
% speed of the cart
ds = diff(s, t);

% link CoG positions as function of cart position and link angles 
% also, get its symbolic derivative
[q, ~] = positions(s, th, l, n);
dq = diff(q, t);

% rotational speed (symbolic derivative of the link angles)
dth = diff(th, t);



%%%%% Kinetic Energy %%%%%
% T is sum of translational and rotational energies
T = .5 * mc * ds^2 ... cart
    + .5 * sum(sum(m .* dq.^2)) ... translation of links
    + .5 * sum(I .* dth.^2); % rotation of links


%%%%% Potential Energy %%%%%
% V is the sum of the y-CoG heights of the links times mg
V = g * sum(m .* ([0 1] * q));


%%%%% (symbolic) Lagrangian of the System %%%%%
L = T - V;




%% Formulate Euler Lagrange Equations of Motion

% unfortunately, MATLAB cannot take partial symbolic derivatives with
% respect to a time derivative, as required by the E-L EoMs.
% 
% therefor, we have to substitute them by a place-holder symbolic variable,
% then take the derivatives, and then substitute back to the original

% Define these variables
syms ds_dummy(t)
syms(sprintfc('dth_dummy_%d(t)', 1:n));
eval("dth_dummy = [" + join(sprintfc('dth_dummy_%d', 1:n), ",") + "];")

% substitute
L = subs(L, ds, ds_dummy);  % cart
L = subs(L, dth, dth_dummy);  % links

% get all n+1 Euler Lagrange equations LHS
ELs  = diff(functionalDerivative(L, ds_dummy), t) ...
            - functionalDerivative(L, s);  % cart
ELth = diff(functionalDerivative(L, dth_dummy), t) ...
            - functionalDerivative(L, th);  % links

% substitute back to original variables
ELs  = subs(ELs, [ds_dummy, dth_dummy], [ds, dth]);
ELth = subs(ELth, [ds_dummy, dth_dummy], [ds, dth]);

% Set equal to external force
ODEs = [ELs; ELth] == [Fc; zeros(n, 1)];  % summarized as one system of ODEs




%% substitute numerical parameters

% substitute numerical parameters into the ODEs
ODEs = subs(ODEs, [m, mc, l, g], [mval, mcval, lval, gval]);




%% Formulate First Order ODEs of the open-loop system
% Currently, our ODEs system is of the following implicit form:
%       f(x, dx/dt ..., d^n x/dt^n, t) == 0
% 
% 1. We can the expand the states to arrive at a larger, non-linear system
% of exclusively first-order ODEs of the form:
%        f(y, dy/dt, t) == 0
% 
% 2. Furthermore, semi-linear means that the highest order derivative terms
% (after expansion of the states simply _the_ derivatives terms) are all
% linear which means the system can be solved for those terms to yield an
% ODE system of the form:
%       M * dy/dt = f(y, t)   where M is the square "mass matrix"
%
% AAAaand, we can integrate this using (for instance) MATLABs ode15s


% 1. Reduce the order by expanding the states (more variables and 
% equations)
[ODEs_1st, vars] = reduceDifferentialOrder(ODEs, origvars);

% 2. The system is still implicit. That means that the derivative
% terms are hidden somewhere inside the f_i(y, dy/dt, t) == 0
% Calculate the mass matrix "massM" and the modified right hand side "f"
[massM, f] = massMatrixForm(ODEs_1st, vars);




%% Linearization of the open-loop system

% write as one-sided ODE system
dvarsdt = massM \ f;

% change variables to x1, x2...
x = sym("x_%d", [2*(n+1) 1]);
dxdt = subs(dvarsdt, vars, x);

% symbolic jacobian function
J = jacobian(dxdt, [x; Fc]);

% numeric jacobian function
Jf_fun = matlabFunction(J(:, 1:end-1), 'vars', [x; Fc]);
Jg_fun = matlabFunction(J(:, end), 'vars', [x; Fc]);




%% Control strategy --> LQR

% input: force at cart

%%% inverted pendulum with n links:
% Optimal control with LQR

% penalty functions
Q = diag([1e2, 2e2 * ones(1, n), 1e1, 1 * ones(1, n)]);  % state penalty
R = 1e0;  % control effort penalty

% summarize target position in vector
target_pos = [s_final th_final ds_final dth_final]';

% get linearization around target position and Fc = 0
linpoint = [target_pos; 0];
c = num2cell(linpoint); % somehow we need it as a cell...

A = Jf_fun(c{:});  % system matrix using the jacobian expression above
B = Jg_fun(c{:});  % control matrix using the jacobian expression above

% solve optimal control problem
[P, K, ~] = icare(A,B,Q,R,[],[],[]); % solve the algebraic riccarti Eq

% gain the system
Fc_LQR = -K * (mod(pi+([s, th, diff([s th], t)]' - target_pos), 2*pi)-pi);

if control_on == 0
    Fc_LQR = 0;
end




%% Close the loop

% substitute our control expression into the ODEs for simulation
% this step can be seen as "closing the loop"
ODEs_controlled = subs(ODEs, Fc, Fc_LQR);

% reformulate the dynamics. This time closed loop
[ODEs_1st_controlled, vars] ...
    = reduceDifferentialOrder(ODEs_controlled, origvars);
[massM_controlled, f_controlled] ...
    = massMatrixForm(ODEs_1st_controlled, vars);




%% Solve the ODEs

% States and Derivatives are still symbolic. Convert to executable matlab
% function:
massM_num = odeFunction(massM_controlled, vars);
f_num = odeFunction(f_controlled, vars);

% summarize initial conditions in a vector
y0 = [s_init, th_init, ds_init, dth_init];

% set up option for the solver and solve (to a really small tolerance)
tic
    opt = odeset('Mass', massM_num, 'RelTol', 1e-5, 'AbsTol' , 1e-5);
    [tSol, ySol] = ode15s(f_num, t_range, y0, opt);
toc




%%
if false 
    % this section still needs work. It will (at some point, hopefully) be
    % able to solve a transitional problems (like swing-up or swing-down)
    % problem based on a pretty recent constrained and overdetermined BVP
    % approach developed by Knut Graichen and Michael Zeitz
    % https://link-springer-com.tudelft.idm.oclc.org/chapter/10.1007/11529798_15
    
    %% solve side-step BVP


    % substitute our control expression into the ODEs for simulation
    % this step can be seen as "closing the loop"
    syms Fc_bvp(t)
    ODEs_controlled = subs(ODEs, Fc, Fc_bvp);

    % ODEs_controlled = [ODEs_controlled; diff(Fc_bvp, t) == diff(Fc_bvp, t)];

    origvarsBVP = [origvars Fc_bvp];


    % reformulate the dynamics. This time closed loop
    [ODEs_1st_controlled, vars_BVP] ...
        = reduceDifferentialOrder(ODEs_controlled, origvarsBVP);
    [massM_BVP, f_BVP] ...
        = massMatrixForm(ODEs_1st_controlled, vars_BVP);

    tlambda = linspace(0, 2, 4);
    lambda_init = ones(size(tlambda));
    solinit = bvpinit([0 2], @bvpinit_fun, lambda_init);

    bvpopts = bvpset('AbsTol', 1e-9, 'RelTol', 1e-6, 'NMax', 2000, 'Stats', 'on');
    sBVP = bvp4c( ...
        @(t, x, lambda) bvp_ODEs(...
            massM_BVP(:, [true(1, n+1) false true(1, n+1)]), f_BVP,...
            vars_BVP, t, x, tlambda, lambda...
            ), ...
        @(xa, xb, lambda) bvp_bvs(xa, xb, lambda), ...
        solinit,...
        bvpopts);

    ySol = sBVP.y'; tSol = sBVP.x'; y0 = sBVP.y(:, 1)';

    bvfun(sBVP.y(:, 1), sBVP.y(:, end))

end




%% plotting

% positions of initial condition
[cgpos, endpoints] = positions(y0(1), y0(2:(n+1)), lval, n);

% set up figure window and plot initial condition
h = figure('Position', [50 50 1600 900]);
hold on
    plt_line = plot(endpoints(1,:), endpoints(2,:), ...
        'b-', 'LineWidth', 4);
    plt_cart = scatter(endpoints(1,1), endpoints(2,1), ...
        500, '+', 'k', 'LineWidth', 3);
    plt_cons = scatter(endpoints(1, 2:end), endpoints(2, 2:end), ...
        700, '.', 'b', 'LineWidth', 1);
    axis equal
    ylim([-1.2, 1.2] * sum(lval));
    xlim([-1, 2.5]);
    grid on
hold off

% iterate over the solution and plot with fixed time increments (the solver
% does vary its timesteps depending on the stiffness of the current system)
tcur = 0;
dt = 0.02;
start_fuse = 1;

while tcur < t_range(2)
    % update next time value
    tcur = tcur + dt;

    % wait a bit
    % empiracal tests shows computation per frame is around 12ms on my
    % old but gold i7 4900MQ if not in gif-mode
    if tcur == dt
        pause(1)
    else
        pause(max(0, dt/animspeed - 0.012))
    end
    
    % interpolate solution
    y_at_t = interp1(tSol, ySol, tcur, 'pchip');
    
    % get new positions 
    [cgpos, endpoints] ...
        = positions(y_at_t(1), y_at_t(2:(n+1)), lval, n);

    % update plot data
    plt_line.XData = endpoints(1, :); plt_line.YData = endpoints(2, :);
    plt_cart.XData = endpoints(1, 1); plt_cart.YData = endpoints(2, 1);
    plt_cons.XData = endpoints(1, 2:end); 
        plt_cons.YData = endpoints(2, 2:end);
    
    % when do we draw?
    drawnow

    if gifoutput
        % Capture the plot as an image 
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind, cm] = rgb2ind(im, 256);

        % Write to the GIF File 

        if start_fuse == 1
            start_fuse = 0;
            imwrite(imind, cm, gifname, 'gif', ...
                'Loopcount', inf, 'DelayTime', 1); 
%             for i = 1:round(1/dt)
%                 imwrite(imind, cm, gifname, 'gif', ...
%                     'WriteMode', 'append', 'DelayTime', dt/animspeed); 
%             end
        else 
            imwrite(imind, cm, gifname, 'gif', ...
                'WriteMode', 'append', 'DelayTime', dt/animspeed); 
        end

    end
    
end





