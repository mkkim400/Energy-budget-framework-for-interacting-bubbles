function [t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,rho_l] = KM_Ida_solver( ...
    p_inf,equilb,sig,rho_l,mu_l,c_l,gam_v,p_init1,p_init2,p_v,R01,R02,D,grow)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Purpose: Solve two-bubble interaction using the Keller–Miksis equations
%  - Bubble 1: radius R01
%  - Bubble 2: radius R02
%  - Bubbles separated by distance D
%  - Includes compressibility, inertia, and bubble-bubble interaction
%
% Output variables:
%   t   : time vector (s)
%   y1  : bubble 1 radius (m)
%   y2  : auxiliary variable f1 = R1^2 * R1dot
%   y3  : bubble 1 pressure (Pa)
%   y4  : bubble 2 radius (m)
%   y5  : auxiliary variable f2 = R2^2 * R2dot
%   y6  : bubble 2 pressure (Pa)
%   y7  : bubble 1 wall acceleration (m/s^2)
%   y8  : bubble 1 pressure time derivative (Pa/s)
%   y9  : bubble 2 wall acceleration (m/s^2)
%   y10 : bubble 2 pressure time derivative (Pa/s)
%   rho_l : liquid density (kg/m^3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial bubble wall velocities
if equilb == 0
    % Equilibrium initial condition (both bubbles at rest)
    R1_dot = 0;
    R2_dot = 0;
elseif equilb == 1
    % Non-equilibrium initial condition (velocity from pressure mismatch)
    R1_dot = (p_init1 - p_inf) / (rho_l * c_l);
    R2_dot = (p_init2 - p_inf) / (rho_l * c_l);
end
t_init = 0;

%% Initial bubble internal pressures
pg1 = p_init1; 
pg2 = p_init2;
pv  = p_v;       % vapor pressure
pb1 = pg1 + pv;  % total pressure in bubble 1
pb2 = pg2 + pv;  % total pressure in bubble 2

%% Characteristic scales (for nondimensionalization)
Lc = R01;                       % length scale = initial radius of bubble 1
tc = Lc * sqrt(rho_l/p_inf);    % characteristic time scale
Vc = Lc / tc;                   % velocity scale
D  = D / Lc;                    % nondimensional bubble separation
p_inf_st = p_inf / (rho_l*Vc^2);% nondimensional ambient pressure
t_final = 3.5*tc;               % final simulation time

%% Dimensionless groups
Re_l  = (rho_l*Vc*Lc)/mu_l;     % Reynolds number
We    = (rho_l*Vc^2*Lc)/sig;    % Weber number
c_st  = c_l / Vc;               % nondimensional sound speed

%% Initialize storage variables for diagnostics
t_tmp   = 0.0; 
tout    = []; 
pinfout = []; 
Twout   = []; 
k1_wout = [];

%% Initial conditions for ODE solver (nondimensionalized)
R1_2dot_save  = 0; 
p1_b_dot_save = 0; 
R2_2dot_save  = 0; 
p2_b_dot_save = 0;

y0(1,1) = R01/R01;             % bubble 1 radius (normalized)
y0(1,2) = R01^2*R1_dot/Lc^2/Vc;% f1 = R^2*Rdot
y0(1,3) = pb1/(rho_l*Vc^2);    % bubble 1 pressure
y0(1,4) = R02/R01;             % bubble 2 radius (normalized)
y0(1,5) = R02^2*R2_dot/Lc^2/Vc;% f2 = R^2*Rdot
y0(1,6) = pb2/(rho_l*Vc^2);    % bubble 2 pressure

%% Integrate the ODE system
tspan   = [t_init t_final/tc]; % nondimensional time interval
options = odeset('RelTol',1e-15,'AbsTol',1e-15);
[t,y]   = ode15s(@derivs,tspan,y0,options);

%% Convert back to dimensional variables
t  = t * tc;                        % time (s)
y1 = y(:,1) * Lc;                   % bubble 1 radius (m)
y2 = y(:,2) * Lc^2 * Vc;            % f1 dimensional
y3 = y(:,3) * rho_l * Vc^2;         % bubble 1 pressure (Pa)
y4 = y(:,4) * Lc;                   % bubble 2 radius (m)
y5 = y(:,5) * Lc^2 * Vc;            % f2 dimensional
y6 = y(:,6) * rho_l * Vc^2;         % bubble 2 pressure (Pa)
y7 = R1_2dot_save * Vc/tc;          % bubble 1 wall acceleration
y8 = p1_b_dot_save * rho_l*Vc^2/tc; % bubble 1 pressure derivative
y9 = R2_2dot_save * Vc/tc;          % bubble 2 wall acceleration
y10= p2_b_dot_save * rho_l*Vc^2/tc; % bubble 2 pressure derivative

%% Nested function: ODE right-hand side
    function f = derivs(t,y)
        % State variables (dimensionless)
        R1_st   = y(1);  % bubble 1 radius
        f1_st   = y(2);  % f1 = R1^2*R1dot
        p1_b_st = y(3);  % bubble 1 pressure
        R2_st   = y(4);  % bubble 2 radius
        f2_st   = y(5);  % f2 = R2^2*R2dot
        p2_b_st = y(6);  % bubble 2 pressure
        
        % Bubble wall velocities (nondim)
        R1_dot_st = f1_st / R1_st^2;
        R2_dot_st = f2_st / R2_st^2;
        
        %% Pressure derivatives inside bubbles (from adiabatic law)
        f(3,1)   = -3.0/R1_st * (gam_v*p1_b_st*R1_dot_st);
        f(6,1)   = -3.0/R2_st * (gam_v*p2_b_st*R2_dot_st);
        p1_b_dot = f(3,1);    
        p2_b_dot = f(6,1);
        
        %% Radius derivatives
        f(1,1) = R1_dot_st;
        f(4,1) = R2_dot_st;
        
        %% f_dot equations (coupled Keller–Miksis form)
        % Build linear system A*Fdot = B
        A(1,:) = [1/R1_st - R1_dot_st/(c_st*R1_st), 1/D];
        A(2,:) = [1/D, 1/R2_st - R2_dot_st/(c_st*R2_st)];
        
        B(1,:) = (1+R1_dot_st/c_st)*(p1_b_st - p_inf_st) + ...
                 R1_st*p1_b_dot/c_st + ...
                 0.5*R1_dot_st^2*(1 - 3*R1_dot_st/c_st);
        B(2,:) = (1+R2_dot_st/c_st)*(p2_b_st - p_inf_st) + ...
                 R2_st*p2_b_dot/c_st + ...
                 0.5*R2_dot_st^2*(1 - 3*R2_dot_st/c_st);
        
        % Solve for f_dot
        F_dot = inv(A) * B.';  % 2x1 vector
        
        f(2,1) = F_dot(1);
        f(5,1) = F_dot(2);

        %% Accelerations (second derivatives of radius)
        R1_2dot = (f(2,1) - 2*R1_st*R1_dot_st^2) / R1_st^2;
        R2_2dot = (f(5,1) - 2*R2_st*R2_dot_st^2) / R2_st^2;
        
        %% Store diagnostic values for later use
        Nt_tmp = numel(tout);
        if (Nt_tmp > 1 && tout(Nt_tmp) == tout(Nt_tmp-1))
            Nt_tmp = Nt_tmp - 1;
        end
        
        if (t ~= t_tmp)
            if (t <= t_tmp)
                tout(Nt_tmp,1)    = t;
                pinfout(Nt_tmp,1) = p_inf;
                R1_2dot_save(Nt_tmp,1)  = R1_2dot;
                p1_b_dot_save(Nt_tmp,1) = p1_b_dot;
                R2_2dot_save(Nt_tmp,1)  = R2_2dot;
                p2_b_dot_save(Nt_tmp,1) = p2_b_dot;
            else
                tout(Nt_tmp+1,1)    = t;
                pinfout(Nt_tmp+1,1) = p_inf;
                R1_2dot_save(Nt_tmp+1,1)  = R1_2dot;
                p1_b_dot_save(Nt_tmp+1,1) = p1_b_dot;
                R2_2dot_save(Nt_tmp+1,1)  = R2_2dot;
                p2_b_dot_save(Nt_tmp+1,1) = p2_b_dot;
            end
        end
        t_tmp = t;    
    end
end
