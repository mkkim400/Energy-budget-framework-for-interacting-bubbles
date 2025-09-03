clear all; clc; close all; format compact;

%=====================================================================
% Run KMsolver: Bubble-Bubble Interaction Simulation
%=====================================================================
% Outputs (from solver):
%   t   : time
%   y1  : bubble radius (bubble 1)
%   y2  : bubble wall velocity
%   y3  : bubble pressure
%
% Possible PR values: 0.05, 0.06, 0.075, 0.1, 0.25, 0.5, 1, 2.5, 5,
%                    10, 25, 50, 100, 250, 500
%=====================================================================

nBub    = 2;      % number of bubbles (2 or 3)
in_plot = 100;    % plotting frequency (unused here)
R01     = 0.5/1e3; % initial radius of bubble 1 (m)

%=====================================================================
% Parameter sweeps
%=====================================================================
DDi = 6:1:40;              % normalized bubble separation factors
RRi = 1.0:0.5:5.0;         % bubble radius ratios (R02/R01)
AA  = [2.4:0.2:2.8 3.4:0.2:4]; % exponents for PR sweep
PRi = 10.^AA;              % pressure ratios
in_leng1 = numel(DDi);
in_leng2 = numel(RRi);
in_leng3 = numel(PRi);

%=====================================================================
% Initial & physical conditions
%=====================================================================
equilb   = 1;       % 0 = equilibrium (Rdot=0), 1 = non-equilibrium
grow     = 1;       % growth (0) or collapse (1)
p_init1  = 1e5;     % initial pressures (Pa)
p_init2  = 1e5;
p_init3  = 1e5;
p_v      = 0;       % vapor pressure (Pa)

rho_l    = 998;     % liquid density (kg/m^3)
sig      = 1e-15;   % surface tension (N/m) - almost zero
mu_l     = 1e-15;   % viscosity (Pa·s) - almost zero
c_l      = 1510.0;  % liquid sound speed (m/s)
gam_v    = 1.4;     % ratio of specific heats

%%% Solver outputs (dimensional):
% y1 = R1(t)        : bubble 1 radius
% y2 = f1_st        : bubble 1 velocity term
% y3 = p1_b_st      : bubble 1 pressure
% y4 = R2(t)        : bubble 2 radius
% y5 = f2_st        : bubble 2 velocity term
% y6 = p2_b_st      : bubble 2 pressure
% y7...y14 = higher derivatives / extra outputs
%=====================================================================

for k = 1:in_leng3   % loop over PR values
    k
    for j = 1:in_leng1   % loop over separation distances
        for i = 1:in_leng2   % loop over radius ratios
            %---------------------------------------------------------
            % Set parameters for this run
            %---------------------------------------------------------
            DD = DDi(j);       % bubble separation factor
            RR = RRi(i);       % bubble radius ratio
            PR = PRi(k);       % pressure ratio

            D     = DD*R01;     % actual bubble separation (m)
            R02   = R01*RR;     % bubble 2 radius (m)
            p_inf = PR*1e5;     % ambient pressure (Pa)

            %---------------------------------------------------------
            % Call solver depending on number of bubbles
            %---------------------------------------------------------
            if nBub == 2
                [t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,rho_l] = ...
                    KM_Ida_solver(p_inf,equilb,sig,rho_l,mu_l,c_l,gam_v,...
                                  p_init1,p_init2,p_v,R01,R02,D,grow);
            elseif nBub == 3
                [t,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,rho_l] = ...
                    KM_Ida_solver_3(PR,equilb,sig,mu_l,c_l,gam_v,...
                                    p_init1,p_init2,p_init3,p_v,...
                                    R01,R02,R03,D,grow);
            end

            %---------------------------------------------------------
            % Quick diagnostic outputs
            %---------------------------------------------------------
            [R_min, in] = min(y1/R01)  %#ok<NOPRT>
            t(in)/(0.915*R01*sqrt(rho_l/(PR*1e5-p_init1)))

            %---------------------------------------------------------
            % Time & radius data for plotting
            %---------------------------------------------------------
            t_plot = t*1e6;      % convert to microseconds
            R_plot = y1*1e3;     % convert to mm

            % Example plots (commented)
            % figure,plot(t_plot,y1/(1e-3),'b-',t_plot,y4/(1e-3),'r-','LineWidth',1.5);
            % xlabel('$t \ ({\rm \mu s})$','interpreter','latex');
            % ylabel('$R \ ({\rm mm})$','interpreter','latex');
            % set(gca,'FontName','Times New Roman','FontSize',20);

            % Save reduced data (every 3rd point)
            t_leng = numel(t_plot);
            out    = [t_plot(1:3:t_leng) R_plot(1:3:t_leng)];

            %---------------------------------------------------------
            % Organize output folder structure
            %---------------------------------------------------------
            if grow==0 && equilb==0
                name1 = ['D:\2_cav\1_b-b_interact\0_KMIda_data\1 grow_gas_equilb\no_pv_mu_sig\DD=',num2str(DD)];
            elseif grow==0 && equilb==1
                name1 = ['D:\2_cav\1_b-b_interact\0_KMIda_data\1 grow_gas_non_equilb\no_pv_mu_sig\DD=',num2str(DD)];
            elseif grow==1 && equilb==0
                name1 = ['D:\2_cav\1_b-b_interact\0_KMIda_data\2 col_gas_equilb\no_pv_mu_sig\DD=',num2str(DD)];
            elseif grow==1 && equilb==1
                name1 = ['D:\2_cav\1_b-b_interact\0_KMIda_data\2 col_gas_non_equilb\no_pv_mu_sig\DD=',num2str(DD)];
            end

            if nBub == 2
                name2 = ['\R01=',num2str(R01),...
                         '_RR=',num2str(RR),...
                         '_PR=',num2str(PR),...
                         '\pb1=',num2str(p_init1/1e3),...
                         '_rho=',num2str(rho_l)];
                mkdir([name1,name2]);
            elseif nBub == 3
                name2 = ['\R01=',num2str(R01),...
                         '_R02=',num2str(R02),...
                         '_PR=',num2str(PR),...
                         '_n=',num2str(nBub),...
                         '\pb1=',num2str(p_init1/1e3),...
                         '_rho=',num2str(rho_l)];
                mkdir([name1,name2]);
            end

            %---------------------------------------------------------
            % Write dimensional results to text files
            %---------------------------------------------------------
            fileID0=fopen([name1,name2,'\R_Q.txt'],'w');
            fprintf(fileID0,'t(s)\t  R(m)\t  R_dot(m/s)\t  R_2dot(m/s^2)\t  p_b(pa)\t  p_b_dot(pa/s)\r\n');

            fileID1=fopen([name1,name2,'\R_Q_non.txt'],'w');
            fprintf(fileID1,'t/t_non\t  t/t_c\t  R/R0\t  R_dot_non\t  R_2dot_non\t  p_b_non\t  p_b_dot_non\r\n');

            if nBub == 2
                output0 = [t y1 y2 y3 y4 y5 y6 y7 y8 y9 y10];
                fprintf(fileID0,'%.8e\t %.8e\t %.8e\t %.8e\t %.8e\t %.8e %.8e %.8e %.8e %.8e %.8e\r\n',output0');
            elseif nBub == 3
                output0 = [t y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14];
                fprintf(fileID0,'%.8e\t %.8e\t %.8e\t %.8e\t %.8e\t %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\r\n',output0');
            end

            %---------------------------------------------------------
            % Non-dimensional scaling
            %---------------------------------------------------------
            Lc    = R01;  % length scale
            t_non = Lc*sqrt(rho_l/(p_inf-p_init1)); % natural time scale
            tc    = 0.915*t_non;                   % collapse time
            Vc    = Lc/tc;                         % velocity scale

            R_non       = y1/Lc;
            R_dot_non   = y2/Vc;
            R_2dot_non  = y4/(Lc/tc^2);
            p_b_non     = y3/(rho_l*Vc^2);
            p_b_dot_non = y5/(rho_l*Vc^2/tc);

            output1 = [t/t_non t/tc R_non R_dot_non R_2dot_non p_b_non p_b_dot_non];
            fprintf(fileID1,'%.8e\t %.8e\t %.8e\t %.8e\t %.8e\t %.8e\t %.8e\r\n',output1');

            fclose(fileID0);
            fclose(fileID1);

            clear output0 output1
        end
    end
end
