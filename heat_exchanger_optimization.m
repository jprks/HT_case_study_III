%{
    Title: Heat Transfer Case Study III - Heat Exchanger Analysis
    Author: James Emerson Parkus
    Date: 04/20/2019
    Purpose: The purpose of this script is to optimize the design of a
    cross-flow unmixed heat exchanger.
%}
clc
clear
close all

%% Given
effect = 0.75; % [-] - Effectiveness of Heat Exchanger
ACH = 0.35; % [-] - Percent of Volume replaced per hour
A_house = 2000; % [ft^2] - Area of the House
A_house = A_house*0.0929; % [m^2] - m^2 conversion
H_house = 10*0.3048; % [m] - Height of the house
c_p = 1007; % [J/(kg K)] - Specific heat of air @ 300K
rho = 1.1614; % [kg/m^3] - Density of air @ 300K
T_ci = 273.15; % [K] - Temperature of Supply air inlet
T_hi = 293.15; % [K] - Temperature of House air inlet
nu = 15.89*10^-6; % [m^2/s] - Kinematic Viscosity of Air
k = 0.026; % [W/(m K)] - Thermal Conductivity of Aluminum 2024 - T6 @ 300K

global h_core w_core l_core Pr
h_core = 9*0.0254; % [m] - Height of HX core
w_core = 6*0.0254; % [m] - Width of HX core
l_core = 6*0.0254; % [m] - Length of HX core
Pr = 0.707; % [-] - Prantl number for air
%% Analysis
V_house = H_house*A_house; % [m^3] - Volume of the house

V_flow = V_house*ACH/3600; % [m^3/s] - Volumetric flowrate
m_flow = V_flow*rho; % [kg/s] - Mass flowrate
C = m_flow*c_p; % [J/(s K)] - Specific heat capacity

q_max = C*(T_hi - T_ci); % [W] - Maximum heat transfer
q_actual = effect*q_max; % [W] - Actual heat transfer

%% find NTU by numerical iteration of Eqn. 11.32 pg. 664
error = 1;
inc = 0.001;
i = 1;
x_guess_1 = 1;
x_guess_2 = 1.5;

x_sol_1 = NTU_correlation(x_guess_1);
x_sol_2 = NTU_correlation(x_guess_2);

del_E_1 = abs(effect - x_sol_1);
del_E_2 = abs(effect - x_sol_2);
x_old = [x_guess_1 x_guess_2];

while error > 0.001 && i < 500000   
    if del_E_1 < del_E_2                % if the delta got larger as x went from 1 to 1.5
        x_new = x_old(i,1) - inc;
        TF = 0;
    else                                % if the delta got larger as x went from 1.5 to 1
        x_new = x_old(i,2) + inc;
        TF = 1;
    end
    
    x_sol_new(i,1) = NTU_correlation(x_new);
    error = abs(effect - x_sol_new(i,1));
    i = i + 1;
    if TF == 0
        x_old(i,1) = x_new;
    else
        x_old(i,2) = x_new;
    end
end

NTU = x_new - inc;

%% Find U_calc and U_required
i = 1;
sep_d(i,1) = 0.500; % [m] - Separation distance of the plates
TF = 0;
while TF == 0 && sep_d(i,1) > 0.0005 && i < 500000
    A_hx = area_hx(sep_d(i,1)); % [m] - Area of the HX
    U_req = NTU*C/A_hx; % [W/(m K)] - Effective convective coefficient
    vel = m_flow/(rho*w_core*sep_d(i,1)); % [m/s] - Flow velocity
    D_h = 2*w_core*sep_d(i,1)/(w_core + sep_d(i,1)); % [m] - Hydraulic Diameter
    
    Re = vel*D_h/nu; % [-] - Reynold's Number of air
    Re_limit = 2300; % [-] - Reynold's Number limit for Laminar flow
    if Re < Re_limit
        x_fd_h = 0.05*Re*D_h;
        x_fd_t = 0.05*Re*D_h*Pr;
        Nu = 3.66;
    else
        x_fd_h = 10*D_h;
        x_fd_t = x_fd_h;
        Nu = Nu_turbulent(Re);
    end
    
    FD(i,1:3) = [x_fd_h x_fd_t Re];
    
    h = Nu*k/D_h; % [W/(m^2 K)] - Convective coefficient
    
    U_calc(i,1) = (2/h)^-1; % [W/(m^2 K)] - Effective convective coefficient
    
    f = 64/Re;
    del_P = f*rho*vel^2*w_core/(2*D_h);
    Power = del_P*V_flow;
    u_error = abs(U_calc(i,1) - U_req);
    Req_U(i,1) = U_req;
    
    if i > 1
        if Req_U(i-1,1) > U_calc(i-1,1) && Req_U(i,1) < U_calc(i,1) || Req_U(i-1,1) < U_calc(i-1,1) && Req_U(i,1) > U_calc(i,1)
            TF = 1;
        else
            sep_d(i+1,1) = sep_d(i,1) - 0.001;
        end
    else
        sep_d(i+1,1) = sep_d(i,1) - 0.001;
    end
    
    sep(i,1) = sep_d(i,1);
    
    Reynolds(i,1) = Re;
    h_conv(i,1) = h;
    velocity(i,1) = vel;
    
    i = i + 1;
end
i = i - 1;
%% Important Data
U_c = [U_calc(i-1,1);U_calc(i,1)];
U_r = [Req_U(i-1,1);Req_U(i,1)];
Ren = [Reynolds(i-1,1);Reynolds(i,1)];
H = [h_conv(i-1,1);h_conv(i,1)];
Velo = [velocity(i-1,1);velocity(i,1)];
Separation = [sep(i-1,1);sep(i,1)];

Important = table(Separation,U_c,U_r,Ren,H,Velo);

%% Plotting
hold on
grid on
mm_sep = sep*1000;
plot(mm_sep,Req_U);
plot(mm_sep,U_calc);
title('Convergence of U_{required} and U_{calculated}');
xlabel('Plate Separation Distance [mm]');
ylabel('U [W/(m^2 K)]');

%% Function
function A = area_hx(d)
global h_core w_core l_core
A = (l_core*h_core*w_core)/d;
end

function x_ntu = NTU_correlation(x)
x_ntu = 1 - exp(x^0.22*(exp(-x^0.78)-1));
end

function Nu_turb = Nu_turbulent(Re)
global Pr
Nu_turb = 0.023*Re^(4/5)*Pr^0.4;
end





