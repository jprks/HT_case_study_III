%{
    Title: Heat Transfer Case Study III - Heat Exchanger Analysis
    Author: James Emerson Parkus
    Date: 04/20/2019
    Purpose: The purpose of this script is to optimize the design of a
    cross-flow unmixed heat exchanger.
%}

%% Given
effect = 0.75; % [-] - Effectiveness of Heat Exchanger
ACH = 0.35; % [-] - Percent of Volume replaced per hour
A_house = 2000; % [ft^2] - Area of the House
A_house = A_house/0.0929; % [m^2] - m^2 conversion
H_house = 10*0.3048; % [m] - Height of the house
c_p = 1007; % [J/(kg K)] - Specific heat of air
rho = 1.1614; % [kg/m^3] - Density of air
T_ci = 273.15; % [K] - Temperature of Supply air inlet
T_hi = 293.15; % [K] - Temperature of House air inlet
nu = 15.89*10^-6; % [m^2/s] - Kinematic Viscosity of Air
Pr = 0.707; % [-] - Prantl number for air

global h_core w_core l_core
h_core = 9*0.0254; % [m] - Height of HX core
w_core = 6*0.0254; % [m] - Width of HX core
l_core = 6*0.0254; % [m] - Length of HX core
%% Analysis
V_house = H_house*A_house; % [m^3] - Volume of the house

V_flow = V_house*ACH; % [m^3/s] - Volumetric flowrate
m_flow = V_flow*rho; % [kg/s] - Mass flowrate
C = m_flow*c_p; % [J/(s K)] - Specific heat capacity

q_max = C*(T_hi - T_ci); % [W] - Maximum heat transfer
q_actual = effect*q_max; % [W] - Actual heat transfer

% find NTU by numerical iteration of Eqn. 11.32 pg. 664


for n = 2:2:10
    sep_d = n*h_core; % [m] - Separation distance of the HX plates
    A_hx = area_hx(sep_d); % [m] - Area of the HX
    
    U_req = NTU*C/A_hx; % [W/(m K)] - Effective convective coefficient
    vel = m_flow/(rho*w_core*sep_d); % [m/s] - Flow velocity
    
    D_h = 2*w_core*sep_d/(w_core + sep_d); % [m] - Hydraulic Diameter
    
    Re_limit = 2300; % [-] - Reynold's Number limit for Laminar flow
    if Re < Re_limit
        x_fd_h = 0.05*Re*D_h;
        x_fd_t = 0.05*Re*D_h*Pr;
    else
        x_fd_h = 10*D_h;
        x_fd_t = x_fd_h;
    end
    
    if x_fd_h > l_core
        % Nu correlation for under-developed flow
    else
        % Nu correlation for fully-developed flow
    end
    
    h = Nu*k/D_h; % [W/(m^2 K)] - Convective coefficient
    
    U_calc = (2/h)^-1; % [W/(m^2 K)] - Effective convective coefficient
    
    f = 64/Re;
    del_P = f*rho*vel^2*w_core/(2*D_h);
    Power = del_P*V_flow;
end

%% Function
function A = area_hx(d)
global h_core w_core l_core
A = (l_core*h_core*w_core)/d;
end

function x_ntu = NTU_correlation(x)
x_ntu = 1 - exp(x^0.22-(exp(-x^0.78-1)));
end







