

num = xlsread('para_1.xlsx','Sheet1');
a = num(1,cn);
tau_0 = 296/200;



beta_v= num(4,cn);  %https://https://www.sciencedirect.com/science/article/pii/S0266352X20302925

beta_d=13.9*(1e-4)-beta_v;

K = num(5,cn); %Unit: N/m^2, P106 Table 3.1 
B = num(6,cn); % Default
G = num(7,cn); % Unit: N/m^2, P106 Table 3.1
v = num(8,cn);  % Default

c_d = num(13,cn); %Unit: J/m^3 K,  P633, Table 11.3

t = num(16:31,cn);
nn = length(t);

%% composite coefficient

alpha_d = K* beta_d; %% defined the equation 11.56
beta_e = a * beta_d+ beta_v;
%beta_e = beta_d;
M = B*K/(a-a^2*B);
m_d = c_d/tau_0;
eta = a*(1-2*v)/(2-2*v);
eta_d = alpha_d *(1-2*v)/(2*(1-v));
%% 
g_1 = eta/G;
g_2 = eta_d/G;
%% Four fluid flow coefficient
%visco_sity = 1e-4;
visco_sity = 1e-2;
kappa = k/visco_sity; % the devidator is the viscosity
kappa_pT = k_pT;
kappa_Tp = k_Tp/(tau_0);
kappa_T = k_T/(tau_0);
B_1 = kappa*M;
B_2 = -kappa_pT*M;
B_3 = g_1*a*M+1;
B_4 = g_2*a*M-beta_e*M;
B_5 = a*M;
B_6 = -kappa_Tp;
B_7 = kappa_T;
B_8 = K*beta_d*g_1*tau_0-beta_e*tau_0;
B_9 = g_2*K*beta_d*tau_0+m_d*tau_0;
B_10 = K*beta_d*tau_0;
BB_check = [B_1 B_2 B_3 B_4 B_5 B_6 B_7 B_8 B_9 B_10];
%% Build up the R matrix 
A_matrix = [B_1 B_2; B_6 B_7];
B_matrix = [B_3 B_4; B_8 B_9];
R_matrix = (B_matrix)\A_matrix;
fluid_diffu= R_matrix(1,1);
fluid_thermal_diffu= R_matrix(1,2);
thermal_diffu = R_matrix(2,2);
thermal_fluid_diffu = R_matrix(2,1);
diffu_ratio_dig = thermal_diffu/fluid_diffu;
diffu_ratio_off_dig = thermal_fluid_diffu/fluid_thermal_diffu;

zeta11_complex = diffu_ratio_dig*(beta_v+beta_d)*a/(beta_d*B)+0.04;
zeta11_simple = diffu_ratio_dig;

zeta12 = diffu_ratio_off_dig;

%zeta11 = 1/zeta11+0.03;
zeta11_simple
%zeta11_complex
%zeta12