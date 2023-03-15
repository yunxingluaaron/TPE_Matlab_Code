syms r R Cjj t m n t0 s
syms RR BV1 BV2 BV3

format long
n = 6;


%%
syms a alpha_d tau_0
syms beta_e beta_d beta_v eta eta_d
syms M K B G v
syms k_Tp k_pT k_T k
syms c_d
%a:alpha, Biot effective stress coefficient, 0<= alpha <=1
%alpha_d: drained thermal elastic effective stress coefficient, alpha_d = K * beta_d
%tau_0: reference absolute temperature

%beta_e: coeffcient of volumetric thermal expansion for variation in fluid content at constant frame volume 
%beta_d: coefficient of drained volumetric thermal expansion of porous medium frame
%beta_v: coefficient of volumetric thermal expansion for variation in fluid content

%eta
%eta_d

%M: pore pressure coefficient in the constitutive equation
%K: drained bulk modulus
%B: Skempton pore pressure coefficient, 0 <= B <=1
%G: shear modulus
%v: Poisson's ratio

%k_Tp: mechano-caloric coefficient
%k_pT:thermo-osmosis coefficient
%k_T: effective thermal conductivity, P625, Table 11.2
%k: permeability coefficient

%m_d: m_d = c_d/tau_0
%c_d: drained specific heat at constant strain,P633, Table 11.3

%% baisc coefficient
syms cn;
cn =1;
num = xlsread('para_1.xlsx','Sheet1');
a = num(1,cn);
tau_0 = num(2,cn); %Unit: K

beta_d= num(3,cn); %https://https://www.sciencedirect.com/science/article/pii/S0266352X20302925
beta_v= num(4,cn);  %https://https://www.sciencedirect.com/science/article/pii/S0266352X20302925

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
k_Tp = num(9,cn); %Unit: m^2/s.K %%https://www.sciencedirect.com/science/article/pii/S0266352X20302925
k_pT = num(10,cn);%https://https://www.sciencedirect.com/science/article/pii/S0266352X20302925
k_T = num(11,cn);%Unit: W/m K
k = num(12,cn); %Unit: m^2/Pa.s Coefficient of permeability P108 Table 3.2 

visco_sity = 1;
%visco_sity = 1e-4;
%gravity_density = 1e4;
kappa = k/visco_sity; % the devidator is the viscosity
kappa_pT = k_pT;
kappa_Tp = k_Tp/(tau_0);
kappa_T = k_T/(tau_0);
%kappa_T = k_T;
%kappa_Tp = k_Tp;


%B_1 = 1-(alpha_d*eta/(m_d*G)+beta_e/m_d)*(M*k_pT*G)/(a*M*eta+G)+kappa_T;

B_1 = kappa*M;
%B_1 = -kappa*M;

B_2 = -kappa_pT*M;
%B_2 = kappa_pT*M;

B_3 = g_1*a*M+1;

B_4 = g_2*a*M-beta_e*M;

B_5 = a*M;

B_6 = -kappa_Tp;
%B_6 = kappa_Tp;

B_7 = kappa_T;
%B_7 = -kappa_T;

B_8 = K*beta_d*g_1*tau_0-beta_e*tau_0;

B_9 = g_2*K*beta_d*tau_0+m_d*tau_0;

B_10 = K*beta_d*tau_0;

BB_check = [B_1 B_2 B_3 B_4 B_5 B_6 B_7 B_8 B_9 B_10];


%% Build up the R matrix 
A_matrix = [B_1 B_2; B_6 B_7];

%A_matrix = [B_1/(B_1*visco_sity) B_2/(B_1*visco_sity); B_6/(B_1*visco_sity) B_7/(B_1*visco_sity)];  % extract the hydraulic coefficienet from the matrix 

B_matrix = [B_3 B_4; B_8 B_9];

%R_matrix = pinv(B_matrix)*A_matrix;
R_matrix = (B_matrix)\A_matrix;

%R_matrix =[R_matrix(1,1)/R_matrix(1,1) R_matrix(1,2)/R_matrix(1,1);R_matrix(2,1)/R_matrix(1,1) R_matrix(2,2)/R_matrix(1,1)];


%R_matrix =[R_matrix_1(1,1)/B_1 R_matrix_1(1,2)/B_1;R_matrix_1(2,1)/B_1 R_matrix_1(2,2)/B_1];

%R_matrix =[R_matrix(1,1) R_matrix(1,2);R_matrix(2,1) R_matrix(2,2)];
%% EigenValue of R
EigenValue_R_new_matrix = eig(R_matrix);


CapitalLambda_1 =EigenValue_R_new_matrix(1);
CapitalLambda_2 =EigenValue_R_new_matrix(2);

if CapitalLambda_1<CapitalLambda_2
   
   temp_Y = CapitalLambda_1;
   CapitalLambda_1 = CapitalLambda_2;
   CapitalLambda_2 = temp_Y;
   
end 

CapitalLambda_1 = abs(CapitalLambda_1);

CapitalLambda_2 = abs(CapitalLambda_2);

if CapitalLambda_1<CapitalLambda_2
   
   temp_Y = CapitalLambda_1;
   CapitalLambda_1 = CapitalLambda_2;
   CapitalLambda_2 = temp_Y;
   
end 
%% Construction of the Transition P function

%Transition_P_12 = (CapitalLambda_1-R_11)/R_12;
%Transition_P_21 = (CapitalLambda_2-R_22)/R_21;  %% most original

Transition_P_12 = (CapitalLambda_2-R_matrix(2,2))/R_matrix(2,1);   %% Bunger Version
Transition_P_21 = (CapitalLambda_1-R_matrix(1,1))/R_matrix(1,2);  %% Bunger Version

%Transition_P_12 = (CapitalLambda_1-R_11)/R_12;  %% Jour well
%Transition_P_21 = (CapitalLambda_2-R_22)/R_21;

Transition_P = [1 Transition_P_12;Transition_P_21 1]; %% Bunger Version

%Transition_P = [1 Transition_P_21;Transition_P_12 1];


Transitio_P_inverse = pinv(Transition_P);