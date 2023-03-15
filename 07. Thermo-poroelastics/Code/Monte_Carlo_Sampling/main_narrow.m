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

num = xlsread('para_overall_loop_reduce_six_highDifference.xlsx','Sheet1');
t = num(16:31,1);
nn = length(t);


Over_all_data={};

NNA_matrix = zeros(1,100000);
%Sample_size = 100;

%Sample_iterations=0;

for iii = 1:300
    

NNA_matrix(iii) = iii;

Sample_rate_1 = unifrnd(0,1);
Sample_rate_2 = unifrnd(0,1);
Sample_rate_3 = unifrnd(0,1);
Sample_rate_4 = unifrnd(0,1);
Sample_rate_5 = unifrnd(0,1);
Sample_rate_6 = unifrnd(0,1);
Sample_rate_7 = unifrnd(0,1);
Sample_rate_8 = unifrnd(0,1);
Sample_rate_9 = unifrnd(0,1);
Sample_rate_10 = unifrnd(0,1);
Sample_rate_11 = unifrnd(0,1);


a = num(1,5)+(num(1,7)-num(1,5))*Sample_rate_1;
                               
beta_v= num(4,5)+(num(4,7)-num(4,5))*Sample_rate_2;

K = num(5,6)+(num(5,7)-num(5,6))*Sample_rate_3;
                                
B = num(6,5)+(num(6,7)-num(6,5))*Sample_rate_4;
                                
G = num(7,6)+(num(7,7)-num(7,6))*Sample_rate_5;
                                
v = num(8,6)+(num(8,7)-num(8,6))*Sample_rate_6;
                              
k_Tp = num(9,5)+(num(9,7)-num(9,5))*Sample_rate_7;
                                
k_pT = num(10,5)+(num(10,7)-num(10,5))*Sample_rate_8;
                                
k_T = num(11,5)+(num(11,7)-num(11,5))*Sample_rate_9;
                                
k = num(12,5)+(num(12,7)-num(12,5))*Sample_rate_10;
                                
c_d = num(13,7)+(num(13,5)-num(13,7))*Sample_rate_11;
                                
tau_0=num(2,5);

beta_d=num(2,5);

%beta_d=13.9*(1e-4)-beta_v;
%beta_d=13.9*(1e-4)-beta_v;
                                
                                
Zeta_One = (k)^2*(k_T)^2/((k_Tp)^2*c_d*beta_v);

Zeta_One_one = (k)*(k_T)/((k_Tp)*c_d*beta_v);
                               
Zeta_Two= a*B/K;
                                
                               % input_matrix(1,cn_0-4)=a;
                               % input_matrix(2,cn_1-4)=beta_v ;
                               % input_matrix(3,cn_2-4)=B;
                               % input_matrix(4,cn_3-4)=K;
                               % input_matrix(5,cn_4-4)=G;
                               % input_matrix(6,cn_5-4)=v;
                               % input_matrix(7,cn_6-4)=k_Tp;
                               % input_matrix(8,cn_7-4)=k_pT;
                               % input_matrix(9,cn_8-4)=k_T;
                               % input_matrix(10,cn_9-4)=k;
                               % input_matrix(11,cn_10-4)=c_d;
                                
                                %input_matrix_this_time ={input_matrix(1,cn_0-4),input_matrix(2,cn_1-4),input_matrix(3,cn_2-4),input_matrix(4,cn_3-4),input_matrix(5,cn_4-4),input_matrix(6,cn_5-4),input_matrix(7,cn_6-4),input_matrix(8,cn_7-4)...
                                   %input_matrix(9,cn_8-4), input_matrix(10,cn_9-4),input_matrix(11,cn_10-4)};
                                   

                                
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
visco_sity =1;
%kappa = k/visco_sity; % the devidator is the viscosity
%kappa_pT = k_pT;
%kappa_T = k_T;
%kappa_Tp = k_Tp;
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

input_matrix_this_time ={a,beta_d,beta_v,K,B,G,v,k_Tp,k_pT,k_T,k,c_d,fluid_diffu,thermal_diffu,fluid_thermal_diffu,thermal_fluid_diffu,diffu_ratio_dig,diffu_ratio_off_dig};

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
%if isreal(CapitalLambda_1) & isreal(CapitalLambda_2) & CapitalLambda_1>0 & CapitalLambda_2 >0 
    
 
    
%    if CapitalLambda_1<CapitalLambda_2
   
%         temp_Y = CapitalLambda_1;
%         CapitalLambda_1 = CapitalLambda_2;
%         CapitalLambda_2 = temp_Y;
%    end 
    
%else 
%    input_matrix_this_time=[];
    
%end 


%% Construction of the Transition P function

Transition_P_12 = (CapitalLambda_2-R_matrix(2,2))/R_matrix(2,1);   %% Bunger Version
Transition_P_21 = (CapitalLambda_1-R_matrix(1,1))/R_matrix(1,2);  %% Bunger Version

Transition_P = [1 Transition_P_12;Transition_P_21 1]; 

Transitio_P_inverse = pinv(Transition_P);

%%Starting Calculations
format long
% Temperature in time domain    
    RR_T = zeros(1,100);
    % Pore pressure in time domain  
    RR_p = zeros(1,100);
    % Radial stress in time domain
    sigma_rr = zeros(1,100);
    
    u_rr = zeros(1,100);
% Store the temperature data point in one matrix

    RRR_T = zeros(nn,100);
    RRR_p_1 = zeros(nn,100);
    RRR_p_2 = zeros(nn,100);
    RRR_p_3 = zeros(nn,100);
    RRR_rr = zeros(nn,100);
    
    c1_check_matrix = zeros(6*nn,100);
    c2_check_matrix = zeros(6*nn,100);
    fs_check_matrix = zeros(6*nn,100);
    
    a7_check_matrix = zeros(6*nn,100);
    a8_check_matrix = zeros(6*nn,100);
    a9_check_matrix = zeros(6*nn,100);
    
    
    AA1 = 2*G*v*eta/(G*(1-2*v))-(2*G-2*G*v)*eta/(G*(1-2*v));
    
    AA2 = 2*G*v*eta_d/(G*(1-2*v))-(2*G-2*G*v)*eta/(G*(1-2*v));
    
    AA3 = (-2*G*v+2*G)*eta/(G*(1-2*v))-a;
    
    AA4 = (-2*G*v+2*G)*eta_d/(G*(1-2*v))-alpha_d;
    
    AA5 = (2*G)/(1-2*v);
    
   %% Loading modes input values:
   Input = xlsread('para_overall_loop_reduce_six_highDifference.xlsx','Sheet1');

   Loading_mode_1_input = [Input(1,3), 0, 0];
   Loading_mode_2_input = [0, Input(2,3), 0];
   Loading_mode_3_input = [0, 0, Input(3,3)];
   

   
  
%% Mode_1 Loading: pore pressure loading
for mmm = 1:1:nn
    
    t0 = t(mmm);
  
    
    for r = 1:100
    
    % Temperature in Laplace domain
    tilde_TT_total = 0;
    % Pore pressure in Laplace domain
    tilde_pp_total  = 0;
    % Radial stress in Laplace domain
  
    tilde_sigma_rr_part_1_total = 0;
    
    tilde_sigma_rr_part_2_total = 0;
    
    tilde_sigma_rr_part_3_total = 0;
    
    tilde_sigma_rr_part_4_total = 0;
    
    tilde_sigma_rr_part_5_total = 0;
    
 
      for j = 1:n
    
        kk = floor((j+1)/2);
        mi = min(j,n/2);
    
        Cjj = 0;
    
            for i = kk:mi
      
                    Cj = ((i^(n/2))*factorial(2*i))/(factorial(n/2-i)*factorial(i)*factorial(i-1)*factorial(j-i)*factorial(2*i-j));
                    Cjj = Cjj + Cj;
            end 
    
                    Cjj = Cjj *((-1)^(j+n/2));
    
                     s = j*log(2)/t0;
    
    
    %% obtain c1 & c2 & fs based on different boundary conditions
     
      
      H_matrix = Transition_P\((B_matrix)\[B_5;B_10]);
     
     
     H_1 = H_matrix(1);
     
 
     H_2 = H_matrix(2);
     
 
     r_at_boundary = 100/250;
     
     B_0_1_boundary = besseli(0,sqrt(s/CapitalLambda_1)*r_at_boundary); % B_0_1: stand for Besseli function, 0 order, contains Lambda 1

     B_0_2_boundary = besseli(0,sqrt(s/CapitalLambda_2)*r_at_boundary);   % B_0_1: stand for Besseli function, 0 order, contains Lambda 2
  
     B_1_1_boundary = besseli(1,sqrt(s/CapitalLambda_1)*r_at_boundary)/sqrt(s/CapitalLambda_1); % B_1_1: stand for Besseli function, 1 order, contains Lambda 1
                         
     B_1_2_boundary = besseli(1,sqrt(s/CapitalLambda_2)*r_at_boundary)/sqrt(s/CapitalLambda_2); % B_1   _1: stand for Besseli function, 1 order, contains Lambda 2
 
 
 
 
    a_1 = Transition_P(1,1)*B_0_1_boundary;
 
    a_2 = Transition_P(1,2)*B_0_2_boundary;
 
    a_3 = Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2;
   
 
 
    a_4 = Transition_P(2,1)*B_0_1_boundary;
 
    a_5 = Transition_P(2,2)*B_0_2_boundary;
 
    a_6 = Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2;
    
    a_7 = (1/r_at_boundary)*B_1_1_boundary*(AA1*Transition_P(1,1)+AA2*Transition_P(2,1))+B_0_1_boundary*(AA3*Transition_P(1,1)+AA4*Transition_P(2,1));
    
    a_8 = (1/r_at_boundary)*B_1_2_boundary*(AA1*Transition_P(1,2)+AA2*Transition_P(2,2))+B_0_2_boundary*(AA3*Transition_P(1,2)+AA4*Transition_P(2,2));
   
    a_9 = -AA5+(0.5*AA1*((Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2))+ 0.5*AA2*((Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2))+ AA3*((Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2))+...
    AA4*((Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2)));
    
    for c_i = 1:6*nn
        
        a7_check_matrix(c_i,r) = a_7;
        a8_check_matrix(c_i,r) = a_8;
        a9_check_matrix(c_i,r) = a_9;
        
    end
 
   
    Coefficient_Matrix = [a_1 a_2 -a_3; a_4 a_5 -a_6; a_7 a_8 -a_9];
    
     
 
    
    % thermal, pore and radial here are three different boundary conditions: 
    % x: thermal conditions;
    % y: pore pressure conditions;
    % z: radial stress boundary conditions;
    
    x = Loading_mode_1_input(1)/s;
    y = Loading_mode_1_input(2)/s;
    z = Loading_mode_1_input(3)/s;
 
    Boundary_Matrix = [x;y;z];
    
    
 
    Coefficient_result_matrix = Coefficient_Matrix\Boundary_Matrix;
    
    c_1 = Coefficient_result_matrix(1);
 
    c_2 = Coefficient_result_matrix(2);
 
    fs = Coefficient_result_matrix(3);
    
    
    
    for c_i = 1:6*nn
        
        c1_check_matrix(c_i,r) = c_1;
        c2_check_matrix(c_i,r) = c_2;
        fs_check_matrix(c_i,r) = fs;
        
    end

    
   
   
   

                     
 %% In Laplace domain, obtian the p_tilde, T_tilde and sigma_rr_tilde
                    r_in_field = r/250;
                   
                    
                    B_0_1_field = besseli(0,sqrt(s/CapitalLambda_1)*r_in_field); % B_0_1: stand for Besseli function, 0 order, contains Lambda 1

                    B_0_2_field = besseli(0,sqrt(s/CapitalLambda_2)*r_in_field);   % B_0_1: stand for Besseli function, 0 order, contains Lambda 2
  
                    B_1_1_field = besseli(1,sqrt(s/CapitalLambda_1)*r_in_field)/sqrt(s/CapitalLambda_1); % B_1_1: stand for Besseli function, 1 order, contains Lambda 1
                         
                    B_1_2_field = besseli(1,sqrt(s/CapitalLambda_2)*r_in_field)/sqrt(s/CapitalLambda_2); % B_1   _1: stand for Besseli function, 1 order, contains Lambda 2
 
    
                   
                 
                     tilde_TT = Transition_P(2,1)*(c_1*B_0_1_field-H_1*fs*CapitalLambda_1) + Transition_P(2,2)*(c_2*B_0_2_field-H_2*fs*CapitalLambda_2);   %% without s     
                        
                    
                       
                     tilde_pp = Transition_P(1,1)*(c_1*B_0_1_field-H_1*fs*CapitalLambda_1) + Transition_P(1,2)*(c_2*B_0_2_field-H_2*fs*CapitalLambda_2);   %% without s 
                   
    
           
   
%% Sigma_rr_first transfer back to time domain

                 a_77 = (1/r_in_field)* B_1_1_field *(AA1*Transition_P(1,1)+AA2*Transition_P(2,1))+B_0_1_field*(AA3*Transition_P(1,1)+AA4*Transition_P(2,1));
    
                 a_88 = (1/r_in_field)* B_1_2_field *(AA1*Transition_P(1,2)+AA2*Transition_P(2,2))+B_0_2_field*(AA3*Transition_P(1,2)+AA4*Transition_P(2,2));
                 
                 
                 a_99 = -AA5+(0.5*AA1*((Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2))+ 0.5*AA2*((Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2))+ AA3*((Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2))+...
    AA4*((Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2)));


%% To calculate the radial displacement

                 a_77_ur = B_1_1_field*(Transition_P(1,1)*eta/G+Transition_P(2,1)*eta_d/G);
    
                 a_88_ur = B_1_2_field*(Transition_P(1,2)*eta/G+Transition_P(2,2)*eta_d/G);
    
                 a_99_ur = 0.5*r_in_field*(eta*(Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2)/G+eta_d*(Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2)/G); % original but results does not make sense
%% Iteration of s and Cjj

        tilde_TT_total = tilde_TT_total + tilde_TT*Cjj;
  
        tilde_pp_total = tilde_pp_total + tilde_pp*Cjj;
        
        tilde_sigma_rr_part_3_total = tilde_sigma_rr_part_3_total + (a_77*c_1+a_88*c_2-a_99*fs)*Cjj; %% sigma_rr
              
        tilde_sigma_rr_part_4_total = tilde_sigma_rr_part_4_total + (a_77_ur*c_1+a_88_ur*c_2-a_99_ur*fs)*Cjj; %% u_rr
      end 
%% Transfer back to time domain

    RR_T(r) = tilde_TT_total*(log(2)/t0);
    RR_p(r) = tilde_pp_total*(log(2)/t0);
    
    sigma_rr(r)= tilde_sigma_rr_part_3_total*(log(2)/t0);
   
    u_rr(r) = tilde_sigma_rr_part_4_total*(log(2)/t0);
    
    RRR_T(mmm,r) = RR_T(r); 
    
    RRR_p_1(mmm,r) = sigma_rr(r)-RR_p(r);
    
    RRR_rr(mmm,r) = u_rr(r);

    end
  %cd 'C:\Users\yul184\Dropbox\P&A Project\10. 2022\4\Sum_traction'
  %save('mode_1.mat','RRR_p')
end 


%% Mode_2 Loading: temperature loading
for mmm = 1:1:nn
    
    t0 = t(mmm);
  
    
    for r = 1:100
    
    % Temperature in Laplace domain
    tilde_TT_total = 0;
    % Pore pressure in Laplace domain
    tilde_pp_total  = 0;
    % Radial stress in Laplace domain
  
    tilde_sigma_rr_part_1_total = 0;
    
    tilde_sigma_rr_part_2_total = 0;
    
    tilde_sigma_rr_part_3_total = 0;
    
    tilde_sigma_rr_part_4_total = 0;
    
    tilde_sigma_rr_part_5_total = 0;
    
 
      for j = 1:n
    
        kk = floor((j+1)/2);
        mi = min(j,n/2);
    
        Cjj = 0;
    
            for i = kk:mi
      
                    Cj = ((i^(n/2))*factorial(2*i))/(factorial(n/2-i)*factorial(i)*factorial(i-1)*factorial(j-i)*factorial(2*i-j));
                    Cjj = Cjj + Cj;
            end 
    
                    Cjj = Cjj *((-1)^(j+n/2));
    
                     s = j*log(2)/t0;
    
    
    %% obtain c1 & c2 & fs based on different boundary conditions
     
      
      H_matrix = Transition_P\((B_matrix)\[B_5;B_10]);
     
     
     H_1 = H_matrix(1);
     
 
     H_2 = H_matrix(2);
     
 
     r_at_boundary = 100/250;
     
     B_0_1_boundary = besseli(0,sqrt(s/CapitalLambda_1)*r_at_boundary); % B_0_1: stand for Besseli function, 0 order, contains Lambda 1

     B_0_2_boundary = besseli(0,sqrt(s/CapitalLambda_2)*r_at_boundary);   % B_0_1: stand for Besseli function, 0 order, contains Lambda 2
  
     B_1_1_boundary = besseli(1,sqrt(s/CapitalLambda_1)*r_at_boundary)/sqrt(s/CapitalLambda_1); % B_1_1: stand for Besseli function, 1 order, contains Lambda 1
                         
     B_1_2_boundary = besseli(1,sqrt(s/CapitalLambda_2)*r_at_boundary)/sqrt(s/CapitalLambda_2); % B_1   _1: stand for Besseli function, 1 order, contains Lambda 2
 
 
 
 
    a_1 = Transition_P(1,1)*B_0_1_boundary;
 
    a_2 = Transition_P(1,2)*B_0_2_boundary;
 
    a_3 = Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2;
   
 
 
    a_4 = Transition_P(2,1)*B_0_1_boundary;
 
    a_5 = Transition_P(2,2)*B_0_2_boundary;
 
    a_6 = Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2;
    
    a_7 = (1/r_at_boundary)*B_1_1_boundary*(AA1*Transition_P(1,1)+AA2*Transition_P(2,1))+B_0_1_boundary*(AA3*Transition_P(1,1)+AA4*Transition_P(2,1));
    
    a_8 = (1/r_at_boundary)*B_1_2_boundary*(AA1*Transition_P(1,2)+AA2*Transition_P(2,2))+B_0_2_boundary*(AA3*Transition_P(1,2)+AA4*Transition_P(2,2));
   
    a_9 = -AA5+(0.5*AA1*((Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2))+ 0.5*AA2*((Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2))+ AA3*((Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2))+...
    AA4*((Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2)));
    
    for c_i = 1:6*nn
        
        a7_check_matrix(c_i,r) = a_7;
        a8_check_matrix(c_i,r) = a_8;
        a9_check_matrix(c_i,r) = a_9;
        
    end
 
   
    Coefficient_Matrix = [a_1 a_2 -a_3; a_4 a_5 -a_6; a_7 a_8 -a_9];
    
     
 
    
    % thermal, pore and radial here are three different boundary conditions: 
    % x: thermal conditions;
    % y: pore pressure conditions;
    % z: radial stress boundary conditions;
    
    x = Loading_mode_2_input(1)/s;
    y = Loading_mode_2_input(2)/s;
    z = Loading_mode_2_input(3)/s;
 
    Boundary_Matrix = [x;y;z];
    
    
 
    Coefficient_result_matrix = Coefficient_Matrix\Boundary_Matrix;
    
    c_1 = Coefficient_result_matrix(1);
 
    c_2 = Coefficient_result_matrix(2);
 
    fs = Coefficient_result_matrix(3);
    
    
    
    for c_i = 1:6*nn
        
        c1_check_matrix(c_i,r) = c_1;
        c2_check_matrix(c_i,r) = c_2;
        fs_check_matrix(c_i,r) = fs;
        
    end

    
   
   
   

                     
 %% In Laplace domain, obtian the p_tilde, T_tilde and sigma_rr_tilde
                    r_in_field = r/250;
                   
                    
                    B_0_1_field = besseli(0,sqrt(s/CapitalLambda_1)*r_in_field); % B_0_1: stand for Besseli function, 0 order, contains Lambda 1

                    B_0_2_field = besseli(0,sqrt(s/CapitalLambda_2)*r_in_field);   % B_0_1: stand for Besseli function, 0 order, contains Lambda 2
  
                    B_1_1_field = besseli(1,sqrt(s/CapitalLambda_1)*r_in_field)/sqrt(s/CapitalLambda_1); % B_1_1: stand for Besseli function, 1 order, contains Lambda 1
                         
                    B_1_2_field = besseli(1,sqrt(s/CapitalLambda_2)*r_in_field)/sqrt(s/CapitalLambda_2); % B_1   _1: stand for Besseli function, 1 order, contains Lambda 2
 
    
                   
                 
                     tilde_TT = Transition_P(2,1)*(c_1*B_0_1_field-H_1*fs*CapitalLambda_1) + Transition_P(2,2)*(c_2*B_0_2_field-H_2*fs*CapitalLambda_2);   %% without s     
                        
                    
                       
                     tilde_pp = Transition_P(1,1)*(c_1*B_0_1_field-H_1*fs*CapitalLambda_1) + Transition_P(1,2)*(c_2*B_0_2_field-H_2*fs*CapitalLambda_2);   %% without s 
                   
    
           
   
%% Sigma_rr_first transfer back to time domain

                 a_77 = (1/r_in_field)* B_1_1_field *(AA1*Transition_P(1,1)+AA2*Transition_P(2,1))+B_0_1_field*(AA3*Transition_P(1,1)+AA4*Transition_P(2,1));
    
                 a_88 = (1/r_in_field)* B_1_2_field *(AA1*Transition_P(1,2)+AA2*Transition_P(2,2))+B_0_2_field*(AA3*Transition_P(1,2)+AA4*Transition_P(2,2));
                 
                 
                 a_99 = -AA5+(0.5*AA1*((Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2))+ 0.5*AA2*((Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2))+ AA3*((Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2))+...
    AA4*((Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2)));


%% To calculate the radial displacement

                 a_77_ur = B_1_1_field*(Transition_P(1,1)*eta/G+Transition_P(2,1)*eta_d/G);
    
                 a_88_ur = B_1_2_field*(Transition_P(1,2)*eta/G+Transition_P(2,2)*eta_d/G);
    
                 a_99_ur = 0.5*r_in_field*(eta*(Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2)/G+eta_d*(Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2)/G); % original but results does not make sense
%% Iteration of s and Cjj

        tilde_TT_total = tilde_TT_total + tilde_TT*Cjj;
  
        tilde_pp_total = tilde_pp_total + tilde_pp*Cjj;
        
        tilde_sigma_rr_part_3_total = tilde_sigma_rr_part_3_total + (a_77*c_1+a_88*c_2-a_99*fs)*Cjj; %% sigma_rr
              
        tilde_sigma_rr_part_4_total = tilde_sigma_rr_part_4_total + (a_77_ur*c_1+a_88_ur*c_2-a_99_ur*fs)*Cjj; %% u_rr
      end 
%% Transfer back to time domain

    RR_T(r) = tilde_TT_total*(log(2)/t0);
    RR_p(r) = tilde_pp_total*(log(2)/t0);
    
    sigma_rr(r)= tilde_sigma_rr_part_3_total*(log(2)/t0);
   
    u_rr(r) = tilde_sigma_rr_part_4_total*(log(2)/t0);
    
    RRR_T(mmm,r) = RR_T(r); 
    
    RRR_p_2(mmm,r) = sigma_rr(r)+RR_p(r);
    
    RRR_rr(mmm,r) = u_rr(r);

    end
  %cd 'C:\Users\yul184\Dropbox\P&A Project\10. 2022\4\Sum_traction'
  %save('mode_2.mat','RRR_p')
end 


%% Superposition
%cd 'C:\Users\yul184\Dropbox\P&A Project\10. 2022\4\Sum_traction'
%Mode_1 = load('mode_1.mat');
%Mode_2 = load('mode_2.mat');
%Mode_3 = load('mode_3.mat');

%Mode1_value = Mode_1.RRR_p;
%Mode2_value = Mode_2.RRR_p;
%Mode3_value = Mode_3.RRR_p;
RRR_p_3 = Loading_mode_3_input(3)*ones(16,100);

Total = RRR_p_1 + RRR_p_2 + RRR_p_3;

%% Picking the max point from the curve (tensile failure)
[value,location]= max(Total,[],2);
Max_tensile = max(value);
%if Max_tensile <=0
  %  Max_tensile=0;
%end
input_matrix_this_time{end+1}=Max_tensile;

Over_all_data{end+1}=input_matrix_this_time;
input_matrix_this_time=[];

%% delete the temperotry file
%cd 'C:\Users\yul184\Dropbox\P&A Project\10. 2022\4\Sum_traction'
%delete mode_1.mat;
%delete mode_2.mat;
%delete mode_3.mat;
clc
             
%Sample_iterations = Sample_iterations+1;         

        
end       
        


