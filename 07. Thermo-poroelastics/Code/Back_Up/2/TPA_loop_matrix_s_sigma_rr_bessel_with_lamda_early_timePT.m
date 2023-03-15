
%% input

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
for iii = 1:3
%% baisc coefficient
syms cn;
cn =1;
num = xlsread('para.xlsx','Sheet1');
a = num(1,cn);
tau_0 = num(2,cn); %Unit: K

beta_d= num(3,cn); %https://https://www.sciencedirect.com/science/article/pii/S0266352X20302925
beta_v= num(4,cn);  %https://https://www.sciencedirect.com/science/article/pii/S0266352X20302925

K = num(5,cn); %Unit: N/m^2, P106 Table 3.1 
B = num(6,cn); % Default
G = num(7,cn); % Unit: N/m^2, P106 Table 3.1
v = num(8,cn);  % Default

c_d = num(13,cn); %Unit: J/m^3 K,  P633, Table 11.3

t = num(16:26,cn);
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


k_tp_matrix =[1e-5, 1e-6,1e-7];


k_T = 10;
k_Tp = k_tp_matrix(iii);
k = 1e-18;
k_pT= 1e-14;



visco_sity = 1e-4;
kappa = k/visco_sity; % the devidator is the viscosity
kappa_pT = k_pT;
%kappa_Tp = k_Tp/(tau_0);
%kappa_T = k_T/(tau_0);
kappa_T = k_T;
kappa_Tp = k_Tp;


%B_1 = 1-(alpha_d*eta/(m_d*G)+beta_e/m_d)*(M*k_pT*G)/(a*M*eta+G)+kappa_T;

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

%CapitalLambda_1 = 30;

%CapitalLambda_2 = 3;
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

%Input is over 
%%



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
    RRR_p = zeros(nn,100);
    RRR_rr = zeros(nn,100);
    
    c1_check_matrix = zeros(6*nn,100);
    c2_check_matrix = zeros(6*nn,100);
    fs_check_matrix = zeros(6*nn,100);
    
    a7_check_matrix = zeros(6*nn,100);
    a8_check_matrix = zeros(6*nn,100);
    a9_check_matrix = zeros(6*nn,100);
    
    %AA1 = 2*G*v*eta/(G*(1-2*v))-(6*G*v+2*G)*eta/(G*(1-2*v));
    
    %AA2 = 2*G*v*eta_d/(G*(1-2*v))-(6*G*v+2*G)*eta_d/(G*(1-2*v));
    
    %AA3 = (6*G*v+2*G)*eta/(G*(1-2*v))-a;
    
    %AA4 = (6*G*v+2*G)*eta_d/(G*(1-2*v))-alpha_d;
    
    %AA5 = (8*G*v+2*G)/(1-2*v);
    
    
    AA1 = 2*G*v*eta/(G*(1-2*v))-(2*G-2*G*v)*eta/(G*(1-2*v));
    
    AA2 = 2*G*v*eta_d/(G*(1-2*v))-(2*G-2*G*v)*eta/(G*(1-2*v));
    
    AA3 = (-2*G*v+2*G)*eta/(G*(1-2*v))-a;
    
    AA4 = (-2*G*v+2*G)*eta_d/(G*(1-2*v))-alpha_d;
    
    AA5 = (2*G)/(1-2*v);
    
   %AA5_1 = (2*G)/(1-2*v);
    
   %AA5_2 = 2*(2*G-2*G*v)/(1-2*v);
   
   
   %% color map setting 
   
   lightBLUE = [0.356862745098039,0.811764705882353,0.956862745098039];
darkBLUE = [0.0196078431372549,0.0745098039215686,0.670588235294118];
   
   
   blueGRADIENTflexible = @(i_color,N_color) lightBLUE + (darkBLUE-lightBLUE)*((i_color-1)/(N_color-1));
    

    
    t0 = 60;
  
    
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
     
      %H_matrix = Transitio_P_inverse*(pinv(B_matrix)*[B_5;B_10]);
      H_matrix = Transition_P\((B_matrix)\[B_5;B_10]);
     
     
     H_1 = H_matrix(1);
     
 
     H_2 = H_matrix(2);
     
 
     r_at_boundary = 100/400;
     
     B_0_1_boundary = besseli(0,sqrt(s/CapitalLambda_1)*r_at_boundary); % B_0_1: stand for Besseli function, 0 order, contains Lambda 1

     B_0_2_boundary = besseli(0,sqrt(s/CapitalLambda_2)*r_at_boundary);   % B_0_1: stand for Besseli function, 0 order, contains Lambda 2
  
     B_1_1_boundary = besseli(1,sqrt(s/CapitalLambda_1)*r_at_boundary)/sqrt(s/CapitalLambda_1); % B_1_1: stand for Besseli function, 1 order, contains Lambda 1
                         
     B_1_2_boundary = besseli(1,sqrt(s/CapitalLambda_2)*r_at_boundary)/sqrt(s/CapitalLambda_2); % B_1   _1: stand for Besseli function, 1 order, contains Lambda 2
 
 
     %B_0_1_boundary = besseli(0,sqrt(s)*CapitalLambda_1*r_at_boundary); % B_0_1: stand for Besseli function, 0 order, contains Lambda 1

     %B_0_2_boundary = besseli(0,sqrt(s)*CapitalLambda_2*r_at_boundary);   % B_0_1: stand for Besseli function, 0 order, contains Lambda 2
  
     %B_1_1_boundary = besseli(1,sqrt(s)*CapitalLambda_1*r_at_boundary)/(CapitalLambda_1*sqrt(s)); % B_1_1: stand for Besseli function, 1 order, contains Lambda 1
                         
     %B_1_2_boundary = besseli(1,sqrt(s)*CapitalLambda_2*r_at_boundary)/(CapitalLambda_2*sqrt(s)); % B_1   _1: stand for Besseli function, 1 order, contains Lambda 2
                     
 
    a_1 = Transition_P(1,1)*B_0_1_boundary;
 
    a_2 = Transition_P(1,2)*B_0_2_boundary;
 
    %a_3 = Transition_P(1,1)*H_1*s+Transition_P(1,2)*H_2*s;
    a_3 = Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2;
   
 
 
    a_4 = Transition_P(2,1)*B_0_1_boundary;
 
    a_5 = Transition_P(2,2)*B_0_2_boundary;
 
    %a_6 = Transition_P(2,1)*H_1*s+Transition_P(2,2)*H_2*s;
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
 
     %Coefficient_Matrix = [a_1 a_2 a_3; a_4 a_5 a_6; a_7 a_8 a_9];
    Coefficient_Matrix = [a_1 a_2 -a_3; a_4 a_5 -a_6; a_7 a_8 -a_9];
    
     
 
    
    % thermal, pore and radial here are three different boundary conditions: 
    % x: thermal conditions;
    % y: pore pressure conditions;
    % z: radial stress boundary conditions;
    
    x = -2.9e7/s;
    y = 0;
    z = 0;
 
    Boundary_Matrix = [x;y;z];
    
    
    
    %Coefficient_result_matrix = inv(Coefficient_Matrix)*Boundary_Matrix;
 
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
                    r_in_field = r/400;
                   
                    
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
   
                 
                 
                 %a_99_ur = 0.5*r_in_field*s*(eta*(Transition_P(2,1)*H_1+Transition_P(2,2)*H_2)/G+eta_d*(Transition_P(1,1)*H_1+Transition_P(1,2)*H_2)/G); % Compraed with original, the fisrt term is removed.
                 
              
             
                 
                 %a_99_ur =-r_in_field + 0.5*r_in_field*(eta*(Transition_P(1,1)*H_1*CapitalLambda_1+Transition_P(1,2)*H_2*CapitalLambda_2)/G+eta_d*(Transition_P(2,1)*H_1*CapitalLambda_1+Transition_P(2,2)*H_2*CapitalLambda_2)/G);
                 
                 
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
    
    %RRR_T(mmm,r) = RR_T(r); 
    
    %RRR_p(mmm,r) = sigma_rr(r)-RR_p(r);
    
    %RRR_rr(mmm,r) = u_rr(r);

    end 
%% Calculate the sigma_rr
hold on 
lightBLUE = [  0.750980392156863,   0.750980392156863 ,  0.750980392156863];
darkBLUE = [0  ,   0 ,    0];
   
   
blueGRADIENTflexible = @(i_color,N_color) lightBLUE + (darkBLUE-lightBLUE)*((i_color-1)/(N_color-1));
Rx = linspace(0,0.25,100);

%plot(Rx,RR_T, 'Color', blueGRADIENTflexible(mmm,nn));
%plot(Rx,RR_p, 'Color', blueGRADIENTflexible(mmm,nn));
%plot(Rx,RR_p,'Color', blueGRADIENTflexible(iii,3));
plot(Rx,RR_T,'Color', blueGRADIENTflexible(iii,3),'LineWidth',2);
%plot(Rx,RR_T);
%plot(Rx,sigma_rr, 'Color', blueGRADIENTflexible(mmm,nn));
%plot(Rx,sigma_rr-RR_p, 'Color', blueGRADIENTflexible(mmm,nn));
%plot(Rx,u_rr, 'Color', blueGRADIENTflexible(mmm,nn));
%title('Bessel with Lamda')
%leg1 = legend('k_Pt=-9','t=18','t=19','t=20','t=21','t=22','t=23','t=24','t=25','t=26','t=27');

%leg1 = legend('k_tp=-5','k_tp=-6','k_tp=-7','k_tp=-0', 'double zero');

%ylim([0 20e6])

end

%leg1 = legend('k_{pt}/k=10^{1}','k_{pt}/k=10^{3}','k_{pt}/k=10^{5}');
leg1 = legend('k_{T}/k_{tp}=10^{1}','k_{T}/k_{tp}=10^{3}','k_{T}/k_{tp}=10^{5}');
rect = [0.25, 0.25, .25, .25];
set(leg1, 'Position', rect)
leg1.FontSize = 20;

ylim([0 1.2])
box on;
