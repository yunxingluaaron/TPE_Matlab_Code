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
   Input = xlsread('para.xlsx','Sheet1');

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
    RR_p(r) = tilde_pp_total*(log(2)/t0)-Input(1,3);
    
    sigma_rr(r)= tilde_sigma_rr_part_3_total*(log(2)/t0);
   
    u_rr(r) = tilde_sigma_rr_part_4_total*(log(2)/t0);
    
    RRR_T(mmm,r) = RR_T(r); 
    
    RRR_p_1(mmm,r) = sigma_rr(r)+ RR_p(r);
    % Mode_1 Loading: pore pressure loading
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
    
    RRR_p_2(mmm,r) = sigma_rr(r)+ RR_p(r);
    
    % Mode_2 Loading: temperature loading
    
    RRR_rr(mmm,r) = u_rr(r);

    end
  %cd 'C:\Users\yul184\Dropbox\P&A Project\10. 2022\4\Sum_traction'
  %save('mode_2.mat','RRR_p')
end 

%% Mode_3 Loading: iso-tropic far field stress loading
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
    
    x = Loading_mode_3_input(1)/s;
    y = Loading_mode_3_input(2)/s;
    z = Loading_mode_3_input(3)/s;
 
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
    
    RRR_p_3(mmm,r) = sigma_rr(r)-RR_p(r);
    
    RRR_rr(mmm,r) = u_rr(r);
    
    end
  %cd 'C:\Users\yul184\Dropbox\P&A Project\10. 2022\4\Sum_traction'
  %save('mode_3.mat','RRR_p')
  
end 


%% Superposition
%cd 'C:\Users\yul184\Dropbox\P&A Project\10. 2022\4\Sum_traction'
%Mode_1 = load('mode_1.mat');
%Mode_2 = load('mode_2.mat');
%Mode_3 = load('mode_3.mat');

%Mode1_value = Mode_1.RRR_p;
%Mode2_value = Mode_2.RRR_p;
%Mode3_value = Mode_3.RRR_p;

%Total = Mode1_value + Mode2_value + Mode3_value;
Total = RRR_p_1+RRR_p_2+RRR_p_3;
%% Picking the max point from the curve


%% Black to grey
lightBLUE = [  0.85980392156863,   0.850980392156863 ,  0.850980392156863];
darkBLUE = [0  ,   0 ,    0];
%% Red Group  
%lightBLUE = [  0.988235294117647,   0.890196078431372,   0.890196078431372];
%darkBLUE = [0.600000000000000,   0.019607843137255,   0.019607843137255]; 
%% Green
%lightBLUE = [ 0.768627450980392 ,  0.929411764705882 ,  0.447058823529412];
%darkBLUE = [ 0.113725490196078 ,  0.380392156862745 ,  0.019607843137255]; 
%% Blue
%lightBLUE = [0.356862745098039,0.811764705882353,0.956862745098039];
%darkBLUE = [0.0196078431372549,0.0745098039215686,0.670588235294118]; 
   
   
blueGRADIENTflexible = @(i_color,N_color) lightBLUE + (darkBLUE-lightBLUE)*((i_color-1)/(N_color-1));
    
Rx = linspace(0,0.25,100);
for n=1:nn
    hold on;
    N = Total(n,:);
    plot(Rx,N,'Color', blueGRADIENTflexible(n,11));
end 

plot([min(xlim()),max(xlim())],[0,0], 'k--')