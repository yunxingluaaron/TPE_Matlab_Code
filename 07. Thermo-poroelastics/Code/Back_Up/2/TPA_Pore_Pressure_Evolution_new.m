format long

    %t_time = t./86400;
    
    size_time = size(t);
    
    size_time = size_time(1);
    
    num_time = xlsread('Book1.xlsx','Sheet1');
    
    
    n_time_point = 450;
    
    
    num_time_second = num_time(1:n_time_point,5);
    
    num_time_day =num_time(1:n_time_point,6);
    
    num_time_log =num_time(1:n_time_point,7);
    
   
    
    
    n_time_exponent_x_point = linspace(0.1,11,n_time_point);
    
    n_time_exponent_x_point = transpose(n_time_exponent_x_point);
    
    n_time_exponent_log = zeros(n_time_point,1);
    
    n_time_exponent_second = zeros(n_time_point,1);
    
    for p = 1:n_time_point
    %n_time_exponent_log(p,1) = 0.3864*n_time_exponent_x_point(p,1)-4.6003;
    %n_time_exponent_log(p,1) = 0.4113*n_time_exponent_x_point(p,1)-4.7079;
    n_time_exponent_log(p,1) = 0.0147*n_time_exponent_x_point(p,1)-0.8976;
    
    end
    
   
    for p = 1:n_time_point
    n_time_exponent_second(p,1) = 10^n_time_exponent_log(p,1)*86400;
    
    end
    
    pore_pressure_array = zeros(n_time_point,1);
     
    temperature_array = zeros(n_time_point,1);
    
    
    n_time_temperature = zeros(n_time_point,1);
    
    
    n_time_oore_pressure = zeros(n_time_point,1);
    
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
   
lightBLUE = [  0.850980392156863,   0.850980392156863 ,  0.850980392156863];
darkBLUE = [0  ,   0 ,    0];
   
   
   blueGRADIENTflexible = @(i_color,N_color) lightBLUE + (darkBLUE-lightBLUE)*((i_color-1)/(N_color-1));
    

for mmm = 1:1:n_time_point
    
    t0 = num_time_second (mmm,1);
  
    
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
     
 
     r_at_boundary = 100/250;
     
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
    % x: pore pressure conditions;
    % y: thermal conditions;
    % z: radial stress boundary conditions;
    
    x = 2.9e7/s;
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
    
    RRR_T(mmm,r) = RR_T(r); 
    
    RRR_p(mmm,r) = sigma_rr(r)+ RR_p(r);
    
    RRR_rr(mmm,r) = RR_p(r);
    
    pore_pressure_array(mmm,1) = RR_p(1);
    
    
    temperature_array(mmm,1) = RR_T(1);
    
    
    n_time_temperature(mmm,1) = t0;
    
    
    n_time_pore_pressure(mmm,1) = t0;

    end 
%% Calculate the sigma_rr
hold on 

%plot(n_time_exponent_second./86400,pore_pressure_array);
%semilogx(n_time_exponent_second,pore_pressure_array)
%plot(Rx,RR_T, 'Color', blueGRADIENTflexible(mmm,nn));
%plot(Rx,-RR_p, 'Color', blueGRADIENTflexible(mmm,nn));
%plot(t,pore_pressure_array, 'Color', blueGRADIENTflexible(mmm,nn));
%plot(Rx,sigma_rr, 'Color', blueGRADIENTflexible(mmm,nn));
%plot(Rx,sigma_rr+RR_p, 'Color', blueGRADIENTflexible(mmm,nn));
%plot(Rx,u_rr, 'Color', blueGRADIENTflexible(mmm,nn));
%title('Bessel with Lamda')
%leg1 = legend('t=17','t=18','t=19','t=20','t=21','t=22','t=23','t=24','t=25','t=26','t=27');

end 
%axis([0.000001 1])
%n_time_exponent_day = n_time_exponent_second./86400;

n_time_exponent_day_p = n_time_pore_pressure./86400;
n_time_exponent_day_T = n_time_temperature./86400;


plot(num_time_day,pore_pressure_array,'k','LineWidth',2);

%plot(n_time_exponent_day_T,temperature_array,'k','LineWidth',2);

grid on
set(gca,'xscale','log')
%ylim([0 3e7])
%xlim([0.5e-5 1e0])
%ylim([0 6000])


ly=ylabel('\Delta p  (Pa)');
lx=xlabel('Time(Day)');
set(gca, 'linewidth', 1.1)
set(gca, 'FontSize', 15)
set(lx,'FontSize',24)
set(ly,'FontSize',24)
box on;