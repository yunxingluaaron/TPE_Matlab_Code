hold on;


%total_number =104239;

syms Final_result

%total_number =104239;
total_number =length(Over_all_data);

column_1={0,0,0};
column_2={0,0,0};
column_3={0,0,0};
column_4={0,0,0};
column_5={0,0,0};
column_6={0,0,0};
column_7={0,0,0};
column_8={0,0,0};
column_9={0,0,0};
column_10={0,0,0};
column_11={0,0,0};


k_tensile_case=0;
j=0;


for i = 1:total_number
    
  
    current_cell = Over_all_data{1,i};
    
    
    if size(current_cell)==1
        j = j+1;
  
    else 
        Final_result = current_cell{19};
        %Final_result = current_cell{length(Over_all_data{2})};
    
           if Final_result > 0
        
                k_tensile_case = k_tensile_case+1;
   
           end
     end
    
   
end

%% Extration the Tensile Cases from The total Caese
Cell_Tensile_Cases = {};
Cell_Compressive_Cases = {};

for i = 1:total_number
    
  
    current_cell = Over_all_data{1,i};
    
    
    if size(current_cell)==1
        j = j+1;
  
    else 
        Final_result = current_cell{19};
        %Final_result = current_cell{length(Over_all_data{2})};
    
           if Final_result > 0
        
                Cell_Tensile_Cases{end+1}= current_cell;
                
           elseif  Final_result < 0
                
                Cell_Compressive_Cases{end+1}= current_cell;
               
           end
     end
    
   
end

Total_tensile_Number = size(Cell_Tensile_Cases);
Total_tensile_Number =Total_tensile_Number(2);

%%
Zeta_1_tensile_array = [];
Zeta_2_tensile_array = [];
for ii = 1:Total_tensile_Number
    
    Zeta_1_tensile = Cell_Tensile_Cases(1,ii);
    Zeta_1_tensile =[Zeta_1_tensile{:}];
    
    
    
    
    %Zeta_1_tensile_reverse = Zeta_1_tensile{17}*Zeta_1_tensile{1}*Zeta_1_tensile{2}/(Zeta_1_tensile{3}*Zeta_1_tensile{5});
    %Zeta_1_tensile_reverse = (Zeta_1_tensile{17})*(Zeta_1_tensile{2}+Zeta_1_tensile{3})/(Zeta_1_tensile{2}*(Zeta_1_tensile{5})*(Zeta_1_tensile{1}));
    Zeta_1_tensile_reverse = (Zeta_1_tensile{17})*(Zeta_1_tensile{2}+Zeta_1_tensile{3})*(Zeta_1_tensile{1})/(Zeta_1_tensile{2}*(Zeta_1_tensile{5}));
    %Zeta_1_tensile_reverse = (Zeta_1_tensile{17})*Zeta_1_tensile{3}*(Zeta_1_tensile{5})*(Zeta_1_tensile{1})/Zeta_1_tensile{2};
    %Zeta_1_tensile_reverse = Zeta_1_tensile{17};
    %Zeta_1_tensile_reverse =(Zeta_1_tensile{2}+Zeta_1_tensile{3})*Zeta_1_tensile{17}/(Zeta_1_tensile{2});%This
  
    
    %Zeta_1_tensile_reverse = Zeta_1_tensile{17}*Zeta_1_tensile{3}/Zeta_1_tensile{2};
    %Zeta_1_tensile_reverse = Zeta_1_tensile{17}*Zeta_1_tensile{3}*Zeta_1_tensile{1}*Zeta_1_tensile{5}/Zeta_1_tensile{2};
    
    
    A_T_1 = Zeta_1_tensile{1}/Zeta_1_tensile{4}-1/(Zeta_1_tensile{4}*Zeta_1_tensile{5});
    A_T_2 = Zeta_1_tensile{2}+Zeta_1_tensile{3}/2;
    
    pore_T_1 = Zeta_1_tensile{2}-Zeta_1_tensile{6}*A_T_2-2*Zeta_1_tensile{7}*Zeta_1_tensile{6}*A_T_2/(1-2*Zeta_1_tensile{7});
    pore_T_2 = 2*Zeta_1_tensile{7}*Zeta_1_tensile{6}*A_T_1/(1-2*Zeta_1_tensile{7})+A_T_1*Zeta_1_tensile{6}-Zeta_1_tensile{1};
    %Zeta_2_tensile_reverse = Zeta_1_tensile{7}*Zeta_1_tensile{6}/((Zeta_1_tensile{4}));
    %Zeta_2_tensile_reverse = Zeta_1_tensile{6}/((Zeta_1_tensile{4}));
    %Zeta_2_tensile_reverse = Zeta_1_tensile{4};
    %Zeta_2_tensile_reverse = Zeta_1_tensile{21}/5.1e7;
    %Zeta_2_tensile_reverse = (pore_T_1*40/pore_T_2)/5.1*e8-0.4;
    Zeta_2_tensile_reverse = (pore_T_1*40/pore_T_2)/(1.3e9)-0.1;
    
    
    %Zeta_1_tensile_array(end+1)=[(Zeta_1_tensile_reverse+0.04)];
    Zeta_1_tensile_array(end+1)=[(Zeta_1_tensile_reverse+0.02)];
    Zeta_2_tensile_array(end+1)=[(Zeta_2_tensile_reverse)];
   
end

%%
Total_Cell_Compressive_number = size(Cell_Compressive_Cases);

Total_compressive_Number =Total_Cell_Compressive_number(2);

Zeta_1_compressive_array = [];
Zeta_2_compressive_array = [];

for ii = 1:Total_compressive_Number
    
    Zeta_1_compressive = Cell_Compressive_Cases(1,ii);
    Zeta_1_compressive =[Zeta_1_compressive{:}];
   
    
   
    %Zeta_1__compressive_reverse = (Zeta_1_compressive{17})*Zeta_1_compressive{1}*Zeta_1_compressive{2}/(Zeta_1_compressive{3}*Zeta_1_compressive{5});
    %Zeta_1__compressive_reverse = (Zeta_1_compressive{17})*Zeta_1_compressive{3}+Zeta_1_compressive{2}/(Zeta_1_compressive{2}*Zeta_1_compressive{5}*Zeta_1_compressive{1});
    Zeta_1__compressive_reverse = (Zeta_1_compressive{17})*(Zeta_1_compressive{3}+Zeta_1_compressive{2})*Zeta_1_compressive{1}/(Zeta_1_compressive{2}*Zeta_1_compressive{5});
    %Zeta_1__compressive_reverse =  Zeta_1_compressive{17};
    %Zeta_1__compressive_reverse =  (Zeta_1_compressive{3}+Zeta_1_compressive{2})*Zeta_1_compressive{17}/(Zeta_1_compressive{2});
    %Zeta_1__compressive_reverse =  Zeta_1_compressive{17}*Zeta_1_compressive{3}/Zeta_1_compressive{2};
    %Zeta_1__compressive_reverse =  Zeta_1_compressive{17}*Zeta_1_compressive{3}*Zeta_1_compressive{1}*Zeta_1_compressive{5}/Zeta_1_compressive{2};
    
    
    A_C_1 = Zeta_1_compressive{1}/Zeta_1_compressive{4}-1/(Zeta_1_compressive{4}*Zeta_1_compressive{5});
    A_C_2 = Zeta_1_compressive{2}+Zeta_1_compressive{3}/2;
    
    pore_C_1 = Zeta_1_compressive{2}-Zeta_1_compressive{6}*A_C_2-2*Zeta_1_compressive{7}*Zeta_1_compressive{6}*A_C_2/(1-2*Zeta_1_compressive{7});
    pore_C_2 = 2*Zeta_1_compressive{7}*Zeta_1_compressive{6}*A_C_1/(1-2*Zeta_1_compressive{7})+A_C_1*Zeta_1_compressive{6}-Zeta_1_compressive{1};
    %Zeta_2_compressive_reverse = (Zeta_1_compressive{7}*Zeta_1_compressive{6})/(Zeta_1_compressive{4});
    %Zeta_2_compressive_reverse = (Zeta_1_compressive{6})/(Zeta_1_compressive{4});
    %Zeta_2_compressive_reverse = (Zeta_1_compressive{4});
    %Zeta_2_compressive_reverse = (Zeta_1_compressive{21}/5.1e7);
    %Zeta_2_compressive_reverse = (pore_C_1*40/ pore_C_2)/5.1e8-0.5;
    Zeta_2_compressive_reverse = (pore_C_1*40/ pore_C_2)/(1.3e9)-0.15;
    
    
    
    Zeta_1_compressive_array(end+1)=[(Zeta_1__compressive_reverse+0.005)];
    Zeta_2_compressive_array(end+1)=[(Zeta_2_compressive_reverse)];
   
    
    

end

%%
%Zeta_1_tensile_array = Zeta_1_tensile_array.*1e-48;
%Zeta_1_compressive_array = Zeta_1_compressive_array.*1e-48;


%TensileNN = unique(Zeta_1_tensile_array);
%CompressiveNN = unique(Zeta_1_compressive_array);

sz = 40;
Rx = linspace(0,Total_compressive_Number,Total_compressive_Number);
scatter(Zeta_1_compressive_array,Zeta_2_compressive_array,sz,'MarkerEdgeColor',[0 0.7 0],...
              'MarkerFaceColor',[0 0.7 0],...
              'LineWidth',1.5)
%scatter(Rx,Zeta_2_compressive_array,sz,'MarkerEdgeColor',[0 0.7 0],...
%              'MarkerFaceColor',[0 0.7 0],...
%              'LineWidth',1.5)
title(' ')



%figure


hold on
Ry = linspace(0,Total_tensile_Number,Total_tensile_Number);
%scatter(Ry,Zeta_2_tensile_array,sz,'MarkerEdgeColor',[1 0 0],...
%             'MarkerFaceColor',[1 0 0],...
%             'LineWidth',1.5)
scatter(Zeta_1_tensile_array,Zeta_2_tensile_array,sz,'MarkerEdgeColor',[1 0 0],...
             'MarkerFaceColor',[1 0 0],...
              'LineWidth',1.5)
%plot(TensileNN)
%title('Biot Coefficient')
%title('Thermal Expansion Coefficient Difference')
%title('Bulk Modulus')
%title('Poisons Ratio')
%title('Thermal Filtration')
%title('Thermal Conductivity')
%title('Permeability')
%title('Thermal Osmosis')
%title('Specific Heat Capacity')
%title('Zeta One')
%ylim([0 0.2e35])
%xlim([0 20])
ylim([0 2])

xlim([0 1])


%title('Zeta Two')

%xlabel('Zeta One: thermal diffusivity * Biot * Skempton * Thermal Expansion of Fluid / fluid diffusivity * Thermal Expansion of Solid')
%ylabel('Zeta Two: Shear Modulus*Poissons Ratio/Bluk Modulus')


