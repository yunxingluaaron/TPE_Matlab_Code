%%
Band_interval = 200;

x_cumulative_tensile = linspace(0,1,Band_interval);

x_cumulative_compressive = linspace(0,1,Band_interval);

%%

Zeta_1_tensile_array_sort = sort(Zeta_1_tensile_array);

Zeta_1_compressive_array_sort = sort(Zeta_1_compressive_array);


%%

s_current_tensile = 0;

s_current_compressive  = 0;

%%
Tensile_distribution = zeros(Band_interval,1);

Compressive_distribution = zeros(Band_interval,1);

for i = 1:Band_interval
        
        s_current_tensile = 0;
        
        for j = 1: Total_tensile_Number
            
            if  i < Band_interval-1 && Zeta_1_tensile_array_sort(j) >= x_cumulative_tensile(i) && Zeta_1_tensile_array_sort(j) <= x_cumulative_tensile(i+1)
        
                s_current_tensile = s_current_tensile+1;
             
            end   
         
        end 
        
        Tensile_distribution(i) = s_current_tensile;
        
        
        s_current_compressive = 0;
        
        for j = 1: Total_compressive_Number
            
            if  i < Band_interval-1 && Zeta_1_compressive_array_sort(j) >= x_cumulative_compressive(i) && Zeta_1_compressive_array_sort(j) <= x_cumulative_compressive(i+1)
        
                s_current_compressive = s_current_compressive+1;
             
            end   
         
        end 
        
         Compressive_distribution(i) = s_current_compressive;
    
end  

Total_distribution = zeros(Band_interval,1);

for i = 1:Band_interval

    Total_distribution(i) = Compressive_distribution(i)+ Tensile_distribution(i);

end

Total_out_put_to_Python = cat(2,Total_distribution,Tensile_distribution);


%plot(x_cumulative_tensile,tensile_accumulative_distritution,'r');

%hold on 

%plot(x_cumulative_compressive,compressive_accumulative_distritution,'color',[0 0.5 0]);


%xlim([0 3])

xlabel('Zeta One: HTND') 
ylabel('Cumulative Distribution Function (CDF)')
hA = get(gca);
%hA.XAxis.MinorTickValues = MinorXTicks;
hA.XAxis.MinorTick='on';
h=legend('Tensile','Compressive');
%h=legend('Tensile Distribution','Compressive Distribution');
set(h,'FontSize',20);
