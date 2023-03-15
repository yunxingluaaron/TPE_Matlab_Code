%X = Cell_Tensile_Cases;
X=Over_all_data;
%X = Cell_Compressive_Cases;
numer = 1642;
Matrix_To_Python = zeros(numer,19);
%Matrix_To_Python = zeros(1996,15);
%Matrix_To_Python = zeros(1246,15);
for i = 1:numer
    
    if length(X{i}) == 19
    
        cell_1 = X{i};

            for j = 1:19

                cell_2 = cell_1{j};

                Matrix_To_Python(i,j)= real(cell_2);
    
            end
    end 
end


for i = 1:numer
   
    if Matrix_To_Python(i,19) > 0
        
        Matrix_To_Python(i,19) = 1;
        
    elseif Matrix_To_Python(i,19) < 0
      
         Matrix_To_Python(i,19) = 0;
    end
        
end
cd 'C:\Users\yul184\Dropbox\P&A Project\15. Publications\07. Thermo-poroelastics\Code\Monte_Carlo_Sampling'

csvwrite('Total.csv',Matrix_To_Python)