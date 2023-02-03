function enzyme_pos1 = enzyme_pos1(enzyme_size,center_rad, cell_size)
        x_center = (cell_size - enzyme_size) * (2 * rand() - 1);
        y_max = ((cell_size - enzyme_size) ^ 2 - x_center ^ 2) ^ 0.5;        
        y_center = sign(2 * rand() - 1) * y_max * rand();
        z_max = ((cell_size - enzyme_size) ^ 2 - x_center ^ 2 - y_center ^ 2) ^ 0.5;        
        z_center = sign(2 * rand() - 1) * z_max * rand();
        enzyme_pos1(1, 1) = x_center;
        enzyme_pos1(1, 2) = y_center;
        enzyme_pos1(1, 3) = z_center;      
        %Check no one falls inside center sphere
        while enzyme_pos1(1,1)^2+enzyme_pos1(1,2)^2+enzyme_pos1(1,3)^2<(center_rad+enzyme_size)^2
            x_center = (cell_size - enzyme_size) * (2 * rand() - 1);
            y_max = ((cell_size - enzyme_size) ^ 2 - x_center ^ 2) ^ 0.5;        
            y_center = sign(2 * rand() - 1) * y_max * rand();
            z_max = ((cell_size - enzyme_size) ^ 2 - x_center ^ 2 - y_center ^ 2) ^ 0.5;        
            z_center = sign(2 * rand() - 1) * z_max * rand();
            enzyme_pos1(1, 1) = x_center;
            enzyme_pos1(1, 2) = y_center;
            enzyme_pos1(1, 3) = z_center;   
        end;
end