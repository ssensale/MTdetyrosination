function cylinder_location = onetubule_setup(center_rad, cell_rad, xcutoff, count, numtubesleft)
    ratio = xcutoff / cell_rad;
    if count <= numtubesleft
        x_begin = (center_rad + (ratio * center_rad)) * rand() - center_rad; 
        y_begin = rand() * sign(2 * rand() - 1) * (center_rad ^ 2 - x_begin ^ 2) ^ 0.5;
        z_begin = sign(2 * rand() - 1) * (center_rad ^ 2 - x_begin ^ 2 - y_begin ^ 2) ^ 0.5;

        beginpoint = cat(1, x_begin, y_begin, z_begin);

        x_end = (cell_rad + (ratio * cell_rad)) * rand() - cell_rad;
        checked = 0;
        
        while ~checked
            if sign(x_end) == sign(x_begin)
                checked = 1;
            else
                x_end = (cell_rad + (ratio * cell_rad)) * rand() - cell_rad;
            end
        end

        y_end = (2 .* rand(1, 1) - 1) * (cell_rad ^ 2 - x_end ^ 2) ^ 0.5;

        if sign(y_end) ~= sign(y_begin)
            y_end = -y_end;
        end

        z_end = sign(2 * rand() - 1) * (cell_rad ^ 2 - x_end ^ 2 - y_end ^ 2) ^ 0.5;

        if sign(z_end) ~= sign(z_begin)
            z_end = -z_end;
        end
        
    else
        x_begin = (center_rad - (ratio * center_rad)) * rand() + (ratio * center_rad); 
        y_begin = rand() * sign(2 * rand() - 1) * (center_rad ^ 2 - x_begin ^ 2) ^ 0.5;
        z_begin = sign(2 * rand() - 1) * (center_rad ^ 2 - x_begin ^ 2 - y_begin ^ 2) ^ 0.5;

        beginpoint = cat(1, x_begin, y_begin, z_begin);

        x_end = (cell_rad - (ratio * cell_rad)) * rand() + (ratio * cell_rad);
        checked = 0;
        
        while ~checked
            if sign(x_end) == sign(x_begin)
                checked = 1;
            else
                x_end = (cell_rad - (ratio * cell_rad)) * rand() + (ratio * cell_rad);
            end
        end

        y_end = (2 .* rand(1, 1) - 1) * (cell_rad ^ 2 - x_end ^ 2) ^ 0.5;

        if sign(y_end) ~= sign(y_begin)
            y_end = -y_end;
        end

        z_end = sign(2 * rand() - 1) * (cell_rad ^ 2 - x_end ^ 2 - y_end ^ 2) ^ 0.5;

        if sign(z_end) ~= sign(z_begin)
            z_end = -z_end;
        end
    end
    
    endpoint = cat(1, x_end, y_end, z_end);
    cylinder_location([1,2,3])=beginpoint;
    cylinder_location([4,5,6])=endpoint;
end