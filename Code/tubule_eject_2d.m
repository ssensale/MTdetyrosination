function new_position = tubule_eject_2d(old_position, new_position, tubule_size, enzyme_size, cylinder1, cylinder2)
    distance = ((new_position(1) - old_position(1)) ^ 2 + ...
        (new_position(2) - old_position(2)) ^ 2 + ...
        (new_position(3) - old_position(3)) ^ 2) ^ 0.5;
    
    launch_point = old_position; 
    
    vector = [rand(), rand(), rand()];
    vector = distance * vector/(vector(1)^2+vector(2)^2+vector(3)^2)^0.5;
    
    new_position = launch_point + vector;
    
    % Evaluate to see if the point ended up in the cylinder
    cylinder = cat(1, cylinder1', cylinder2');
    distancia = dcylindersphere(cylinder', tubule_size, new_position', enzyme_size);
    
    % If so, flip vector and send the enzyme in the other direction
    if distancia < tubule_size
        vector = -1 .* vector;
        new_position = launch_point + vector;
    end

end