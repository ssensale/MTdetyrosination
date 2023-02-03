function new_position = cell_rebound(old_position, newposition, cell_rad, enzyme_size)
   
    vector = newposition - old_position;
    vector_mag = (vector(1) ^ 2 + vector(2) ^ 2 + vector(3) ^ 2) ^ 0.5;
    LINE   = [old_position(1) old_position(2) old_position(3) vector(1) vector(2) vector(3)];
    SPHERE = [0 0 0 cell_rad];
    GC = intersectLineSphere(LINE, SPHERE);
    vectorA=GC(1,:)-old_position;
    vectorB=GC(2,:)-old_position;
    if vectorA(1)*vector(1)+vectorA(2)*vector(2)+vectorA(3)*vector(3)<0
        move=vector-vectorB;
        new_position=GC(2,:)-move;
    else
        move=vector-vectorA;
        new_position=GC(1,:)-move;        
    end;
    
    if new_position(1)^2+new_position(2)^2+new_position(3)^2>(cell_rad-enzyme_size)^2
        new_position=old_position;
    end
end