function q = project_point_to_line_segment(A,B,p)
  % returns q the closest point to p on the line segment from A to B 

  % vector from A to B
  AB = (B-A);
  % squared distance from A to B
  AB_squared = AB(1)*AB(1)+AB(2)*AB(2)+AB(3)*AB(3);
  if(AB_squared == 0)
    % A and B are the same point
    q = A;
  else
    % vector from A to p
    Ap = (p-A);
        % Consider the line extending the segment, parameterized as A + t (B - A)
    % We find projection of point p onto the line. 
    % It falls where t = [(p-A) . (B-A)] / |B-A|^2
    t = (Ap(1)*AB(1)+Ap(2)*AB(2)+Ap(3)*AB(3))/AB_squared;
    if (t < 0.0) 
      % "Before" A on the line, just return A
      q = A;
    else if (t > 1.0) 
      % "After" B on the line, just return B
      q = B;
    else
      % projection lines "inbetween" A and B on the line
      q = A + t * AB;
    end
    end
end;