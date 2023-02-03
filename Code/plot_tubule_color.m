function plot_tubule_color = plot_tubule_color(cylinder_location,tubule_size,modified,growing,discretization)

vector=cylinder_location([4,5,6])-cylinder_location([1,2,3]);


%coords(abs(growing)+1,1) coords(abs(growing)+1,2) coords(abs(growing)+1,3)
%cylinder_location(1) + vector(1)*abs(growing)/discretization;

drawCylinder3([cylinder_location(1)+vector(1)*abs(growing)/discretization cylinder_location(2)+vector(2)*abs(growing)/discretization cylinder_location(3)+vector(3)*abs(growing)/discretization cylinder_location(1)+vector(1) cylinder_location(2)+vector(2) cylinder_location(3)+vector(3) tubule_size],12,'open');
hold on
drawCylinder1([cylinder_location(1) cylinder_location(2) cylinder_location(3) cylinder_location(1)+vector(1)*abs(growing)/discretization cylinder_location(2)+vector(2)*abs(growing)/discretization cylinder_location(3)+vector(3)*abs(growing)/discretization tubule_size],12,'open')
hold on
if not(isempty(modified))
    for i=1:size(modified,1)
        drawCylinder2([cylinder_location(1)+vector(1)*(modified(i)-1)/discretization cylinder_location(2)+vector(2)*(modified(i)-1)/discretization cylinder_location(3)+vector(3)*(modified(i)-1)/discretization cylinder_location(1)+vector(1)*(modified(i))/discretization cylinder_location(2)+vector(2)*(modified(i))/discretization cylinder_location(3)+vector(3)*(modified(i))/discretization tubule_size],12,'open')
        hold on
    end;
end;



end
