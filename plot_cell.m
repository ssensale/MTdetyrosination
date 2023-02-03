function plot_cell = plot_cell(center_rad, cell_rad)


[X1,Y1,Z1]=sphere(10);
XA1=X1*center_rad;
YA1=Y1*center_rad;
ZA1=Z1*center_rad;
[X2,Y2,Z2]=sphere(20);
XB1=X2*cell_rad;
YB1=Y2*cell_rad;
ZB1=Z2*cell_rad;

%figure;
surface(XA1,YA1,ZA1,'FaceAlpha',0.1)
hold on
surface(XB1,YB1,ZB1,'FaceAlpha',0.1)
hold on
view(3)
end