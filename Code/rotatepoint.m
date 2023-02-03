function point2 = rotatepoint(point,angle,versor,point3)
%rotates following right hand rule
    x=point(1);y=point(2);z=point(3);
    a=point3(1);b=point3(2);c=point3(3);
    u=versor(1);v=versor(2);w=versor(3);
    angle=angle;
    point2(1)=(a*(v*v+w*w)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(angle))+x*cos(angle)+sin(angle)*(-c*v+b*w-w*y+v*z);
    point2(2)=(b*(u*u+w*w)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(angle))+y*cos(angle)+(c*u-a*w+w*x-u*z)*sin(angle);
    point2(3)=(c*(u*u+v*v)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(angle))+z*cos(angle)+(-b*u+a*v-v*x+u*y)*sin(angle);
end