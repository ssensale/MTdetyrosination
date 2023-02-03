clc;clear;reset(RandStream.getGlobalStream, sum(100*clock));

%%%%%Initial Conditions%%%%%
center_rad=1;               % Diameter of center sphere
num_tubules=1;              % Number of Tubules
tubule_size=2;            % Diameter of Tubules
cell_rad=25;                % Diameter of cell
num_enzymes=1;              % Number of Enzymes
enzyme_size=2;            % Diameter of Enzyme
D3D = 1.067;                    % translational diffusion coefficient in 3D
D2D = 10*0.0079;                    % translational diffusion coefficient in 2D
T = 273+37;                    % temperature in kelvin
k_b = (1.380649 * 10^-23)*1E6*1E6;  % Boltzmann's Constant in microm, kg, s, K
kT = k_b * T;               % thermal energy in joules
zeta3D = kT/D3D;            % drag coefficient in 3D
zetaANGLE = kT/D2D;         % angular drag coefficient in 2D
zetaHEIGHT= kT/D2D;         % height drag coefficient in 2D
dt=3;                       % timestep
rc=sqrt(D3D*dt)*3;                       % verlet cutoff
binding_probability=1;                 % binding probability
unbinding_probability=0;%1-exp(-dt*0.365);  % unbinding probability 
probabilitymodifiestubule=1-exp(-dt*0.74);     %probability it will modify the tubule
rescue=0.1299;
rescuemod=0.1370;
catas=0.026;%0.0606;
catasmod=0;0.0319;
growth=0.192;%0.33;
growthmod=0.32;
shrink=0.54;
shrinkmod=0.35;
discretization=100;%6250;               %discretization of the tubules
plots=0;
still_prob = 1;
modified=cell(num_tubules,1);
growing=discretization*ones(num_tubules,1);
%cylinder_locations = onetubule_setup(center_rad, cell_rad);
cylinder_locations = [0,0,1,0,0,25];
vectors_cylinders = tubule_vector(cylinder_locations,num_tubules);
enzyme_locations=enzyme_setup(num_enzymes, enzyme_size, center_rad, cell_rad);
for i=1:num_tubules
    dx(i)=vectors_cylinders(i,4)/discretization;
end;

%%%%%Initial Conformation%%%%%

accept = 1;
%enzyme_locations=[0 (tubule_size+enzyme_size)*.5 8];
%enzyme_locations=[0 10 0];

if plots
    plot_tubule(center_rad, cell_rad,cylinder_locations,tubule_size,discretization,modified,growing)
end;
%plot_enzymes(num_enzymes,enzyme_locations,enzyme_size)
fprintf('###Initial Conformations Built###\n')

%%%%%Make Movements%%%%%
accepted=0;
counter=0;
check2D=ones(num_enzymes,1);
while not(accepted)
    
    [modified,growing]=signal(modified,growing,num_tubules,rescue,rescuemod,catas,catasmod,growth,growthmod,shrink,shrinkmod,dt,discretization,dx);
    [verletW,verletC] = verlets(num_enzymes,enzyme_locations,enzyme_size,rc,cell_rad,center_rad,num_tubules,cylinder_locations,tubule_size,growing,discretization,vectors_cylinders);
    for i=1:num_enzymes
       if not(check2D(i))
            [enzyme_locations(i,:),check2D(i)]=brownian3d(num_tubules,tubule_size, i, enzyme_locations(i,:),zeta3D,kT,dt,center_rad,enzyme_size,cell_rad,verletW,verletC,cylinder_locations,binding_probability,discretization,vectors_cylinders,growing);
       else
            [enzyme_locations(i,:),check2D(i)]=brownian2d_mobile(vectors_cylinders(check2D(i),:),tubule_size, i, enzyme_locations(i,:),zetaANGLE,zetaHEIGHT,kT,dt,center_rad,enzyme_size,cell_rad,cylinder_locations(check2D(i),:),unbinding_probability,check2D(i),discretization,growing(check2D(i)),dx); 
       end
    end
    fprintf('Step %d\n',counter)    
    counter=counter+1;
    
    if plots
        clf
            plot_cell(center_rad, cell_rad)
            hold on    
            plot_tubule(center_rad, cell_rad,cylinder_locations,tubule_size,discretization,modified,growing);
            hold on
            plot_enzymes(num_enzymes,enzyme_locations,enzyme_size)            
            drawnow()

    pause(0.1)
end;
end