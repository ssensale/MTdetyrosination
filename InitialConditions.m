%%%%%Initial Conditions%%%%%
dt=0.01;                       % timestep
stepsverlet1=25/dt;
totalsteps=5000/dt;
cell_rad=0.5*25;
center_rad=0.5; 
discretization=(cell_rad-center_rad)*1000/20;
kmod=0.74;
kunb=0.33;
kunbslide=0.63694;
kbin=1000;
% Radius of center sphere
D3D=1.067;
D2D=0.0079;
discretization=(cell_rad-center_rad)*1000/0.6;
D3D=1.067;
D1D=0.0079;
still_prob=0.5;
tubule_size=0.5*25/1000;            % Radius of Tubules

highstability=1;
nodynamics=0;
num_tubules=50;
num_enzymes=100;

kbin=1000;
enzyme_size=0.5*7/1000; 

if highstability==1
	rescue=0.175;
	catas=0.026;
	growth=0.192;
	growth0=0.192;
	shrink=0.218;
    concentrationoftubulin=35;%micromolar
else
	rescue=0.086;
	catas=0.053;
	growth=0.142;
	growth0=0.142;
	shrink=0.188;
    concentrationoftubulin=25;%micromolar
end
rescuemod=rescue/0.95;
growthmod=growth/0.93;
growthmod0=growth0/0.93;

catasmod=catas/(1.89);
shrinkmod=shrink/(1.54);
CUTOFF=(discretization/50)/5;







T = 273+37;                    % temperature in kelvin
k_b = (1.380649 * 10^-23)*1E6*1E6;  % Boltzmann's Constant in microm, kg, s, K
kT = k_b * T;               % thermal energy in joules
zeta3D = kT/D3D;            % drag coefficient in 3D
zetaANGLE = kT/0.0079;         % angular drag coefficient in 2D
zetaHEIGHT= kT/D2D;         % height drag coefficient in 2D
rc=sqrt(D3D*dt)*3;                       % verlet cutoff
binding_probability=1-exp(-dt*kbin);                 % binding probability
unbinding_probability=1-exp(-dt*kunb);  % unbinding probability 
unbinding_probabilityslide=1-exp(-dt*kunbslide);  % unbinding probability 
probabilitymodifiestubule=1-exp(-dt*kmod);     %probability it will modify the tubule


modified=cell(num_tubules,1);
%%%%%Initial Conformation%%%%%
cylinder_locations = tubule_setup(num_tubules, tubule_size, center_rad, cell_rad, cell_rad, 1);
vectors_cylinders=tubule_vector(cylinder_locations,num_tubules);
enzyme_locations=enzyme_setup(num_enzymes, enzyme_size, center_rad, cell_rad);
fprintf('###Initial Conformations Built###\n')
for i=1:num_tubules
    dx(i)=vectors_cylinders(i,4)/discretization;
end

if nodynamics
    growing=(discretization)*ones(num_tubules,1);%state of the tubule, +1 grows from 1, -1 shrinks from 1 - also saves location of end
    largos=(cell_rad-center_rad)*ones(num_tubules,1);
else
    growing=ones(num_tubules,1);%state of the tubule, +1 grows from 1, -1 shrinks from 1 - also saves location of end
    largos=zeros(num_tubules,1);
end



counterenzymes=zeros(num_enzymes,1);

accepted=0;
counter=0;
check2D=zeros(num_enzymes,1);
nucleation=0.0005;
pnuc=1-exp(-nucleation*dt);

name='Cylinder.dat';
id=(1:num_tubules)';
cylinder_loc=cat(2,id,cylinder_locations);
writematrix(cylinder_loc,name,'Delimiter','tab');

concentrationdependent=1;
if concentrationdependent
    kon=1/60;
    Vcell=(4/3)*3.1416*(cell_rad^3-center_rad^3);
    TubulinTotalNumb=concentrationoftubulin*602*Vcell;
end

counter2=0;
frequency2=stepsverlet1;


countofbindsF=0;
countofbindsM=0;

