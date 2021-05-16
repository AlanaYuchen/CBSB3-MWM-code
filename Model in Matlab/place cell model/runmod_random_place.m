timeset = [];
N_pc = 211; %Population of place cells [100..300]
N_ac = 36; %Population of action cells [25..50]

plot_trajectories = 0; %yes - 1, no - 0
plot_cognitive_maps = 0; %yes - 1, no - 0
pln = plot_trajectories + plot_cognitive_maps;

pool_diameter = 1.4; %Maze diameter (\phi) in metres (m)
platform_radius = 0.06; %Platform radius (m)

% parameters
sigma_pcB = [0.1 0.2]; 
sigma_acB = [1 2]; 

betaB = [0.5 12];
alphaB = [0.005 0.02];
gamaB =[0.75 0.95];

VdecayB = [0.8 0.93]; %Velocity decay
ac_constB = [0.01 0.03]; %acceleration const
WnoiseB = [0.0001 0.0007]; %Weight noise
WmultB = [0.05 0.15];  %Weight multiplier
hitwallB = [0 5]; %punishment for hitting the wall
speedB = [0.1 0.25];
etdecayB = 0.75; % LESS THAN GAMA£¡

% Wmult = zeros(8,4);
% 
% sigma_pc = Wmult;
% sigma_ac = Wmult;
% 
% etdecay = Wmult;
% beta = Wmult;
% alpha = Wmult;
% gamma = Wmult;
% 
% Vdecay = Wmult;
% ac_const = Wmult;
% Wnoise = Wmult;
% 
% hitwall = Wmult;
% speed = Wmult;

Npsets = 1000;
params = zeros(12,Npsets); % 12 parameters

params(1,:) = rand(1,Npsets)*(sigma_pcB(2)-sigma_pcB(1))+sigma_pcB(1);
params(2,:) = rand(1,Npsets)*(sigma_acB(2)-sigma_acB(1))+sigma_acB(1);
params(3,:) = rand(1,Npsets)*(betaB(2)-betaB(1))+betaB(1);
params(4,:) = rand(1,Npsets)*(alphaB(2)-alphaB(1))+alphaB(1);
params(5,:) = rand(1,Npsets)*(gamaB(2)-gamaB(1))+gamaB(1);
params(6,:) = rand(1,Npsets)*(VdecayB(2)-VdecayB(1))+VdecayB(1);
params(7,:) = rand(1,Npsets)*(ac_constB(2)-ac_constB(1))+ac_constB(1);
params(8,:) = rand(1,Npsets)*(WnoiseB(2)-WnoiseB(1))+WnoiseB(1);
params(9,:) = rand(1,Npsets)*(WmultB(2)-WmultB(1))+WmultB(1);
params(10,:) = rand(1,Npsets)*(hitwallB(2)-hitwallB(1))+hitwallB(1);
params(11,:) = rand(1,Npsets)*(speedB(2)-speedB(1))+speedB(1);
params(12,:) = rand(1,Npsets).*(params(5,:,:,:)-etdecayB)+etdecayB;

Nruns = 10; %how many runs to run if not plotting anything
Ntrials = 4;  %number of trials per day
Ndays = 8;  %number of days



if (pln > 0.5) %if any plots
    clf
    PMs = zeros(8,Ndays,Ntrials);  %performance measures to compute: latency, distance
    %time in target quadrant, opposite quadrant, and wall zone
else
    PMs = zeros(8,Ndays,Ntrials,Nruns,Npsets);  %multiple runs
end

% speed per step from Hanbing
if (pln > 0.5) %if any plots
    AMs = cell(Ndays,Ntrials);  %performance measures to compute: speed per step
else
    AMs = cell(Ndays,Ntrials,Nruns);  %multiple runs and multiple steps
end

 tic
 
for pset = 1:Npsets
    
    pset
    toc;
    t1=toc
sigma_pc = params(1,pset);
sigma_ac = params(2,pset);

etdecay = params(12,pset);
beta = params(3,pset);
alpha = params(4,pset);
gamma = params(5,pset);

Vdecay = params(6,pset);
ac_const = params(7,pset);
Wnoise = params(8,pset);

Wmult = params(9,pset);
hitwall = params(10,pset);
speed = params(11,pset);
    

% params(1,:,:,:) = rand(1,Npsets,Ndays,Ntrials)*(sigma_pcB(2)-sigma_pcB(1))+sigma_pcB(1);
% params(2,:,:,:) = rand(1,Npsets,Ndays,Ntrials)*(sigma_acB(2)-sigma_acB(1))+sigma_acB(1);
% params(3,:,:,:) = rand(1,Npsets,Ndays,Ntrials)*(betaB(2)-betaB(1))+betaB(1);
% params(4,:,:,:) = rand(1,Npsets,Ndays,Ntrials)*(alphaB(2)-alphaB(1))+alphaB(1);
% params(5,:,:,:) = rand(1,Npsets,Ndays,Ntrials)*(gamaB(2)-gamaB(1))+gamaB(1);
% params(6,:,:,:) = rand(1,Npsets,Ndays,Ntrials)*(VdecayB(2)-VdecayB(1))+VdecayB(1);
% params(7,:,:,:) = rand(1,Npsets,Ndays,Ntrials)*(ac_constB(2)-ac_constB(1))+ac_constB(1);
% params(8,:,:,:) = rand(1,Npsets,Ndays,Ntrials)*(WnoiseB(2)-WnoiseB(1))+WnoiseB(1);
% params(9,:,:,:) = rand(1,Npsets,Ndays,Ntrials)*(WmultB(2)-WmultB(1))+WmultB(1);
% params(10,:,:,:) = rand(1,Npsets,Ndays,Ntrials)*(hitwallB(2)-hitwallB(1))+hitwallB(1);
% params(11,:,:,:) = rand(1,Npsets,Ndays,Ntrials)*(speedB(2)-speedB(1))+speedB(1);
% params(12,:,:,:) = rand(1,Npsets,Ndays,Ntrials).*(params(5,:,:,:)-etdecayB)+etdecayB;


% then for the real model

%Platform coordinates:
platform_x = cos(-pi/4)*pool_diameter/4; %x coordinate
platform_y = sin(-pi/4)*pool_diameter/4; %y coordinate

%Starting locations of the modeled animal (4 different ones):
strad = pool_diameter/2*0.9; %10% of maze radius to the wall
starting_xs = strad * [cos(pi/6) cos(pi/3) cos(7*pi/6) cos(4*pi/3)]; %x coordinates
starting_ys = strad * [sin(pi/6) sin(pi/3) sin(7*pi/6) sin(4*pi/3)]; %y coordinates

th = 0:pi/50:2*pi; %for plotting circles :)

if (pln > 0.5)

%Generate initial weights
weights = rand(N_pc,N_ac)*Wmult;

%Generate place cells
PC_x = zeros(1,N_pc); %1xN_pc matrix containing the x = 0 coordinate for each place cell
PC_y = zeros(1,N_pc); %1xN_pc matrix containing the y = 0 coordinate for each place cell
for i = 1:N_pc %For each place cell:
    PC_x(i) = (rand - 0.5)*pool_diameter; %Random positions of place cells
    PC_y(i) = (rand - 0.5)*pool_diameter;
    while (PC_x(i)^2 + PC_y(i)^2 > (pool_diameter/2)^2) %Checks for out of bounds
        PC_x(i) = (rand - 0.5)*pool_diameter;
        PC_y(i) = (rand - 0.5)*pool_diameter;
    end
end
    
for day = 1:Ndays
    for trial = 1:Ntrials
        
        whichplatform = 2*randi(2) - 3; %variable
        %whichplatform = 1; %fixed
        
        idx = randi(4); %randomly choose one of 4 starting locations
        starting_x = starting_xs(idx);
        starting_y = starting_ys(idx);
        

        [wres, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants,...
            latency, speed_std,speed_ps,mean_angle,time_step] = ...
        run_trial (weights, params(9,pset), params(1,pset),params(2,pset), PC_x, PC_y, params(6,pset), params(7,pset), params(3,pset), params(12,pset), ...
    params(4,pset), params(5,pset), params(8,pset), whichplatform*platform_x, whichplatform*platform_y, starting_x, starting_y, params(11,pset), params(10,pset));
        
        
        %run trial
        weights = wres;
               
        PMs(1,day,trial) = latency;
        PMs(2,day,trial) = dist;
        PMs(3,day,trial) = quadrants(4)*100; % target
        PMs(4,day,trial) = quadrants(2)*100; % opposite
        PMs(5,day,trial) = wall_zone*100;
        PMs(6,day,trial) = speed_std*100;
        PMs(7,day,trial) = mean_angle;
        PMs(8,day,trial) = time_step;
        AMs{day,trial}= speed_ps;
        
        %record performance measures
        
        if (plot_trajectories)
        subplot(Ndays, Ntrials*pln,(day-1)*Ntrials*pln+(trial-1)*pln+1);
        hold on
        %plot the trajectory
        for i = 1:(length(track_x))-1
            line(track_x(i:i+1),track_y(i:i+1),'Color',[i/length(track_x),0,1-i/length(track_x)]);
        end
        %plot the maze and platform
        plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),'k');
        plot(whichplatform*platform_x+platform_radius*cos(th),whichplatform*platform_y+platform_radius*sin(th),'k');
        end
        
        if (plot_cognitive_maps)
        subplot(Ndays, Ntrials*pln,(day-1)*Ntrials*pln+trial*pln);
        hold on
        %plot the cognitive map
        for x = -(pool_diameter/2):(pool_diameter/6):(pool_diameter/2)
            for y = -(pool_diameter/2):(pool_diameter/6):(pool_diameter/2)
                if (x^2 + y^2 <= (pool_diameter/2)^2)
                    for k = 1:N_ac
                        PC_activation = zeros(1,N_pc);
                        for i = 1:N_pc
                            PC_activation(i) = exp(-((x - PC_x(i))^2 + (y...
                            - PC_y(i))^2)/(2*sigma_pc^2));
                        end
                %Calculate AC activation (i.e. value of the action)
                        AC_activation = zeros(1,N_ac);
                        for i = 1:N_ac
                            for j = 1:N_pc
                                AC_activation(i) = AC_activation(i) + ...
                                    PC_activation(j)*weights(j,i);
                            end
                        end
                        x2 = x + (AC_activation(k)/10)*cos(k/N_ac*2*pi);
                        y2 = y + (AC_activation(k)/10)*sin(k/N_ac*2*pi);
                        hold on;
                        line([x x2],[y y2],'Color',[k/N_ac 0 1-k/N_ac]);
                    end
                end
            end
        end
        %plot the maze and platform
        plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),'k');
        plot(whichplatform*platform_x+platform_radius*cos(th),whichplatform*platform_y+platform_radius*sin(th),'k'); 
        end

    end
end
else
% run multiple times without plotting!=====================================

for rep = 1:Nruns
%Generate initial weights for each run
weights = rand(N_pc,N_ac)*Wmult;

%Generate place cells for each run
PC_x = zeros(1,N_pc); %1xN_pc matrix containing the x = 0 coordinate for each place cell
PC_y = zeros(1,N_pc); %1xN_pc matrix containing the y = 0 coordinate for each place cell
for i = 1:N_pc %For each place cell:
    PC_x(i) = (rand - 0.5)*pool_diameter; %Random positions of place cells
    PC_y(i) = (rand - 0.5)*pool_diameter;
    while (PC_x(i)^2 + PC_y(i)^2 > (pool_diameter/2)^2) %Checks for out of bounds
        PC_x(i) = (rand - 0.5)*pool_diameter;
        PC_y(i) = (rand - 0.5)*pool_diameter;
    end
end

for day = 1:Ndays
    for trial = 1:Ntrials
        
        idx = randi(4); %randomly choose one of 4 starting locations
        starting_x = starting_xs(idx);
        starting_y = starting_ys(idx);
        
        whichplatform = 2*randi(2) - 3; %variable
        %whichplatform = 1; %fixed
        
    
        [wres, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants,...
            latency, speed_std, speed_ps,mean_angle,time_step] = ...
        run_trial (weights, params(9,pset), params(1,pset),params(2,pset), PC_x, PC_y, params(6,pset), params(7,pset), params(3,pset), params(12,pset), ...
    params(4,pset), params(5,pset), params(8,pset), whichplatform*platform_x, whichplatform*platform_y, starting_x, starting_y, params(11,pset), params(10,pset));

        %run trial
        weights = wres;    
        PMs(1,day,trial,rep,pset) = latency;
        PMs(2,day,trial,rep,pset) = dist;      
        if(whichplatform == 1)       
            PMs(3,day,trial,rep,pset) = quadrants(4)*100; %target 
            PMs(4,day,trial,rep,pset) = quadrants(2)*100; %opposite
        elseif (whichplatform == -1)
            PMs(3,day,trial,rep,pset) = quadrants(2)*100; %target
            PMs(4,day,trial,rep,pset) = quadrants(4)*100; %opposite
        end
        PMs(5,day,trial,rep,pset) = wall_zone*100;
        PMs(6,day,trial,rep,pset) = speed_std*100;
        PMs(7,day,trial,rep,pset) = mean_angle;
        PMs(8,day,trial,rep,pset) = time_step;
        
        %AMs{day,trial,rep} = speed_ps;
       
    end
end
end
end
toc;
t2 = toc;
timeset = [timeset t2-t1];

end



%record performance measures

PMs = real(PMs);
dims = size(PMs);


outf = [];
if (length(dims) == 3)
    for i = 1:dims(1)
        for j = 1:dims(2)
            for k = 1:dims(3)
                outf = [outf; i j k PMs(i,j,k)];
            end
        end
    end
elseif (length(dims) == 5)     
    PMs = squeeze(mean(PMs(:,:,:,:,:),4));
    dims = size(PMs);
    for i = 1:dims(1)
        for j = 1:dims(2)
            for k = 1:dims(3)
                for l = 1:dims(4)
                    outf = [outf; i j k l PMs(i,j,k,l)];

                end
            end
        end
    end
end

PMmodel_var_100psets = PMs;
%save PMmodel_var_100psets
%save(sprintf('PMs5'), 'PMs')

csvwrite('output_variable_1000trial.csv',outf);
csvwrite('output_variable_param_1000trial.csv',params);