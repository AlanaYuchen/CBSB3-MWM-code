%Estimate parameters for all trials given mouse data, for all trials
function [best_Wmult, best_sigma_pc, best_sigma_ac, best_Vdecay, ...
    best_ac_const, best_beta, best_etdecay, best_alpha, ...
    best_gamma, best_noise, best_wpun, best_gof, best_PMs, goferr_set, params, PMs] = estpar_iter3 (mousedata, allvars,strain,temp)

% mousedata - 6 measures x 8 days x 4 trials
N_pc = 211; %Population of place cells
N_ac = 36;
pool_diameter = 1.4; %Maze diameter (\phi) in metres (m)
platform_radius = 0.06; %Platform radius
etdecaymult = 0.98;

sigma_pcB = [0.1 0.2]; % old [0.1 0.2];
sigma_acB = [1.05 1.05]; % old[1 2];
VdecayB = [0.8 0.93]; %Velocity decay
ac_constB = [0.01 0.03]; %acceleration const
betaB = [0.5 12]; %Exploration-exploitation factor (\beta)
alphaB = [0.005 0.02]; %Learning rate (\alpha)
gammaB = [0.75 0.95]; %Discounting factor (\gamma)
etdecayB = [0.75 0.75];%Eligibility trace decay (lambda) old [0.75..0.95] LESS THAN GAMMA!
WnoiseB = [0.00055 0.00055]; %Weight noise
WmultB = [0.095 0.095];  %Weight multiplier
hitwallB = [3.75 3.75]; %punishment for hitting the wall

best_Wmult = zeros(8,4);
best_sigma_pc = best_Wmult;
best_sigma_ac = best_Wmult;
best_Vdecay = best_Wmult;
best_ac_const = best_Wmult;
best_beta = best_Wmult;
best_etdecay = best_Wmult;
best_alpha = best_Wmult;
best_gamma = best_Wmult;
best_noise = best_Wmult;
best_wpun = best_Wmult;
best_gof = best_Wmult;

best_PMs = zeros(8,4,8);

Nsets = 5000;
Nruns = 5;


goferr_set = zeros(Nsets,8,4);

params = zeros(11,Nsets,8,4);
PM = zeros(8,Nruns,Nsets,8,4);
PMs=zeros(8,Nsets,8,4);


%The ranges of parameters
rangesize = zeros(11,1);

rangesize(1) = sigma_pcB(2)-sigma_pcB(1);
rangesize(2) = sigma_acB(2)-sigma_acB(1);
rangesize(3) = VdecayB(2)-VdecayB(1);
rangesize(4) = ac_constB(2)-ac_constB(1);
rangesize(5) = betaB(2)-betaB(1);
rangesize(6) = alphaB(2)-alphaB(1);
rangesize(7) = gammaB(2)-gammaB(1);

params(1,:,:,:) = rand(1,Nsets,8,4)*(sigma_pcB(2)-sigma_pcB(1))+sigma_pcB(1); %1*8*4*Nsets
params(2,:,:,:) = rand(1,Nsets,8,4)*(sigma_acB(2)-sigma_acB(1))+sigma_acB(1);
params(3,:,:,:) = rand(1,Nsets,8,4)*(VdecayB(2)-VdecayB(1))+VdecayB(1); % constant
params(4,:,:,:) = rand(1,Nsets,8,4)*(ac_constB(2)-ac_constB(1))+ac_constB(1);  % constant
params(5,:,:,:) = rand(1,Nsets,8,4)*(betaB(2)-betaB(1))+betaB(1);
params(6,:,:,:) = rand(1,Nsets,8,4)*(alphaB(2)-alphaB(1))+alphaB(1);
params(7,:,:,:) = rand(1,Nsets,8,4)*(gammaB(2)-gammaB(1))+gammaB(1);
params(8,:,:,:) = rand(1,Nsets,8,4)*(etdecayB(2)-etdecayB(1))+etdecayB(1);
params(9,:,:,:) = rand(1,Nsets,8,4)*(WnoiseB(2)-WnoiseB(1))+WnoiseB(1);
params(10,:,:,:) = rand(1,Nsets,8,4)*(WmultB(2)-WmultB(1))+WmultB(1);
params(11,:,:,:) = rand(1,Nsets,8,4)*(hitwallB(2)-hitwallB(1))+hitwallB(1); 

%Platform coordinates:
platform_x = cos(-pi/4)*pool_diameter/4; %x coordinate
platform_y = sin(-pi/4)*pool_diameter/4; %y coordinate
        
%Starting location of the modeled animal (4 DIFFERENTS):
strad = pool_diameter/2*0.9; % 90% of distance to the wall
starting_xs = strad * [cos(pi/6) cos(pi/3) cos(7*pi/6) cos(4*pi/3)]; %x coordinates
starting_ys = strad * [sin(pi/6) sin(pi/3) sin(7*pi/6) sin(4*pi/3)]; %y coordinates

weights = rand(Nruns,N_pc,N_ac);
 

PC_x = zeros(Nruns,N_pc); %1xN_pc matrix containing the x = 0 coordinate for
                      %each place cell
PC_y = zeros(Nruns,N_pc); %1xN_pc matrix containing the y = 0 coordinate for
                      %each place cell

for repeat = 1:Nruns % place cell position differs by run (replication)
    for i = 1:N_pc %For each place cell:
        PC_x(repeat,i) = (rand - 0.5)*pool_diameter; %Random positions of place cells
        PC_y(repeat,i) = (rand - 0.5)*pool_diameter;
        while (PC_x(repeat,i)^2 + PC_y(repeat,i)^2 > (pool_diameter/2)^2) %Checks for out of bounds
            PC_x(repeat,i) = (rand - 0.5)*pool_diameter;
            PC_y(repeat,i) = (rand - 0.5)*pool_diameter;
        end
    end
end

close all
no_iter = 3;
weightsall = zeros(no_iter,Nruns,N_pc,N_ac);
for s = 1:no_iter
    weightsall(s,:,:,:) = weights;
end
    
prop_range = [1 0.2 0.04]; %ranges of parameters for each iteration
color_histogram = ['r', 'g', 'b'];

for day = 1:8
    for trial = 1:4
        day
        trial
        figure;
        bestgof = 1e10;
        Nbest = 100; 
        Nreps = Nsets / Nbest; %10 here
        gofminsofar =ones(Nsets*no_iter,1)*bestgof;
        
        %Iterations
        for iter = 1:no_iter
            %Setting parameter sets for this iteration
            if (iter == 1)
                params2 = squeeze(params(:,:,day,trial)); %params2 = (11xNsets)
            else
                bests = params2(:,indexes); %bests = (11xNbest) %indexes come from the end of iteration
                to_add = (rand(11,Nsets)*2-1).*rangesize * prop_range(iter);%to_add = (11xNsets) ???appropriate noiose
                bests = repmat(bests, 1, Nreps); %Nreps psets for each of the original best Nbest psets (10 for 30)     
                params2 = to_add + bests; %params2 = (11xNsets)                             
            end
            
            %Finding gofs for this iteration
            for pset = 1:Nsets
                pset
                iter

                sigma_pc = params2(1,pset);
                sigma_ac = params2(2,pset);
                Vdecay = params2(3,pset);
                ac_const = params2(4,pset);
                beta = params2(5,pset);
                lrate = params2(6,pset);
                discf = params2(7,pset);
                etdecay = params2(8,pset);
                noise = params2(9,pset);
                Wmult = params2(10,pset);
                wall_pun = params2(11,pset);

                speed = mousedata(8,day,trial)*0.01;

                % initialize weight0, not really sure what weights
                if (day == 1 && trial == 1)
                    %wcurr = weights*params2(10,pset);
                    wcurr = squeeze(weightsall(iter,:,:,:))*params2(10,pset);
                else
                    %wcurr = weights;
                    wcurr = squeeze(weightsall(iter,:,:,:));
                end
                wres = wcurr;
                resPMs = zeros(8,Nruns);

                for repeat = 1:Nruns
                    stidx = randi(4);
                    starting_x = starting_xs(stidx);
                    starting_y = starting_ys(stidx);


                    time1 = tic;
                    [wres1, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants,...
                    latency, speed_std, mean_angle,time_step] = ...
                    run_trial (squeeze(wcurr(repeat,:,:)), Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, ...
                    lrate, discf, noise, platform_x, platform_y, starting_x, starting_y, speed, wall_pun);

                    wres(repeat,:,:) = reshape(wres1,1,N_pc,N_ac);
                    resPMs(1,repeat) = latency; % size: 8*repeat
                    resPMs(2,repeat) = dist*100;
                    resPMs(3,repeat) = quadrants(4);
                    resPMs(4,repeat) = quadrants(2);
                    resPMs(5,repeat) = wall_zone;
                    resPMs(6,repeat) = speed_std*100;
                    resPMs(7,repeat) = mean_angle;
                    resPMs(8,repeat) = time_step;

                    PM(1,repeat,pset,day,trial) = latency;
                    PM(2,repeat,pset,day,trial) = dist*100;
                    PM(3,repeat,pset,day,trial) = quadrants(4);
                    PM(4,repeat,pset,day,trial) = quadrants(2);
                    PM(5,repeat,pset,day,trial) = wall_zone;
                    PM(6,repeat,pset,day,trial) = speed_std*100;
                    PM(7,repeat,pset,day,trial) = mean_angle;
                    PM(8,repeat,pset,day,trial) = time_step;

                end
                PMs(:,pset,day,trial) = median(PM(:,:,pset,day,trial),2);

                medianPMs = median(resPMs,2); % take median on repeats size: 8

                % calculate chi based on the median over repeats
                goferr = sum((medianPMs(1:7,:)-squeeze(mousedata(1:7,day,trial))).^2./squeeze(allvars(1:7,day,trial)).^2);
                goferr_set(pset,day,trial) = goferr;
                if (goferr < bestgof)
                    bestgof = goferr;
                    bestw = wres;
                    bestpset = pset;
                    bestPMs = medianPMs;
                    
                    %weights = bestw; %This was here accidentally so it
                    %made the iteration progression better than expected I guess
                    
                     best_sigma_pc(day,trial) = params2(1,bestpset);
                     best_sigma_ac(day,trial) = params2(2,bestpset);
                     best_Vdecay(day,trial) = params2(3,bestpset);
                     best_ac_const(day,trial) = params2(4,bestpset);
                     best_beta(day,trial) = params2(5,bestpset);
                     best_alpha(day,trial) = params2(6,bestpset);
                     best_gamma(day,trial) = params2(7,bestpset);
                     best_etdecay(day,trial) = params2(8,bestpset);
                     best_noise(day,trial) = params2(9,bestpset);
                     best_Wmult(day,trial) = params2(10,bestpset);
                     best_wpun(day,trial) = params2(11,bestpset);
                    best_gof(day,trial) = bestgof;
                    best_PMs(day,trial,:) = reshape(bestPMs,1,1,8);
                end
            end 
   
            %Best Nbest psets from the Nsets in this iteration
            [gofs, indexes] = mink(goferr_set(:,day,trial), Nbest);
            
            %Best Nbest psets from all sets tried so far
            %gofminsofar = [gofminsofar; goferr_set(:,day,trial)]
            gofminsofar(((iter-1)*Nsets + 1):iter*Nsets) = goferr_set(:,day,trial); %Nsets*iter many gofs
            [gofs1,indexes_best] = mink(gofminsofar, Nbest);
            
            histogram(real(gofs1),'NumBins', 10,'FaceColor', color_histogram(iter))
            title(['Mouse S',num2str(strain),'T',num2str(temp),'Day ',num2str(day),', trial ',num2str(trial)]);
            hold on
            
            weightsall(iter,:,:,:) = bestw;
        end
        
        legend('Iter 1',' Iter 2',' Iter 3')
        
        disp(['Day ',num2str(day),', trial ',num2str(trial),' - ',num2str(bestgof)])
        %Best PMs and parameter set for this trial
        params2(:,bestpset)' %(11x1)
        bestPMs' %(8x1)
    end
end
