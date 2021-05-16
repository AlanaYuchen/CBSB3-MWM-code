# ======================= run trial function =======================
run_trial <- function(weights0, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, lrate, discf, noise, platform_x, platform_y, starting_x, starting_y, speed, wall_pun)
{
  # FIXED PARAMETERS OF THE EXPERIMENT
  
  pool_diameter <- 1.4 #Maze diameter in metres (m)
  platform_radius <- 0.06 #Platform radius
  
  N_pc <- 211 #Population of place cells
  N_ac <- 36 #Population of action cells
  which <- 0
  
  dist <- 0
  wall_zone <- 0
  quadrants <- c(0,0,0,0) #Percentage spent on each quadrant
  
  weights <- weights0 #Initialize modifiable weights
  
  el_tr <- matrix(rep(0, N_pc*N_ac), nrow = N_pc) #Initialize eligibility traces matrix
  
  #Initialize trajectories
  track_x <- starting_x #Current position of trajectory is equal to the starting location of the animal
  track_y <- starting_y
  vel_x <- 0
  vel_y <- 0
  
  # NAVIGATION LOOP
  while ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 > platform_radius^2) # while out of platform
  {
    weights <- weights*(1-noise) + matrix(runif(N_pc*N_ac), nrow = N_pc)*Wmult*noise
    
    #Calculate PC activation
    PC_activation <- rep(0, N_pc)
    for (i in 1:N_pc){
      PC_activation[i] <- exp(-((track_x[length(track_x)] - PC_x[i])^2 + (track_y[length(track_y)] - PC_y[i])^2)/(2*sigma_pc^2))
    }
    
    #Calculate AC activation (i.e. value of the action, Q)
    if (length(track_x) > 1){
      prevQ <- AC_activation[which] #Displays the Q value before movement
    }
    #print(paste('PC_actiavtion',PC_activation))
    #print(paste('weight',weights))
    AC_activation <- PC_activation %*% weights
    #print(paste('AC_activation',AC_activation))
    #print(paste('AC_activation1',sum(is.na(AC_activation))>0))
    #print(paste('PC_activation',sum(is.na(PC_activation))>0))
    #print(paste('weight',sum(is.na(weights))>0))
    #print(paste('beta',beta))
    #Make an action
    ACsel <- AC_activation^beta
    #print(paste('ACsel1',sum(is.na(ACsel))>0))
    #print(paste('ACsel',ACsel))
    #print(paste('sumACsel',sum(ACsel)))
    ACsel <- ACsel / sum(ACsel)
    #print(paste('ACsel2',sum(is.na(ACsel))>0))
    ASrand <- runif(1)
    which <- 1
    ASsum <- ACsel[1]
    #print(paste('AC_activation2',sum(is.na(AC_activation))>0))
    #print(paste('which',which))
    #print(paste('ACsel',length(ACsel)))
    #print(paste('assum',ASsum))
    #print((paste('ASrand',ASrand)))
    while (which < N_ac & ASsum < ASrand){
      which <- which + 1
      ASsum <- ASsum + ACsel[which]
    }
    
    #Eligibility traces
    el_tr <- el_tr * etdecay
    
    for (j in 1:N_ac){
      itmp <- min(abs(j-which), N_ac-abs(j-which))
      actgaus <- exp(-(itmp*itmp)/(2*sigma_ac*sigma_ac))
      el_tr[,j] <- el_tr[,j] + actgaus*AC_activation[j]*t(t(PC_activation))
    }
    
    vel_x = c(vel_x, (vel_x[length(vel_x)]+ac_const*cos(which/N_ac*2*pi))*Vdecay)
    vel_y = c(vel_y, (vel_y[length(vel_y)]+ac_const*sin(which/N_ac*2*pi))*Vdecay)
    #velocity per time step (not second)
    track_x = c(track_x, track_x[length(track_x)]+vel_x[length(vel_x)])
    track_y = c(track_y, track_y[length(track_y)]+vel_y[length(vel_y)])
    
    #Check if not out of bounds, reset location & speed if so
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > (pool_diameter/2)^2)
    {
      ratio = (track_x[length(track_x)]^2 + track_y[length(track_y)]^2)/((pool_diameter/2)^2)
      track_x[length(track_x)] = track_x[length(track_x)]/sqrt(ratio)
      track_y[length(track_y)] = track_y[length(track_y)]/sqrt(ratio)
      vel_x[length(vel_x)] = track_x[length(track_x)] - track_x[length(track_x)-1]
      vel_y[length(vel_y)] = track_y[length(track_y)] - track_y[length(track_y)-1]
    }
    
    
    if (length(track_x) > 2)
    { if ((track_x[length(track_x)]  - platform_x)^2 + (track_y[length(track_y)]  - platform_y)^2 < platform_radius^2)
    { rew = 10 } #found platform - reward
      else if (track_x[length(track_x)]^2+track_y[length(track_y)]^2 > (0.99*pool_diameter/2)^2)
      { rew = -wall_pun } #hit wall - punishment
      else
      { rew = 0 } #didn't find - no reward
      
      currQ = AC_activation[which]
      tderr = rew + discf*currQ - prevQ #temporal difference error
      weights = pmax(weights + lrate*tderr*el_tr, 0)
    }
    
    laststep = sqrt((track_x[length(track_x)]-track_x[length(track_x)-1])^2 + (track_y[length(track_y)]-track_y[length(track_y)-1])^2)
    dist = dist + laststep
    
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > 0.8*(pool_diameter/2)^2)
    { wall_zone = wall_zone + 1 }
    else if (track_x[length(track_x)] > 0 && track_y[length(track_y)] > 0)
    { quadrants[1] = quadrants[1] + 1 }
    else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] > 0)
    { quadrants[2] = quadrants[2] + 1 }
    else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] < 0)
    { quadrants[3] = quadrants[3] + 1 }
    else
    { quadrants[4] = quadrants[4] + 1 }
    
    if (length(track_x) > 100) # evaluate latency only after 100+ steps to be accurate
    { speed_ts = mean(sqrt((vel_x[-1]^2+vel_y[-1]^2))) # speed in meters/time step
    latency = (length(track_x)-1) * speed_ts / speed # convert to seconds
    if (latency > 60) # if more than a minute, stop
    { break }
    }
    
  }
  
  latency <- length(track_x)-1 # latency in time steps
  wall_zone <- wall_zone/latency
  quadrants <- quadrants/latency
  speed_ts <- mean(sqrt((vel_x[-1]^2+vel_y[-1]^2))) # speed in meters/time step
  # speed per action step from Hanbing
  speed_ps = (vel_x[-1]^2+vel_y[-1]^2)^0.5
  
  # time step
  time_step = speed_ts/speed
  
  # mean turning angle 
  vel = as.matrix(data.frame(vel_x,vel_y))
  angle=c()
  for (steps in 2:(length(vel_x)-1)){
    A = vel[steps,]
    B = vel[steps+1,]
    angle = c(angle, acos((A%*%B)[1,1])/norm(as.matrix(A)*norm(as.matrix(B))))
  }
  angle = angle * 180 / pi
  mean_angle = mean(angle)
  
  # speed standard deviation
  speed_std = sd((vel_x[-1]^2+vel_y[-1]^2)^0.5)
  speed_std = speed_std / time_step
  
  latency <- latency * speed_ts / speed # latency in seconds
  return(list(weights, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, latency, speed_std, speed_ps,mean_angle,time_step))
}

# =========================== main =============================== 
timestart=Sys.time()

N_pc <- 211 #Population of place cells [100..300]
N_ac <- 36 #Population of action cells [25..50]

plot_trajectories <- 0 #yes - 1, no - 0
plot_cognitive_maps <- 0 #yes - 1, no - 0
pln <- plot_trajectories + plot_cognitive_maps
Nruns <- 1 #how many runs to run if not plotting anything

pool_diameter <- 1.4 #Maze diameter (\phi) in metres (m)
platform_radius <- 0.06 #Platform radius (m)
sigma_pc <- c(0.1,0.2) #place cell sigma (standard deviation), in meters [0.05..0.2]
sigma_ac <- c(1,2) #action cell sigma (standard deviation), in action cells [1..3]

etdecay <- 0.75 #Eligibility trace decay (lambda) [0.75..0.95] LESS THAN GAMMA!
beta <- c(0.5,12) #Exploration-exploitation factor (\beta) [0.5..12]
alpha <- c(0.005,0.02) #Learning rate (\alpha) [0.005..0.02]
gamma <- c(0.75,0.95) #Discount factor (\gamma) [0.75..0.95]

Vdecay <- c(0.8,0.93) #velocity decay [0.75..0.95]
ac_const <- c(0.01,0.03) #acceleration const [0.01..0.03]
Wnoise <- c(0.0001,0.0007) #Weight noise [0.0001..0.0007]
Wmult <- c(0.05,0.15) #Weight multiplier [0.05..0.15]
hitwall <- c(0,5)#punishment for hitting the wall [0..1]
speed <- c(0.1,0.25) #mouse speed (m/s) [0.1..0.25]

Ntrials <- 4 #number of trials per day
Ndays <- 8 #number of days
track_x_sum=vector(mode='list',length = Ntrials*Ndays)
track_y_sum=vector(mode='list',length = Ntrials*Ndays)

Npsets = 1
params<-array(rep(0,12*Npsets),c(12,Npsets))

params[1,] =runif(Npsets)*(sigma_pc[2]-sigma_pc[1])+sigma_pc[1]
params[2,] =runif(Npsets)*(sigma_ac[2]-sigma_ac[1])+sigma_ac[1]
params[3,] =runif(Npsets)*(beta[2]-beta[1])+beta[1]
params[4,] =runif(Npsets)*(alpha[2]-alpha[1])+alpha[1]
params[5,] =runif(Npsets)*(gamma[2]-gamma[1])+gamma[1]
params[6,] =runif(Npsets)*(Vdecay[2]-Vdecay[1])+Vdecay[1]
params[7,] =runif(Npsets)*(ac_const[2]-ac_const[1])+ac_const[1]
params[8,] =runif(Npsets)*(Wnoise[2]-Wnoise[1])+Wnoise[1]
params[9,] =runif(Npsets)*(Wmult[2]-Wmult[1])+Wmult[1]
params[10,] =runif(Npsets)*(hitwall[2]-hitwall[1])+hitwall[1]
params[11,] =runif(Npsets)*(speed[2]-speed[1])+speed[1]
params[12,] =runif(Npsets)*(params[5,]-etdecay)+etdecay

#performance measures to compute: latency, distance, time in target quadrant, opposite quadrant, and wall zone
if (pln > 0.5) #if any plots
{ PMs <- array(rep(0,8*Ndays*Ntrials), c(8,Ndays,Ntrials))
  AMs <- array(rep(0,Ndays*Ntrials),c(Ndays,Ntrials))
} else {
  PMs <- array(rep(0,8*Ndays*Ntrials*Nruns), c(8,Ndays,Ntrials,Nruns,Npsets))
  AMs <- array(rep(0,Ndays*Ntrials*Nruns),c(Ndays,Ntrials,Nruns))}
#multiple runs

for (pset in 1:Npsets){
  #print(paste('No.',pset))
  sigma_pc = params[1,pset]
  sigma_ac = params[2,pset]
  beta = params[3,pset]
  alpha = params[4,pset]
  gamma = params[5,pset]
  Vdecay = params[6,pset]
  ac_const = params[7,pset]
  Wnoise = params[8,pset]
  Wmult = params[9,pset]
  hitwall = params[10,pset]
  speed = params[11,pset]
  etdecay = params[12,pset]
  
  #Platform coordinates:
  platform_x <- cos(-pi/4)*pool_diameter/4 #x coordinate
  platform_y <- sin(-pi/4)*pool_diameter/4 #y coordinate
  
  #Starting locations of the modeled animal (4 different ones):
  strad <- pool_diameter/2*0.85 #15% of maze radius to the wall
  starting_xs <- strad * c(cos(pi/6), cos(pi/3), cos(7*pi/6), cos(4*pi/3)) #x coordinates
  starting_ys <- strad * c(sin(pi/6), sin(pi/3), sin(7*pi/6), sin(4*pi/3)) #y coordinates
  
  th <- (0:100)/50*pi #for plotting circles :)
  
  if (pln > 0.5) {
  
  #Generate initial weights
  weights <- matrix(runif(N_pc*N_ac), nrow = N_pc)*Wmult
  
  #Generate place cells
  PC_x <- rep(0,N_pc) #1xN_pc matrix containing the x = 0 coordinate for each place cell
  PC_y <- rep(0,N_pc) #1xN_pc matrix containing the y = 0 coordinate for each place cell
  for (i in 1:N_pc) {
    #For each place cell
    PC_x[i] <- (runif(1) - 0.5)*pool_diameter#Random positions of place cells
    PC_y[i] <- (runif(1) - 0.5)*pool_diameter
    while ((PC_x[i]^2 + PC_y[i]^2 > (pool_diameter/2)^2)){
      #Checks for out of bounds
      PC_x[i] <- (runif(1) - 0.5)*pool_diameter
      PC_y[i] <- (runif(1) - 0.5)*pool_diameter
    }
    }
  par(mfrow=c(Ndays/4,Ntrials*pln/4))
  
  
  
  for (day in 1:Ndays) {
      idxs = sample(4) #randomly choose 4 starting locations
      for (trial in 1:Ntrials){
        
        # fix or variable platform
          #whichplatform = 1 # fix
          whichplatform=sample(c(1,-1),1) # variable
          
          idx <- idxs[trial] #take each location
          starting_x <- starting_xs[idx]
          starting_y <- starting_ys[idx]
  
          modresults <- run_trial (weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, whichplatform*platform_x, whichplatform*platform_y, starting_x, starting_y, speed, hitwall)
          #run trial
          weights <- modresults[[1]]
          track_x <- modresults[[2]]
          track_y <- modresults[[3]]
          vel_x <- modresults[[4]]
          vel_y <- modresults[[5]]
          track_x_sum[[(day-1)*4+trial]]=track_x
          track_y_sum[[(day-1)*4+trial]]=track_y
  #        weights <- wres
          
          PMs[1,day,trial] <- modresults[[9]] #latency
          PMs[2,day,trial] <- modresults[[6]] #dist
          if (whichplatform==1){
            PMs[3,day,trial] <- modresults[[8]][4]*100 #target quadrant
            PMs[4,day,trial] <- modresults[[8]][2]*100 #opposite quadrant
            }else{
            PMs[3,day,trial] <- modresults[[8]][2]*100 #target quadrant
            PMs[4,day,trial] <- modresults[[8]][4]*100 #opposite quadrant
            }
          PMs[5,day,trial] <- modresults[[7]]*100 #wall zone
          PMs[6,day,trial] <- modresults[[10]]*100 # speed_std
          PMs[7,day,trial] <- modresults[[12]] # mean_angle
          PMs[8,day,trial] <- modresults[[13]] # time_step
          #AMs[day,trial] <- modresults[[11]] # speed_ps
 
          #record performance measures
          
          if (plot_trajectories)
          {
            #plot the maze
            plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),type = "l", xlab = paste("day ",day,", trial ",trial), ylab = "trajectory")
            #plot the trajectory
            lines(track_x, track_y, type = "l")
            #plot the platform
            lines(platform_x*whichplatform+platform_radius*cos(th),whichplatform*platform_y+platform_radius*sin(th),type = "l")
          }
  
          if (plot_cognitive_maps){
            #plot the maze
            plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),type = "l",xlab = paste("day ",day,", trial ",trial), ylab = "cognitive map")
            #plot the cognitive map
            for (x in (-3:3)*(pool_diameter/6)){
              for (y in (-3:3)*(pool_diameter/6)){
                  if (x^2 + y^2 <= (pool_diameter/2)^2){
                      x2 = x
                      y2 = y
                      for (k in 1:N_ac){
                          PC_activation <- rep(0,N_pc)
                          for (i in 1:N_pc){
                              PC_activation[i] <- exp(-((x - PC_x[i])^2 + (y - PC_y[i])^2)/(2*sigma_pc^2))
                          }
                  #Calculate AC activation (i.e. value of the action)
                          AC_activation <- rep(0,N_ac)
                          for (i in 1:N_ac){
                              for (j in 1:N_pc){
                                  AC_activation[i] <- AC_activation[i] + PC_activation[j]*weights[j,i]
                              }
                          }
                           x2 <- c(x2, x + (AC_activation[k]/10)*cos(k/N_ac*2*pi))
                           y2 <- c(y2, y + (AC_activation[k]/10)*sin(k/N_ac*2*pi))
          #                 line([x x2],[y y2],'Color',[k/N_ac 0 1-k/N_ac])
                      }
                      lines(x2,y2,type = "l",col = "blue")
                 }
              }
            }
            # plot the platform
            lines(whichplatform*platform_x+platform_radius*cos(th),whichplatform*platform_y+platform_radius*sin(th),type = "l")
          }
      }
  }
  
  } else {
  # run multiple times without plotting!
  #print('else')
  for (reps in 1:Nruns){
    #print('Nruns',reps)
  #Generate initial weights for each run
  weights <- matrix(runif(N_pc*N_ac), nrow = N_pc)*Wmult
  
  #Generate place cells for each run
  PC_x <- rep(0,N_pc) #1xN_pc matrix containing the x = 0 coordinate for each place cell
  PC_y <- rep(0,N_pc) #1xN_pc matrix containing the y = 0 coordinate for each place cell
  for (i in 1:N_pc) {
    #For each place cell:
    PC_x[i] <- (runif(1) - 0.5)*pool_diameter#Random positions of place cells
    PC_y[i] <- (runif(1) - 0.5)*pool_diameter
    while ((PC_x[i]^2 + PC_y[i]^2 > (pool_diameter/2)^2)){
      #Checks for out of bounds
      PC_x[i] <- (runif(1) - 0.5)*pool_diameter
      PC_y[i] <- (runif(1) - 0.5)*pool_diameter
    }
  }
  
  for (day in 1:Ndays){
    idxs = sample(4)  #randomly choose 4 starting locations
      for (trial in 1:Ntrials){
        print(paste('Day:',day,'Trial:',trial))
      # fix or variable platform
        # whichplatform = 1 # fix
        whichplatform=sample(c(1,-1),1) # variable
        
        idx <- idxs[trial] #take each location
        starting_x <- starting_xs[idx]
        starting_y <- starting_ys[idx]
  
        modresults <- run_trial (weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, whichplatform*platform_x, whichplatform*platform_y, starting_x, starting_y, speed, hitwall)
        #run trial
        weights <- modresults[[1]]
        
        PMs[1,day,trial,reps,pset] <- modresults[[9]] #latency
        PMs[2,day,trial,reps,pset] <- modresults[[6]] #dist
        if (whichplatform==1){
          PMs[3,day,trial,reps,pset] <- modresults[[8]][4]*100 #target quadrant
          PMs[4,day,trial,reps,pset] <- modresults[[8]][2]*100 #opposite quadrant
        }else{
          PMs[3,day,trial,reps,pset] <- modresults[[8]][2]*100 #target quadrant
          PMs[4,day,trial,reps,pset] <- modresults[[8]][4]*100 #opposite quadrant
        }
        PMs[5,day,trial,reps,pset] <- modresults[[7]]*100 #wall zone
        #record performance measures
        PMs[6,day,trial,reps,pset] <- modresults[[10]]*100 # speed_std
        PMs[7,day,trial,reps,pset] <- modresults[[12]] # mean_angle
        PMs[8,day,trial,reps,pset] <- modresults[[13]] # time_step
        #AMs[day,trial,reps] <- modresults[[11]] # speed_ps
       }
  }
  }
  }
}
timeend = Sys.time()
# ================== writing output file ==================
dims=dim(PMs)
outf=matrix(ncol = length(dims)+1)
if (length(dims)==5){
  for (a in 1:dims[1]){
    for (b in 1:dims[2]){
      for (c in 1:dims[3]){
        for (d in 1:dims[4]){
          for (e in 1:dims[5]){
            outf=rbind(outf,c(a,b,c,d,e,PMs[a,b,c,d,e]))
          }
        }
      }
    }
  }
  outf=outf[-1,]
  out = data.frame(outf)
  colnames(out)<-c('PM','Day','Trial','reps','pset')
  
}else if(length(dims)==3){
  for (a in 1:dims[1]){
    for (b in 1:dims[2]){
      for (c in 1:dims[3]){
        outf=rbind(outf,c(a,b,c,PMs[a,b,c]))
      }
    }
  }
  outf=outf[-1,]
  out = data.frame(outf)
  colnames(out)<-c('PM','Day','Trial')
}

# set your working directory
# setwd("")
# write.csv(out,'output.csv')