run_trial <- function(weights0, Wmult, sigma_dc, sigma_ac, DC, Vdecay, ac_const, beta, etdecay, lrate, discf, noise, platform_x, platform_y, starting_x, starting_y, speed, wall_pun)
{
  # FIXED PARAMETERS OF THE EXPERIMENT
  
  pool_diameter <- 1.4 #Maze diameter in metres (m)
  platform_radius <- 0.06 #Platform radius
  
  N_dc <- 10 #Population of distance cells
  N_ac <- 36 #Population of action cells
  which <- 0
  
  dist <- 0
  wall_zone <- 0
  quadrants <- c(0,0,0,0) #Percentage spent on each quadrant
  
  weights <- weights0 #Initialize modifiable weights
  
  el_tr <- matrix(rep(0, N_dc*N_ac), nrow = N_dc) #Initialize eligibility traces matrix
  
  #Initialize trajectories
  track_x <- starting_x #Current position of trajectory is equal to the starting location of the animal
  track_y <- starting_y
  vel_x <- 0
  vel_y <- 0
  
  # NAVIGATION LOOP
  while ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 > platform_radius^2)
  {
    weights <- weights*(1-noise) + matrix(runif(N_dc*N_ac), nrow = N_dc)*Wmult*noise
    
    #current distance to the wall
    dist.to.wall = pool_diameter/2-sqrt(track_x[length(track_x)]^2+track_y[length(track_y)]^2)
    
    #Calculate DC activation
    DC_activation <- rep(0, N_dc)
    for (i in 1:N_dc){
      DC_activation[i] <- exp(-(dist.to.wall -DC[i])^2/(2*sigma_dc^2))
    }
    
    #Calculate AC activation (i.e. value of the action, Q)
    if (length(track_x) > 1){
      prevQ <- AC_activation[which] #Displays the Q value before movement
    }
    
    AC_activation <- DC_activation %*% weights
    
    #Make an action
    ACsel <- AC_activation^beta
    ACsel <- ACsel / sum(ACsel)
    ASrand <- runif(1)
    which <- 1
    ASsum <- ACsel[1]
    while (which < N_ac && ASsum < ASrand){
      which <- which + 1
      ASsum <- ASsum + ACsel[which]
    }
    
    #Eligibility traces
    el_tr <- el_tr * etdecay
    
    for (j in 1:N_ac){
      itmp <- min(abs(j-which), N_ac-abs(j-which))
      actgaus <- exp(-(itmp*itmp)/(2*sigma_ac*sigma_ac))
      el_tr[,j] <- el_tr[,j] + actgaus*AC_activation[j]*t(t(DC_activation))
    }
    
    #moving direction
    # dir_x = 
    if(track_x[length(track_x)]>=0){
      central.angle = atan(track_y[length(track_y)]/(track_x[length(track_x)]))
    } else {
      central.angle = pi+atan(track_y[length(track_y)]/(track_x[length(track_x)]))
    }
    moving.dir = pi+central.angle+which/N_ac*2*pi
    
    
    
    vel_x = c(vel_x, (vel_x[length(vel_x)]+ac_const*cos(moving.dir))*Vdecay)
    vel_y = c(vel_y, (vel_y[length(vel_y)]+ac_const*sin(moving.dir))*Vdecay)
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
  
  latency <- latency * speed_ts / speed # latency in seconds
  return(list(weights, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, latency))
}


N_dc <- 10 #Population of distance cells [100..300]
N_ac <- 36 #Population of action cells [25..50]

plot_trajectories <- 1 #yes - 1, no - 0
plot_cognitive_maps <- 1 #yes - 1, no - 0
pln <- plot_trajectories + plot_cognitive_maps
Nruns <- 1 #how many runs to run if not plotting anything

pool_diameter <- 1.4 #Maze diameter (\phi) in metres (m)
platform_radius <- 0.06 #Platform radius (m)
sigma_dc <- 0.1 #place cell sigma (standard deviation), in meters [0.05..0.2]
sigma_ac <- 2 #action cell sigma (standard deviation), in action cells [1..3]

etdecay <- 0.83 #Eligibility trace decay (lambda) [0.75..0.95] LESS THAN GAMMA!
beta <- 6 #Exploration-exploitation factor (\beta) [0.5..12]
alpha <- 0.005 #Learning rate (\alpha) [0.005..0.02]
gamma <- 0.95 #Discount factor (\gamma) [0.75..0.95]

Vdecay <- 0.82 #velocity decay [0.75..0.95]
ac_const <- 0.02 #acceleration const [0.01..0.03]
Wnoise <- 0.0004 #Weight noise [0.0001..0.0007]
Wmult <- 0.1 #Weight multiplier [0.05..0.15]
hitwall <- 0.9 #punishment for hitting the wall [0..1]
speed <- 0.175 #mouse speed (m/s) [0.1..0.25]

Ntrials <- 4 #number of trials per day
Ndays <- 80 #number of days


#performance measures to compute: latency, distance, time in target quadrant, opposite quadrant, and wall zone
if (pln > 0.5) #if any plots
{ PMs <- array(rep(0,5*Ndays*Ntrials), c(5,Ndays,Ntrials))
} else {
  PMs <- array(rep(0,5*Ndays*Ntrials*Nruns), c(5,Ndays,Ntrials,Nruns)) }
#multiple runs


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
weights <- matrix(runif(N_dc*N_ac), nrow = N_dc)*Wmult

#Generate distance cells
DC <- rep(0,N_dc)

for (i in 1:N_dc) {
  #For each place cell:
  DC[i] <- runif(1)*(pool_diameter/2) #Random positions of distance cells

}

par(mfrow=c(5,4))

for (day in 1:Ndays) {
  idxs = sample(4) #randomly choose 4 starting locations
  for (trial in 1:Ntrials){
    
    idx <- idxs[trial] #take each location
    starting_x <- starting_xs[idx]
    starting_y <- starting_ys[idx]
    
    modresults <- run_trial (weights, Wmult, sigma_dc, sigma_ac, DC, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x, platform_y, starting_x, starting_y, speed, hitwall)
    #run trial
    weights <- modresults[[1]]
    track_x <- modresults[[2]]
    track_y <- modresults[[3]]
    vel_x <- modresults[[4]]
    vel_y <- modresults[[5]]
    
    #        weights <- wres
    
    PMs[1,day,trial] <- modresults[[9]] #latency
    PMs[2,day,trial] <- modresults[[6]] #dist
    PMs[3,day,trial] <- modresults[[8]][4]*100 #target quadrant
    PMs[4,day,trial] <- modresults[[8]][2]*100 #opposite quadrant
    PMs[5,day,trial] <- modresults[[7]]*100 #wall zone
    #record performance measures
    
    if (plot_trajectories & day%%2==0 & trial==4)
    {
      #plot the maze
      plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),type = "l", xlab = paste("day ",day,", trial ",trial), ylab = "trajectory")
      #plot the trajectory
      lines(track_x, track_y, type = "l")
      #plot the platform
      lines(platform_x+platform_radius*cos(th),platform_y+platform_radius*sin(th),type = "l")
    }
    
    if (plot_cognitive_maps & day%%2==0 & trial==4){
      #plot the maze
      plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),type = "l",xlab = paste("day ",day,", trial ",trial), ylab = "cognitive map")
      #plot the cognitive map
      for (x in (-3:3)*(pool_diameter/6)){
        for (y in (-3:3)*(pool_diameter/6)){
          if (x^2 + y^2 <= (pool_diameter/2)^2){
            x2 = x
            y2 = y
            dist.to.wall = pool_diameter/2-sqrt(x2^2+y2^2)
            for (k in 1:N_ac){
              DC_activation <- rep(0,N_dc)
              for (i in 1:N_dc){
                DC_activation[i] <- exp(-(dist.to.wall -DC[i])^2/(2*sigma_dc^2))
              }
              #Calculate AC activation (i.e. value of the action)
              AC_activation <- rep(0,N_ac)
              for (i in 1:N_ac){
                for (j in 1:N_dc){
                  AC_activation[i] <- AC_activation[i] + DC_activation[j]*weights[j,i]
                }
              }
              #direction of AC activation
              if(x>=0){
                central.angle = atan(y/x)
              } else {
                central.angle = pi+atan(y/x)
              }
              
              moving.dir = pi+central.angle+k/N_ac*2*pi
              
              x2 <- c(x2, x + (AC_activation[k]/10)*cos(moving.dir))
              y2 <- c(y2, y + (AC_activation[k]/10)*sin(moving.dir))
              #                 line([x x2],[y y2],'Color',[k/N_ac 0 1-k/N_ac])
            }
            lines(x2,y2,type = "l",col = "blue")
          }
        }
      }
      # plot the platform
      lines(platform_x+platform_radius*cos(th),platform_y+platform_radius*sin(th),type = "l")
    }
  }
}


}else {
  # run multiple times without plotting!
  
  for (reps in 1:Nruns){
    #Generate initial weights for each run
    weights <- matrix(runif(N_dc*N_ac), nrow = N_dc)*Wmult
    
    #Generate distance cells for each run
    DC <- rep(0,N_dc) #1xN_dc matrix containing the each place cell
    for (i in 1:N_dc) {
      #For each place cell:
      DC[i] <- runif(1)*(pool_diameter/2)#Random positions of distance cells
      
    }
    
    for (day in 1:Ndays){
      idxs = sample(4)  #randomly choose 4 starting locations
      for (trial in 1:Ntrials){
        
        idx <- idxs[trial] #take each location
        starting_x <- starting_xs[idx]
        starting_y <- starting_ys[idx]
        
        modresults <- run_trial (weights, Wmult, sigma_dc, sigma_ac, DC, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x, platform_y, starting_x, starting_y, speed, hitwall)
        #run trial
        weights <- modresults[[1]]
        
        PMs[1,day,trial,reps] <- modresults[[9]] #latency
        PMs[2,day,trial,reps] <- modresults[[6]] #dist
        PMs[3,day,trial,reps] <- modresults[[8]][4]*100 #target quadrant
        PMs[4,day,trial,reps] <- modresults[[8]][2]*100 #opposite quadrant
        PMs[5,day,trial,reps] <- modresults[[7]]*100 #wall zone
        #record performance measures
      }
    }
  }
}

library(tidyr)
library(ggplot2)


if(pln > 0.5){
  latency = PMs[1,,]
  dist = PMs[2,,]
  target_quadrant = PMs[3,,]
  opposite_quadrant = PMs[4,,]
  wall_zone = PMs[5,,]
  
  
  day = c(1:80)
  latency = cbind(latency,day)
  latency = rbind(latency[,c(5,1)],latency[,c(5,2)],latency[,c(5,3)],latency[,c(5,4)])
  colnames(latency)=c('day','latency')
  latency = data.frame(latency)
  
  ggplot(latency,aes(day,latency))+geom_point()+ geom_smooth()
  
}else{
  latency = PMs[1,,,1]
  dist = PMs[2,,,1]
  target_quadrant = PMs[3,,,1]
  opposite_quadrant = PMs[4,,,1]
  wall_zone = PMs[5,,,1]
  
  
  day = c(1:80)
  latency = cbind(latency,day)
  latency = rbind(latency[,c(5,1)],latency[,c(5,2)],latency[,c(5,3)],latency[,c(5,4)])
  colnames(latency)=c('day','latency')
  latency = data.frame(latency)
  
  # 
  # dist = cbind(dist,day)
  # latency = rbind(dist[,c(5,1)],dist[,c(5,2)],dist[,c(5,3)],dist[,c(5,4)])
  # colnames(dist)=c('day','distance')
  # dist = data.frame(dist)
  # 
  # 
  # target_quadrant = cbind(target_quadrant,day)
  # target_quadrant = rbind(target_quadrant[,c(5,1)],target_quadrant[,c(5,2)],target_quadrant[,c(5,3)],target_quadrant[,c(5,4)])
  # colnames(target_quadrant)=c('day','target_quadrant')
  # target_quadrant = data.frame(target_quadrant)
  # 
  # 
  # opposite_quadrant = cbind(opposite_quadrant,day)
  # opposite_quadrant = rbind(opposite_quadrant[,c(5,1)],opposite_quadrant[,c(5,2)],opposite_quadrant[,c(5,3)],opposite_quadrant[,c(5,4)])
  # colnames(opposite_quadrant)=c('day','opposite_quadrant')
  # opposite_quadrant = data.frame(opposite_quadrant)
  # 
  # 
  # wall_zone = cbind(wall_zone,day)
  # wall_zone = rbind(wall_zone[,c(5,1)],wall_zone[,c(5,2)],wall_zone[,c(5,3)],wall_zone[,c(5,4)])
  # colnames(wall_zone)=c('day','wall_zone')
  # wall_zone = data.frame(wall_zone)
  
  
  
  
  ggplot(latency,aes(day,latency))+geom_point()+geom_smooth()
  # ggplot(dist,aes(day,distance))+geom_point()+geom_smooth()
  # ggplot(target_quadrant,aes(day,target_quadrant))+geom_point()+geom_smooth()
  # ggplot(opposite_quadrant,aes(day,opposite_quadrant))+geom_point()+geom_smooth()
  # ggplot(wall_zone,aes(day,wall_zone))+geom_point()+geom_smooth()
}






# par(mfrow=c(1,1))
# plot(latency,latency~day)
