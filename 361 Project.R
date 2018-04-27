# library(ggplot2)
# library(remotes)
# Sys.setenv(PATH = paste("C:/PROGRA~1/ImageMagick-7.0.7-Q16",
#                        Sys.getenv("PATH"), sep = ";"))
# library(animation)
# ani.options(convert = 'C:/PROGRA~1/ImageMagick-7.0.7-Q16/convert.exe', interval=0.1)

# ***** GLOBAL VARIABLES *****#
# *** Adjustable values for simulation ***

library(doParallel)

#Particle swarm optimization formula weights
local_weight = 50
global_weight = 1


num_drones = 500 #number of drones in simulation
speed = 1 #radius of neighborhood that drones are allowed to move in
vision <<- 2 # "" are allowed to scan at each time-step
harvest_rate = 1 #no. of plastic particles one drone can pick up per time-step

dimension = 25 #simulation space is a dimensionxdimension matrix
sparseness = 5 #every 1:sparseness patches in simulation has plastic particles
max_particle_conc = 10 #maximum plastic particles in a patch

# *** For use by program ***
global_best <<- c(0,0) #Location of the patch with largest concentration
g_max_conc <<- 0

field = matrix(rep(0, dimension**2),nrow = dimension, ncol=dimension) #representation of simulation space
directions = matrix(c(1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1),ncol=2, byrow=TRUE) #matrix representation of 8 directions

system_cores = detectCores() #max cores for parellelization

#***** FUNCTIONS *****#
#Graphs the field with the drones in it at a certain time-step
graph = function(drones,i){
  
  #Get all the non-zero values out of the field along with all of their associated coord
  nonzero = which(field != 0, arr.ind = T)
  x = nonzero[,"row"]
  y = nonzero[,"col"]
  values = field[nonzero]
  vis = data.frame(x,y,values)
  
  #Graph the plastic particle locations with larger concentrations being bigger
  plot(x, y, xlim=c(1,dimension), ylim=c(1,dimension), col="#68CBFF", pch = 20, cex = values/max_particle_conc*7, main = paste(i))
  
  #Add text to each point showing its concentration
  if(length(values)!=0) { text(x, y, labels=values) }
  #Plot the drones
  points(drones$coord_x, drones$coord_y, pch=6, col=(1:num_drones))
}

#Returns TRUE if vector <x,y> is within bounds of our simulation
#Coordinate system has <1,1> at the bottom left corner ie. 1st cartesian quadrant only
in_bounds = function(x,y){
  return ((x >= 1 && x <= dimension) && (y >= 1 && y <= dimension))
}

#Selects a new direction for the drone to move randomly
random_move = function(drone){
  repeat{
    dx = sample(-1:1, 1)
    dy = sample(-1:1, 1)
    temp = c(drone$coord_x,drone$coord_y) + c(dx,dy) * speed
    if(c(dx,dy)!=c(0,0) && in_bounds(temp[1],temp[2])){
      drone$coord_x = temp[1]
      drone$coord_y = temp[2]
      return(drone)
    }
  }
}

#Rounds to a 2D vector corresponding to an angle n*PI/4 where n is an integer
round_direction = function(vector){
  deg = atan2(vector[2],vector[1]) #Get the vector in degrees
  octant = round(deg/(pi/4), digits = 0)
  return (directions[(octant%%8)+1,])
}


#Calculate new direction vector for a drone using particle swarm optimization
pso_velocity = function(drone){
  new_v = c(0,0)
  r1 = runif(1, min=0, max=1)
  r2 = runif(1, min=0, max=1)
  local_best = c(drone$local_best_x, drone$local_best_y)
  coord = c(drone$coord_x, drone$coord_y)
  
  delta_local = local_best - coord 
  delta_global = global_best - coord
  
  if(all(global_best == c(0,0))){ r2 = 0 }#don't factor in deltas if the best DNE
  if(all(local_best == c(0,0))){ r1 = 0 }
  
  #This is the core formula by Eberhart and Kennedy
  new_v = local_weight * r1 * delta_local + global_weight * r2 * delta_global
  return (round_direction(new_v))
}

#Moves the drone to new coordinate based on PSO
#NOTE: If speed is too high, it might cause the drone to "overshoot". Need to account for this!
move = function(drone){
  #Annoying, but need to map parts of dataframe to vector to simplify vector math
  local_best = c(drone$local_best_x, drone$local_best_y)
  coord = c(drone$coord_x, drone$coord_y)
  
  if(all(local_best == c(0,0)) && all(global_best == c(0,0))){ #Robot doesn't see anything
    return  (random_move(drone))
  }
  else{ #Use particle swarm optimization with weights to adjust direction vector
    direction = pso_velocity(drone)
  }
  
  repeat{ #Make sure the drone isn't going off the field (very unlikely)
    temp_coord <- coord + direction * speed
    if(in_bounds(temp_coord[1],temp_coord[2])){
      #Again, annoying to reassign to dataframe, but must be done 
      drone$dir_x <- direction[1] 
      drone$dir_y <- direction[2] 
      drone$coord_x <- temp_coord[1] 
      drone$coord_y <- temp_coord[2] 
      return (drone)
    }
    speed = speed - 1
  }
}

#Takes a drone object, scans its surroundings for the highest concentration patch, and
#updates its "local_best" variable accordingly
scan = function(drones){
  
  dimension=dimension
  field = field
  
  #Scan the surroundings of each drone in parallel
  r = foreach(kr=iter(drones,by="row"), .combine = rbind) %dopar%{
    local_best = c(0,0)
    max_conc = 0
    vision = 4
    for(i in -vision:vision){ #loop though "offset x", stay within boundaries of vision 
      for(j in -vision:vision){ #loop through "offset y", ""
        x = kr$coord_x + i #apply the offsets to its current position
        y = kr$coord_y + j
        
        if(((x >= 1 && x <= dimension) && (y >= 1 && y <= dimension))){ #check if current (x,y) is in simulation boundary
          if(field[x,y] > max_conc){ #Update local_best and its associated concentration 
            # if the current patch is better than what we have seen so far
            max_conc <- field[x,y]
            local_best = c(x,y)
          }
        }
      }
    }
    
    #Update drone local best using temp variable
    kr$local_best_x <- local_best[1] 
    kr$local_best_y <- local_best[2] 
    kr$concentration <- max_conc
    return(kr)
  }
  #Assign global best values
  max = r[which.max(df_drones$concentration),]
  g_max_conc <<- max$concentration
  global_best <<- c(max$local_best_x,max$local_best_y)
  return (r)
}

#Initialize the field with random concentrations
initialize = function(sparseness, max_particle_conc){
  for(i in 1:(dimension**2/sparseness)){
    go = TRUE
    while(go){ #Need this while loop to make sure we don't "repeat fill" a patch
      x = sample(1:dimension, 1)
      y = sample(1:dimension, 1)
      if(field[x,y] == 0){
        field[x,y] <- sample(1:max_particle_conc, 1)
        go = FALSE
      }
    }
  }
  return (data.frame(field))
}

#Creates a dataframe of drones with randomly initialized coordinates
init_drones = function(){
  coord_x = vector("numeric", length = num_drones)
  coord_y = vector("numeric", length = num_drones)
  dir_x = vector("numeric", length = num_drones)
  dir_y = vector("numeric", length = num_drones)
  local_best_x = vector("numeric", length = num_drones)
  local_best_y = vector("numeric", length = num_drones)
  concentration = vector("numeric", length = num_drones)
  
  for(i in 1:num_drones){
    coord_x[i] = sample(1:dimension, 1)
    coord_y[i] = sample(1:dimension, 1)
    dir_x[i] = 0
    dir_y[i] = 0
    local_best_x[i] = 0
    local_best_y[i] = 0
  }
  
  return (data.frame(coord_x,coord_y,dir_x,dir_y,local_best_x,local_best_y,concentration))
}

#************************************************************************************#
#***** ACTUAL EXECUTION OF PROGRAM *****#
#Create an array of drones

#Set up parellelization
cl <<- makeCluster(detectCores())
registerDoParallel(cl)

#Used for testing multiple global/local weight ratios
ratios = c(0,0.5,1,2,1000)
mean_iterations = vector("numeric",length=length(ratios))
mean_ideal = vector("numeric",length=length(ratios))

#Repeat simulations per ratio
reps = 5

#Used for saving the results of many simulations
final_data = data.frame(ratios, mean_iterations, mean_ideal)

#Save the mean particles at each timestep for many simulations
max_timesteps = 1000
final_time_data = matrix(0,nrow=max_timesteps,ncol=reps)
colnames(final_time_data) = c(paste("R",ratios))
final_time_data = data.frame(final_time_data)

for (l in 1:length(ratios)){#Iterate through ratios
  
  #Set-up for multiple simulations at this ratio
  mean = 0
  iterations = vector(length=reps)
  ideal_it = vector(length=reps)
  history = data.frame(iterations,ideal_it)
  
  #For tracking the particles at each timestep
  time_data = matrix(0,nrow=max_timesteps,ncol=reps)
  colnames(time_data) = c(paste("sim",1:reps))
  time_data = data.frame(time_data)
  
  local_weight <<-ratios[l]
  
  for(j in 1:reps){
    
    #Initialize the field and drones
    field <<- initialize(sparseness, max_particle_conc)
    df_drones <- init_drones()
    
    #Get the total number of particles in the field in parallel to set up initial values
    total_particles = foreach(b=iter(field,by="row"), .combine = rbind) %dopar%{
      sum(b) 
    }
    
    total_particles = sum(total_particles)
    acceptable <- total_particles*0.05 #The "acceptable" no. of particles left to consider the simulation finished
    ideal_iter <- acceptable/(num_drones*harvest_rate) #Number of iterations if the drones could teleport and see everything
    
    i = 0
    #saveGIF({
    repeat{
      
      #Get the total number of particles in the field in parallel to check for simulation completion
      r = foreach(b=iter(field,by="row"), .combine = rbind) %dopar%{
        sum(b) 
      }
      total = sum(r)
      time_data[i+1,j] = total/total_particles
      if(total < acceptable){ #Everything is harvested
        print(paste("Ideal iterations: ", ideal_iter))
        #Save the simulation results to dataframe
        history[j,]$iterations <- i
        history[j,]$ideal_it <- ideal_iter
        break
      }
      
      print(paste("Iteration",i,"% Left: ",total/total_particles))
      #Step where all drones scan their surroundings
      df_drones <- scan(df_drones)
      
      #Get the patches that need to be harvested in parallel
      to_harvest <- foreach(drone = iter(df_drones,by="row"), .combine=rbind)%dopar%{
        if(field[drone$coord_x, drone$coord_y] > 0){ #Harvest this patch
          return(c(drone$coord_x,drone$coord_y))
        }
      }
      
      #Process the new locations of drones in parallel
      df_drones <- foreach(drone = iter(df_drones,by="row"), .combine=rbind)%dopar%{
        if(field[drone$coord_x, drone$coord_y] == 0){ #Nothing to harvest
          return (move(drone))
        }
        return(drone)
      }
      
      #Harvest all the particles that was found in the previous step
      if(!is.null(to_harvest)){
        for(k in 1:dim(to_harvest)[1]){
          x = to_harvest[k,1]
          y = to_harvest[k,2]
          if(field[x,y]>0){
            field[x,y] <- field[x,y] - 1
          }
        }
      }
      
      i=i+1
      graph(df_drones,i) #Graph the current state of the simulation
      
    }
    #}) 
  }
  #Add the results of the simulation at this ratio to the final dataframe
  final_data[l,]$mean_iterations = mean(history$iterations)
  final_data[l,]$mean_ideal = mean(history$ideal_it)
  final_time_data[l] = apply(time_data,1,mean)
}

#save the simulation to a CSV file
mean(history$iterations)
name = paste("simulation",Sys.time(),".csv")
write.csv(final_time_data, file="simulation_time.csv")
write.csv(final_data, file = name)
plot(final_data$ratios, final_data$mean_iterations/final_data$mean_ideal)

stopCluster(cl) #important to run this line if program stops early! prevents memory leaks