
# library(ggplot2)
# library(remotes)
# Sys.setenv(PATH = paste("C:/PROGRA~1/ImageMagick-7.0.7-Q16",
#                        Sys.getenv("PATH"), sep = ";"))
# library(animation)
# ani.options(convert = 'C:/PROGRA~1/ImageMagick-7.0.7-Q16/convert.exe', interval=0.1)

# ***** GLOBAL VARIABLES *****#
#Adjustable values for simulation
reps = 10

local_weight = 10
global_weight = 1

num_drones = 50
speed = 1
vision = 4
harvest_rate = 1

dimension = 20
sparseness = 5
max_particle_conc = 10


#For use by program
global_best <<- c(0,0)
g_max_conc <<- 0

field = matrix(rep(0, dimension**2),nrow = dimension, ncol=dimension)
saved_field = matrix(rep(0, dimension**2),nrow = dimension, ncol=dimension)

directions = matrix(c(1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1),ncol=2, byrow=TRUE)

#***** FUNCTIONS *****#
#Graphs the field with the drones in it at a certain time-step
graph = function(drones,i){
  #Get all the non-zero values out of the field along with all of their associated coord
  nonzero = which(field != 0, arr.ind = T)
  x = nonzero[,"row"]
  y = nonzero[,"col"]
  values = field[nonzero]
  vis = data.frame(x,y,values)
  
  
  #In-progress better visualization using ggplot and gganimate
  # theme_set(theme_bw())
  # p <- ggplot(data = vis, aes(x=x, y=y, colour=values, frame = i)) + geom_point(size = values) +
  #   annotate("point", x=drone1@coord[1], y=drone1@coord[2], shape = 6) +
  #   annotate("text", x=x, y=y, label=values, size = 5) + scale_colour_gradientn(colours=rainbow(4))
  # gganimate(p)
  
  #Graph the plastic particle locations with larger concentrations being bigger
  plot(x, y, xlim=c(1,dimension), ylim=c(1,dimension), col="#68CBFF", pch = 20, cex = values/2, main = paste(i))
  
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

#Harvests the plastic at the current location at the drone's collection speed
#May need to check for negative values - "overshooting"
harvest = function(drone){
  x = drone$coord_x
  y = drone$coord_y
  field[x,y] = field[x,y] - harvest_rate
  if(c(x,y) == global_best && field[x,y] == 0){
    global_best <<- c(0,0)
    g_max_conc <<- 0
  }
  return (field)
}

#Decide whether the drone should move or stay and harvest the patch
#NOTE: If it sees a better patch while it's harvesting, should it move? 
#This doesn't need to be a function unless we're answering that question
# action = function(drone){
#   return (drone)
# }

#Takes a drone object, scans its surroundings for the highest concentration patch, and
#updates its "local_best" variable accordingly
scan = function(drones){
  
  for(k in 1:num_drones){ #scan for each drone
    
    #initialize some variables before inner loop for scanning
    local_best = c(0,0)
    max_conc <<- 0
    
    for(i in -vision:vision){ #loop though "offset x", stay within boundaries of vision 
      for(j in -vision:vision){ #loop through "offset y", ""
        x = drones[k,]$coord_x + i #apply the offsets to its current position
        y = drones[k,]$coord_y + j
        
        if(in_bounds(x,y)){ #check if current (x,y) is in simulation boundary
          saved_field[x,y] <<- field[x,y]
          #dfsaved = update_saved()
        
          if(field[x,y] > max_conc){ #Update local_best and its associated concentration 
            # if the current patch is better than what we have seen so far
            max_conc <<- field[x,y]
            local_best = c(x,y)
          }
        }
      }
    }
    g_max_conc <<- max(saved_field)
    global_best <<- which(saved_field == max(saved_field), arr.ind = T)
    global_best <<- global_best[1,]
    while (max_conc == 0 && g_max_conc == 0){
      return(random_move(drones))
    }
    
    #Assign temp local best to the drone object itself
    drones[k,]$local_best_x <- local_best[1] 
    drones[k,]$local_best_y <- local_best[2] 
  }
  
  return (drones)
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
  return (field)
}

#Creates a dataframe of drones with randomly initialized coordinates
init_drones = function(){
  coord_x = vector("numeric", length = num_drones)
  coord_y = vector("numeric", length = num_drones)
  dir_x = vector("numeric", length = num_drones)
  dir_y = vector("numeric", length = num_drones)
  local_best_x = vector("numeric", length = num_drones)
  local_best_y = vector("numeric", length = num_drones)
  
  for(i in 1:num_drones){
    coord_x[i] = sample(1:dimension, 1)
    coord_y[i] = sample(1:dimension, 1)
    dir_x[i] = 0
    dir_y[i] = 0
    local_best_x[i] = 0
    local_best_y[i] = 0
  }
  
  return (data.frame(coord_x,coord_y,dir_x,dir_y,local_best_x,local_best_y))
}

#************************************************************************************#
#***** ACTUAL EXECUTION OF PROGRAM *****#
#Create an array of drones

history = c(rep(NA, reps))
for (x in 1:reps){

field = initialize(sparseness, max_particle_conc)
df_drones = init_drones()

i = 0
#saveGIF({
repeat{
  
  if(sum(field != 0) == 0){ #Everything is harvested
    break
  }
  
  df_drones <- scan(df_drones)
  for(n in 1:num_drones){
    if(field[df_drones[n,]$coord_x, df_drones[n,]$coord_y] == 0){ #Nothing to harvest
      df_drones[n,] <- move(df_drones[n,])
    }
    else{ #Harvest the plastic at current patch
      field <- harvest(df_drones[n,])
    }
  }
  i=i+1
  graph(df_drones,i)
  print(i)
}
history[x] = i
#})
}
print(paste("The mean time steps for this optimization setting is: ", mean(history)))
