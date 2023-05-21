include("gift_wrapping_polygonal_norm.jl")


#############
## sample code for computing 
## set of efficient points 
## following the algorithm of Pelegrin-Fernandez 


########################
## pick norm and display it  
########################

## these are the extremal points of the unit ball
directions = 
[ 
[1;0],
[0;1]
];
directions = normalize_dirs(directions); #symmetrize and normalize
cones = [ [directions[i], directions[i+1]] for i in 1:(length(directions)-1)]

#this decides which facets become `round' 
#isFlats = [1 0 1 0]; #1 for the facets of the ball which are flat
isFlats = [1 1 1 1]; #l^1 norm 

#display the unit ball 
plot_cones(cones,isFlats, plotInterior = true, tt=5) #tt is the line width 


########################
## uniformly sample from a diamond 
########################
num_points = 10^2; #set the number of points
points = [[0.0,0.0] for i in 1:num_points]; #initialize the array
for i in 1:num_points
	x,y = 2*rand()-1,2*rand()-1;
	while(abs(x) + abs(y) > 1) 
		x,y = 2*rand()-1,2*rand()-1;
	end
	points[i] = [x,y]; 
end
points = points[sortperm(map(x->x[2], points), rev = true)]


########################
## draw the diamond 
########################

xs = [0; 1; 0; -1; 0];
ys = [-1; 0; 1; 0; -1]
plot(xs,ys, color = "white", fill = (0, :lightgray), label = false, aspect_ratio = true, axis = :off, ticks = :none)

#####################
#compute the set of efficient points
#####################

#press enter to compute next hull
#type in stop and press enter to exit the loop
markersize = 2; #change to 0 to not display the points

small_p = copy(points)
for num_iters in 1:num_points
	plt, bdry_points = get_hull(small_p, cones, isFlats, true);  
	plot_points!(bdry_points,color = "yellow", markersize = markersize)
	display(plt)
	small_p = setdiff(small_p, bdry_points); 
	small_p = small_p[sortperm(map(x->x[2], small_p), rev = true)]
	if(length(small_p) <= 1) break; end
	if(readline() == "stop") break; end
end
plt2 = plot!()
