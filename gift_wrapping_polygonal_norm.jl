using Plots,LinearAlgebra
include("helper_functions.jl")



#######################
## Description 
##		computes the set of efficient points with respect to 
##		the norm ball defined by cones, isFlats
## Inputs
## 		points: set of points 
#		cones: array of extremal directions 
##		isFlats: vector of 1s/0s designating a cone as flat - 1 or round - 0
##		isPlot: true if hull is to be plotted as well
## Output
##      returns plot object and points on the boundary of the set
function get_hull(points, cones, isFlats, isPlot = true)
	bdry_points = []
	x0 = points[1]; #assumes presorted and starts with point that has maximum y variable
	x = x0; 
	dir_ind = 1; 
	for num_iters in 1:length(points)
		dir_ind, y = step_2!(x, dir_ind, points, cones, bdry_points, isFlats)
		if(isPlot) step_3(x,y,isFlats[dir_ind],cones[dir_ind]...); end
		x = y; 
		if(x == x0) break; end 
	end

	return plot!(), bdry_points

end



#######################
## Description 
##		plots a collection of cones with varying round and flat parts
##		round parts are red, flat parts are blue
## Inputs
## 		cones/isFlats: define a norm ball as above
##		plotInterior: true if partition of ball is also to be plotted
##		tt: line width of each cone
## Output
##      returns plot object 

function plot_cones(cones, isFlats; plotInterior= true, tt = 1)
	plot(legend = false, aspect_ratio = true, axis = :off, ticks = :none)
	for i in 1:length(cones)
		d1,d2 = cones[i];
		x1 = [0;0]
		x2 = d1;
		x3 = d2; 
		if(plotInterior)
			plot_line!(x1,x2, width = tt)
			plot_line!(x1,x3, width = tt)
		end
		if(isFlats[i] == 1)
			plot_line!(x2,x3, color = "blue", width = tt)
		else
			#draw quadratic bezeir curve
			ts = range(0,1,length=100)
			xs = []
			ys = []
			X = (x1+x2+x3)/3
			for t in ts
				z = X + (1-t)^2*(X-x2) + t^2*(X-x3)
				z = z + (x2+x3)/3
				push!(xs, z[1])
				push!(ys,z[2])
			end
			plot!(xs,ys, color = "red", linewidth = tt)
		end
	end
	return plot!()
end



#######################
## Description 
##		generates a uniform random point on 
##		a ball of radius, radius
## Inputs
## 		radius: a real number specifying the radius  
## Output
##      returns an array of two points on the ball  
function rand_circ(radius = 1)
	tt = [randn();randn()]; tt = tt./norm(tt,2); #random direction
	rr = rand()^(1/2);
	return radius*rr.*(tt);
end





