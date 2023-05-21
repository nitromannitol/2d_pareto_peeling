#helper functions for running the Pelegrin-Fernandez gift wrapping algorithm

## gift wrapping algorithm to compute the boundary
## of a polygonal norm
function plot_points(points; color = "red", markersize = 1)
	if(length(points) == 0) 
		plot(label = false, aspect_ratio = true, axis = :off, ticks = :none)
	else
		px = map(x->x[1], points); 
		py = map(x->x[2], points); 
		scatter(px, py, color = color, label = false, aspect_ratio = true, axis = :off, ticks = :none, markersize = markersize)
	end
end
function plot_points!(points; color = "red", markersize=1)
	if(length(points) == 0) 
		plot!(label = false, aspect_ratio = true, axis = :off, ticks = :none)
	else
		px = map(x->x[1], points); 
		py = map(x->x[2], points); 
		scatter!(px, py, color = color, label = false, aspect_ratio = true, axis = :off, ticks = :none, markersize = markersize)
	end
end
function plot_line!(x1,x2; color = "black", width = 1)
	plot!([x1[1], x2[1]], [x1[2], x2[2]], color = color, label = false, width = width)
end

function plot_directions(directions)
	plot(label = false, aspect_ratio = true, axis = :off, ticks = :none)
	for i in 1:length(directions)
		x1 = [0;0]; x2 = directions[i]
		plot_line!(x1,x2)
	end
	return plot!()
end



#projects y onto line 
#x + t d
function proj(y,d,x)
	z = (y-x); 
	#normalize direction
	d = d/norm(d,2); 

	#project onto line
	z = z - (z'*d)*d; 
	return z
end

## distance from y to the 
## line x + t*d
## assumes |d| > 0
function d_e(y, d, x)
	return norm(proj(y,d,x),2)
end

#check if x is above a line
function isAbove(z, x, d)
	return sign(det( [d x-z]))>=0;
end

#check if x is strictly above a line
function isOn(z, x, d)
	return isapprox(det( [d x-z]),0);
end

#get points within a cone
function get_points(x, cone, points)
	d1,d2 = cone
	dir_points = Array{Array{Float64,1},1}()
	for i in 1:length(points)
		p = points[i]; 
		if(p == x) continue; end
		if(isAbove(p,x,d1) && !isAbove(p,x,d2) || isOn(p,x,d1) || isOn(p,x,d2))
			push!(dir_points, p)
		end
	end
	return dir_points
end




### computes the next point in the gift wrapping
## algorithm
function step_2!(x, dir_ind, points, cones, bdry_points, isFlats)
	if(dir_ind > length(cones)) 
		dir_ind = 1; 
	end
	cone_points = Array{Array{Float64,1},1}()
	for i in (dir_ind):length(cones)
		cone_points = get_points(x, cones[i],points)
		if(length(cone_points) > 0) dir_ind = i; break; end 
	end
	if(length(cone_points) == 0) 
		for i in 1:dir_ind
			cone_points = get_points(x, cones[i],points)
			if(length(cone_points) > 0) dir_ind = i; break; end 
		end
	end

	
	if(length(cone_points) == 0) error("did not find points"); end; #this should not be possible
	d1,d2 = cones[dir_ind]
	isFlat = isFlats[dir_ind]

	if(isFlat==1)

		#compute distances to the line
		dist_vec = fill(0.0, length(cone_points));
		for i in 1:length(cone_points)
			dist_vec[i] = d_e(cone_points[i], d1, x);
		end

		cone_points = cone_points[dist_vec.==minimum(dist_vec)]


		push!(bdry_points, cone_points...)
	elseif(isFlat == 0)
		dist_vec = fill(0.0, length(cone_points));
		for i in 1:length(cone_points)
			z1 = d_e(cone_points[i], d1, x);
			z2 = d_e(cone_points[i], d2, x);
			dist_vec[i] = z1/z2;
		end
		cone_points = cone_points[dist_vec.==minimum(dist_vec)]
		if(isFlat==0)
			push!(bdry_points, cone_points...)
		end
	end
	
	return dir_ind, cone_points[1]
end

### plots the lines
function step_3(x,y, isFlat, d1, d2)

	if(isFlat == 1)
		## store shortest route from x to y going parallel to di+1 then then di
		x1 = x
		
		a,b = d2;
		c,d = d1; 

		t = inv([d2 d1])*(x-y)
		#these are the same
		x2 = x - t[1]*d2
		x2 = y +t[2]*d1

		x3 = y; 
		plot_line!(x1,x2, color = "blue")
		plot_line!(x2,x3, color = "blue")
	else
		plot_line!(x,y,color="red")
	end
end



function signed_angle(a,b)
	an_ = acos(dot(a,b)/(norm(a)*norm(b)))
	if( a[1]*b[2] - a[2]*b[2] < 0)
		return -an_
	else
		return an_
	end

end

#convert a set of directions to complex coordinates 
function compute_angles(directions)
	c_dirs = [x[1]+im*x[2] for x in directions]
	return angle.(c_dirs);

end

#symmetrizes the set of directions
#normalizes
#returns it in clockwise order
function normalize_dirs(directions)
	directions = map(x-> x./norm(x,2), directions)
	directions = [directions..., .- directions...]
	directions = unique(directions)
	angles = compute_angles(directions)
	directions = directions[sortperm(angles, rev = true)]
	angles = angles[sortperm(angles, rev=true)] 
	ind = 1
	directions = [directions[ind:end]..., directions[1:ind-1]...]

	return [directions..., directions[1]]

end



#### functions for generating random points 

##generate a random point on the ball
function rand_circ(radius = 1)
	tt = [randn();randn()]; tt = tt./norm(tt,2); #random direction
	rr = rand()^(1/2);
	return radius*rr.*(tt);
end

## generates random points from a grid 
function gen_rand_points(N, num_points)
	xs = range(-1,1,length=N);
	ys = range(-1,1,length=N); 

	#select points from which to find the hull
	points = Array{Array{Float64,1},1}()
	px = []; py = []
	for i in 1:num_points
		while(true)
			ii = rand(1:N)
			jj = rand(1:N)
			x = xs[ii];
			y = ys[jj];
			push!(points,[x; y])
			push!(px, x)
			push!(py, y);
		end
	end
	## sort points by y coordinate 
	points = points[sortperm(py, rev = true)]
	px = px[sortperm(py, rev = true)]
	py = py[sortperm(py, rev = true)]
	return points; 
end

## generates random points from a domain
function getSample(num_points, domain_name); 
	if(domain_name == "square")
		points = [[2*rand()-1, 2*rand()-1] for i in 1:num_points]
	elseif(domain_name == "circ")
		points = [rand_circ() for i in 1:num_points]
	elseif(domain_name == "diamond")
		points = [[0.0,0.0] for i in 1:num_points]
		for i in 1:num_points
			x,y = 2*rand()-1,2*rand()-1;
			while(abs(x) + abs(y) > 1) 
				x,y = 2*rand()-1,2*rand()-1;
			end
			points[i] = [x,y]; 
		end

	elseif(domain_name == "q_circ")
		 points = [rand_circ() for i in 1:num_points]
		 for i in 1:num_points
		 	x,y = rand_circ();
		 	while(x <= 0 || y <= 0)
		 		x,y = rand_circ();
		 	end
		 	points[i] = [x,y]; 
		 end
	elseif(domain_name == "l")
		points = [[0.0,0.0] for i in 1:num_points]
		for i in 1:num_points
			x,y = 2*rand()-1,2*rand()-1;
			while(x >=0 && y>=0) 
				x,y = 2*rand()-1,2*rand()-1;
			end
			points[i] = [x,y]; 
		end
	elseif(domain_name == "rounded_l")
		points = [[0.0,0.0] for i in 1:num_points]
		for i in 1:num_points
			x,y = 2*rand()-1,2*rand()-1;
			while(x >=0 && y>=0 || x <=0 && y <= 0 && sqrt(x^2+y^2) > 1)
				x,y = 2*rand()-1,2*rand()-1;
			end
			points[i] = [x,y]; 
		end

	else
		error("domain name not recognized")
	end


	## sort points by y coordinate 
	return points[sortperm(map(x->x[2], points), rev = true)]
end


function drawDomain(domain_name)
	if(domain_name == "square")
		xs = [-1; 1; 1; -1]
		ys = [1; 1; -1; -1]

	elseif(domain_name == "circ")
		thetas = range(0,2*pi,length=1000)
		xs = fill(0.0,length(thetas)); ys = fill(0.0,length(thetas));
		for i in 1:length(thetas)
			theta = thetas[i]; 
			xs[i] = cos(theta); 
			ys[i] = sin(theta);
		end
	elseif(domain_name == "diamond")
		xs = [0; 1; 0; -1; 0];
		ys = [-1; 0; 1; 0; -1]
	elseif(domain_name == "q_circ")
		thetas = range(0,2*pi,length=1000)
		xs = fill(0.0,length(thetas)); ys = fill(0.0,length(thetas));
		for i in 1:length(thetas)
			theta = thetas[i]; 
			if(cos(theta) <=0 || sin(theta) <= 0) continue; end
			xs[i] = cos(theta); 
			ys[i] = sin(theta);
		end
	elseif(domain_name == "l")
		xs = [-1; 0; 0; 1; 1; 0.0; -1]
		ys = [1;1;0;0;-1;-1.0; -1.0]
	elseif(domain_name == "rounded_l")
		xs = [-1; 0; 0; 1; 1; 0.0]
		ys = [1;1;0;0;-1;-1.0]
		thetas = reverse(range(pi,3*pi/2,length=1000))
		for i in 1:length(thetas)
			theta = thetas[i];
			push!(xs, cos(theta)); push!(ys, sin(theta));
		end
		## add circle
		push!(xs, -1); push!(ys, -1); 

	else
		error("domain name not recognized")
	end
	return plot(xs,ys, color = "white", fill = (0, :lightgray), label = false, aspect_ratio = true, axis = :off, ticks = :none);
end


## test sampling functions
# for domain_name in ["square", "circ", "diamond", "q_circ", "l", "rounded_l"]
# 	println(domain_name)
# 	plt = drawDomain(domain_name);
# 	points = getSample(1090, domain_name);
# 	display(plot_points!(points,color = "black", markersize = 1))
# 	if(readline() == "stop") break; end
# end
