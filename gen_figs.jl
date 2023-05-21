##############################
## generate figures in Hamilton-Jacobi scaling limits of Pareto Peeling in 2D 
##############################
include("gift_wrapping_polygonal_norm.jl")
using Random
Random.seed!(1234)




## set parameters
global markersize = 0; #don't plot points

#takes quite awhile about 6000 seconds each
#global num_points = 10^5; #number of points to peel
#global layer_num = 30; #only display every layer_num layers


## a bit faster but still slow
# global num_points = 10^4; #number of points to peel
# global layer_num = 10; #only display every layer_num layers

## super fast but low res
 global num_points = 10^3; #number of points to peel
 global layer_num = 5; #only display every layer_num layers


## change to where you want the pictures to be saved
global file_dir = "./"




######################
## norm ball pictures
######################


l1_directions = 
[ 
[1;0],
[0;1]
];


linf_directions = 
[ 
[1;1],
[-1;1]
];

hex_directions = 
[
[1; sqrt(3)],
[2; 0],
[1; -sqrt(3)]
];

direction_vec = [l1_directions, linf_directions, hex_directions];

if(!isdir("$(file_dir)/norm_balls"))
	mkdir("$(file_dir)/norm_balls")
	mkdir("$(file_dir)/norm_balls/isPartition")
	mkdir("$(file_dir)/norm_balls/noPartition")
end

for i in 1:3, j in 1:2, domain_name in ["circ", "diamond", "square", "q_circ"]
	directions = direction_vec[i]
	directions = normalize_dirs(directions); #symmetrize and normalize
	cones = [ [directions[i], directions[i+1]] for i in 1:(length(directions)-1)]

	isFlats = [[1 1 1 1], [0 1 0 1]][j];
	if(i == 3)
		if(j == 1)
			isFlats = [1 1 1 1 1 1]
		elseif(j == 2)
			isFlats = [0 1 1 0 1 1]
		end
	end

	plt1 = plot_cones(cones,isFlats, plotInterior = false, tt=5); #tt is the line width 
	plt2 = plot_cones(cones,isFlats, plotInterior = true, tt=5); 

	if(j == 1) isFlat_str = "isFlat" else isFlat_str = "isRound"; end
	savefig(plt1, "$(file_dir)/norm_balls/noPartition/norm$(i)_$(isFlat_str)")
	savefig(plt2, "$(file_dir)/norm_balls/isPartition/norm$(i)_$(isFlat_str)")

	########################
	## uniformly sample from the domain 
	########################
	Random.seed!(1)
	points = getSample(num_points, domain_name); ## note this is presorted






	## compute hulls
	small_p = copy(points)
	if(markersize > 0)
		plt = drawDomain(domain_name)
		plot_points!(points,color = "black", markersize = markersize);
	else
		plt = drawDomain(domain_name);
	end

	## plot seperately viscosity solutions and supersolutions for specific choices of norm and domain 
	if(!(isdir("$(file_dir)/$(domain_name)"))) mkdir("$(file_dir)/$(domain_name)"); end

	if(i == 1 && j == 2 && domain_name == "diamond")
		x1s = range(-1,1,length=1000);
		x2s = range(-1,1,length=1000);
		FA(x1,x2) = 1*(abs(x1) + abs(x2) < 0.99)*(1 - abs(x1+x2));
		plt = contour!(x1s,x2s,FA, c=RGB(60/255,60/255,60/255),colorbar = false, levels = 0:0.1:1)
		savefig(plot(plt, dpi=512), "$(file_dir)/$(domain_name)/norm$(i)_$(isFlat_str)_sol");
	elseif(i == 2 && j == 1 && domain_name == "square")

		x1s = range(-1,1,length=1000);
		x2s = range(-1,1,length=1000);
		FB(x1,x2) = 1 - max(abs(x1), abs(x2));
		plt = contour!(x1s,x2s,FB, c=RGB(60/255,60/255,60/255),colorbar = false, levels = 7)
		savefig(plot(plt, dpi=512), "$(file_dir)/$(domain_name)/norm$(i)_$(isFlat_str)_sol");

	end
	#continue;


	println("Computing norm $(i) $(isFlat_str) on $(domain_name)")

	@time for num_iters in 1:num_points
		if(num_iters%layer_num == 0) 		
			plt, bdry_points = get_hull(small_p, cones, isFlats, true);  
		else
			plt, bdry_points = get_hull(small_p, cones, isFlats, false);  
		end
		small_p = setdiff(small_p, bdry_points); 
		small_p = small_p[sortperm(map(x->x[2], small_p), rev = true)]
		if(length(small_p) <= 1) break; end
	end
	if(!(isdir("$(file_dir)/$(domain_name)"))) mkdir("$(file_dir)/$(domain_name)"); end




	savefig(plot(plt, dpi=512), "$(file_dir)/$(domain_name)/norm$(i)_$(isFlat_str)");


	## superimopose viscosity solutions and supersolutions for specific choices of norm and domain 
	 if(i == 1 && j == 2 && domain_name == "diamond")
	 	x1s = range(-1,1,length=1000);
	 	x2s = range(-1,1,length=1000);
	 	FA(x1,x2) = 1*(abs(x1) + abs(x2) <= 1)*(1 - abs(x1+x2));
	 	plt = contour!(x1s,x2s,FA, c=RGB(60/255,60/255,60/255),colorbar = false, levels = 7)
		savefig(plot(plt, dpi=512), "$(file_dir)/$(domain_name)/norm$(i)_$(isFlat_str)_sol_overlay");
	 elseif(i == 2 && j == 1 && domain_name == "square")

	 	x1s = range(-1,1,length=1000);
	 	x2s = range(-1,1,length=1000);
	 	FB(x1,x2) = 1 - max(abs(x1), abs(x2));
	 	plt = contour!(x1s,x2s,FB, c=RGB(60/255,60/255,60/255),colorbar = false, levels = 7)
		savefig(plot(plt, dpi=512), "$(file_dir)/$(domain_name)/norm$(i)_$(isFlat_str)_sol_overlay");
	end




end

#run in bash to trim all the files in the subdirectory
#find -maxdepth 20 -name \*.png -exec bash -c 'convert -trim -resize 712x712 "$0" "$0"' {} \;
#find -maxdepth 20 -name \*.png -exec bash -c 'echo "$0" ' {} \;
