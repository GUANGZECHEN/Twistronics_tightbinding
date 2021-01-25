#using Geometry_twisted
include("Geometry_twisted.jl")

n_a1=50
n_a2=50
m=10
r=1
r0=[0,0,0]
d=1
lattice,mode="honeycomb","uniform_hop"
# m and r determine twist angle

# build the geometry, it contains position of all sites in the primitive unitcell 
g=Geometry_twisted.get_twisted_lattice_2(lattice,n_a1,n_a2,m,r,d,r0,mode)
println(g.twist_angle*180/pi) # print twist_angle
println(size(g.sites,1))      # print number of sites

Geometry_twisted.remove_layer_g(g,0) # remove layer 0
println(size(g.sites,1))
Geometry_twisted.plot_R(g.sites) # plot the sites 
