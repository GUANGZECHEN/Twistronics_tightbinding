using LinearAlgebra
using DelimitedFiles
using Statistics
#using Plots
using DelimitedFiles

using PyCall
#using PyPlot

plt = pyimport("matplotlib.pyplot")

mutable struct get_geometry_twisted
    lattice
    sites
    unit_vector
    dimension
    angle
end

function get_theta(m,r,a1,a2)
    theta=acos((3*m^2+3*m*r+r^2/2)/(3*m^2+3*m*r+r^2))
    A1=m*a2+(m+r)*a1
    A2=-(m+r)*a2+(2*m+r)*a1                                      # the a1,a2 in literature is the opposite definition from here
    return theta, A1, A2
end

function get_lattice_3D_rotated(lattice,n_a1,n_a2,z,theta,r0)   #theta in units of 1, z: layer-index  #add confinement in unitcell, r0: parallel alignment of layer 
    if lattice=="triangular"
        a1=[1,0,0]
        a2=[1/2,-sqrt(3)/2,0]
        R=geometry_simple_3D_rotated(n_a1,n_a2,a1,a2,z,theta,r0)
    elseif lattice=="square"
        a1=[1,0,0]
        a2=[0,1,0]
        R=geometry_simple_3D_rotated(n_a1,n_a2,a1,a2,z,theta,r0)
    elseif lattice=="honeycomb"
        a1=[sqrt(3),0,0]
        a2=[sqrt(3)/2,-3/2,0]
        a3=[sqrt(3)/2,-1/2,0]
        R=geometry_bipartite_3D_rotated(n_a1,n_a2,a1,a2,a3,z,theta,r0)
    else
        println("invalid lattice")
    end
    
    return R,a1,a2
end

function rotation_along_z_3D(theta)   
    U=[[cos(theta),sin(theta),0] [-sin(theta),cos(theta),0] [0,0,1]]     #in this way U=anti-clock rotation of theta
    return U
end

function geometry_simple_3D_rotated(n,m,a1,a2,z,theta,r0)
    U=rotation_along_z_3D(theta)
    a1=U*a1
    a2=U*a2
    R0=[0,0,z]-(n-1)/2*a1-(m-1)/2*a2+r0                                  # s.t. the second layer is rotated w.r.t. [0,0]
    R=[]
    for i=0:n-1
        for j=0:m-1
            push!(R,i*a1+j*a2+R0)
        end
    end
    return R
end

function geometry_bipartite_3D_rotated(n,m,a1,a2,a3,z,theta,r0)
    U=rotation_along_z_3D(theta)
    a1=U*a1
    a2=U*a2
    a3=U*a3
    R0=[0,0,z]-(n-1)/2*a1-(m-1)/2*a2+r0
    R=[]
    for i=0:n-1
        for j=0:m-1
            for k=0:1
                push!(R,i*a1+j*a2+k*a3+R0)
            end
        end
    end
    return R
end

function get_inter_vector_3D(n,m,a1,a2)
    inter_vector=[[0,0,0],n*a1,-n*a1,m*a2,-m*a2,m*a2-n*a1,n*a1-m*a2,m*a2+n*a1,-n*a1-m*a2]
    return inter_vector
end

function nearest(site_a,site_b)
    return 0.1<norm(site_a-site_b)<1.1
end

function next_nearest(site_a,site_b)
    return 1.1<norm(site_a-site_b)<1.9
end


function next_next_nearest(site_a,site_b)
    return 1.9<norm(site_a-site_b)<2.1
end

function plot_R(R)
    plt.figure(figsize=(6,6),dpi=80)
    N=size(R,1)
    x=[]
    y=[]
    x2=[]
    y2=[]
    for i=1:N
        if R[i][3]==0
            push!(x,R[i][1])
            push!(y,R[i][2])
        else
            push!(x2,R[i][1])
            push!(y2,R[i][2])
        end
    end
    #println(size(x,1),size(y,1))
    plt.scatter(x,y,color="red",s=10)
    plt.scatter(x2,y2,color="blue",s=10)
    plt.axis([-40,40,-40,40])
    plt.xlabel("")
    plt.ylabel("")
    plt.show()
end

function add_vacancy(R,vacant_site)
    deleteat!(R,vacant_site)
    return R
end

function merge_R(R,R2)
    N=size(R2,1)
    for i=1:N
        push!(R,R2[i])
    end
    return R
end

function check_site_in_unit_cell(site,A1,A2)      #check whether site is in unitcell [-1/2,1/2]*A1+[-1/2,1/2]*A2, site can be 3D
    num=A1[1]*A2[2]-A1[2]*A2[1]
    num2=site[1]*A1[2]-site[2]*A1[1]
    num3=site[1]*A2[2]-site[2]*A2[1]              
    alpha=num3/num
    beta=-num2/num                                # site=alpha*a1+beta*a2
    if -1/2<round(alpha,digits=10)<=1/2 && -1/2<round(beta,digits=10)<=1/2
        return true
    else 
        return false
    end
end

function get_unit_cell(R,A1,A2)
    R_unitcell=[]
    N=size(R,1)
    n1=0
    n2=0
    for i=1:N
        if check_site_in_unit_cell(R[i],A1,A2)
            push!(R_unitcell,R[i])
            if R[i][3]==0
                n1=n1+1
            elseif R[i][3]==1
                n2=n2+1
            end
        end
    end

    return R_unitcell,n1,n2
end

function check_R(R_unitcell,inter_vector,R)
    n_inter=size(inter_vector,1)
    R=[round.(x,digits=5) for x in R]
    #println(R)
    N=size(R_unitcell,1)
    for ii=2:n_inter
        r=inter_vector[ii]
        for i=1:N
            R0=R_unitcell[i]+r
            R0=round.(R0,digits=5)
            #println(R0)
            if !(R0 in R)
                println(R0)
                #println(R)
                return false
            end
        end
    end

    return true
end

function plot_R_unitcell(R_unitcell,inter_vector)
    n_inter=size(inter_vector,1)
    N=size(R_unitcell,1)
    for ii=2:n_inter
        r=inter_vector[ii]
        for i=1:N
            push!(R_unitcell,R_unitcell[i]+r)
        end
    end
    plot_R(R_unitcell)
end

function get_twisted_lattice(BC,lattice,n_a1,n_a2,m,r,d,r0)
    R,a1,a2=get_lattice_3D_rotated(lattice,n_a1,n_a2,0,0,[0,0,0])
    theta, A1, A2=get_theta(m,r,a1,a2)
    R2,a1_2,a2_2=get_lattice_3D_rotated(lattice,n_a1,n_a2,d,theta,r0)
    R=merge_R(R,R2)
    R_unitcell,n1,n2=get_unit_cell(R,A1,A2)
    #plot_R(R)
    if !(n1==n2)
        println("number of sites on different layers mismatch: ", n1," ",n2)
    end

    if BC=="PBC"
        inter_vector=get_inter_vector_3D(1,1,A1,A2)
    elseif BC=="OBC"
        inter_vector=[[0,0,0]]
    else
        println("invalid bundary condition, using OBC by default")
        inter_vector=[[0,0,0]]
    end

    if !(check_R(R_unitcell,inter_vector,R))
        println("wrong unitvector, superstructure not periodic; or R not large enough")
    end
    
    return R_unitcell, inter_vector, A1, A2, theta
end

function get_reciprocal_vector(A1,A2) 
    A3=[0,0,1]
    B1=2*pi*cross(A2,A3)/dot(A1,cross(A2,A3))
    B2=2*pi*cross(A1,A3)/dot(A2,cross(A1,A3))
    b1=[B1[1],B1[2]]
    b2=[B2[1],B2[2]]
    return b1,b2
end
#n_a1=200
#n_a2=200
#m=20
#r=1
#r0=[0,0.5,0]
#d=1
#geometry="triangular"
#BC="PBC"

#R_unitcell, inter_vector=get_twisted_lattice(BC,geometry,n_a1,n_a2,m,r,d,r0)

#plot_R_unitcell(R_unitcell,inter_vector)




