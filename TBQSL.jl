using LinearAlgebra
using DelimitedFiles
using Statistics
#using Plots
using DelimitedFiles
using Random
using SparseArrays
using Arpack

include("Geometry_twisted.jl")
include("tight-binding_twisted.jl")
using PyCall
plt = pyimport("matplotlib.pyplot")

using Dates

function main()
    n_a1=99
    n_a2=99                # these two create a large system that allows choosing a unitcell inside
    m=34
    r=2
    r0=[0,0,0]              # alignment of two layers
    d=1                     # distance between two layers, in units of lattice constant
    lattice,mode="triangular","pi_flux"
    BC="PBC"

    R_unitcell, inter_vector, A1, A2, theta=get_twisted_lattice(BC,lattice,n_a1,n_a2,m,r,d,r0)
    b1,b2=get_reciprocal_vector(A1,A2)
    K_Gamma=2*b1+b2         # vector from one Gamma point to the other, along Gamma-K-M-K'-Gamma, only for angle(a1,a2)=60degree
    #plot_R_unitcell(R_unitcell,inter_vector)
    println(theta*180/pi)
    println(size(R_unitcell,1))

    t_xy=-1
    t_z=0
    lambda=10
    
    n_k=100
    momenta=[]
    energies=[]
    n_bands=40
    t0=now()
    H_inter=get_H_inter_twisted_bilayer(R_unitcell,inter_vector,t_xy,t_z,d,lambda,mode,theta,r0)
    #plot_H(H_inter[1],R_unitcell)
    t1=now()
    println("time to construct H: ", t1-t0)
    for i=0:n_k
        println("summing over k point: ",i)
        k=i/n_k*K_Gamma
        
        H=get_Hk(R_unitcell,inter_vector,k,H_inter)
  
        eigvals=eigs(H,nev=n_bands,sigma=0.01,which=:LM,maxiter=10000)[1]
        #F=eigen(Matrix(H))
        #eigvals=F.values
        #eigvals=real(eigvals)
        #n_bands=size(eigvals,1)
        for i_band=1:n_bands
            push!(momenta,i)
            push!(energies,eigvals[i_band])
        end
    end
    t3=now()
    println("time to compute eigs: ", t3-t1)

    plt.figure(figsize=(6,6),dpi=80)
    plt.scatter(momenta, energies)
    plt.axis([-1,n_k+1,-0.25,0.25])
    plt.xlabel("k")
    plt.ylabel("E")
    plt.show()
end

main()
