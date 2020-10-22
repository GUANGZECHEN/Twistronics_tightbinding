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
include("valley_operator_twisted.jl")
using PyCall
plt = pyimport("matplotlib.pyplot")
xlsxwriter=pyimport("xlsxwriter")
colormap=pyimport("matplotlib.colors")

using Dates

function main(m,t_z)
    n_a1=99
    n_a2=99                # these two create a large system that allows choosing a unitcell inside, pay attention to choose them large enough
    #m=50

    r=1
    r0=[0,0,0]              # alignment of two layers
    d=1                     # distance between two layers, in units of lattice constant
    #lattice,mode="honeycomb","uniform_hop"
    lattice,mode="triangular","pi_flux"
    #lattice,mode="triangular","uniform_hop"
    BC="PBC"

    R_unitcell, inter_vector, A1, A2, theta=get_twisted_lattice_2(BC,lattice,n_a1,n_a2,m,r,d,r0)
    b1,b2=get_reciprocal_vector(A1,A2)
    #theta2=pi/3      #there are 3 possible ways that two pi-flux can overlap
    #theta=theta+theta2
    
    K_Gamma=b2         # vector from one Gamma point to the other, along Gamma-K-M-K'-Gamma, only for angle(a1,a2)=60degree
    k0=b1/2+b2/4
    
    #plot_R_unitcell(R_unitcell,inter_vector)
    println(theta*180/pi)
    println(size(R_unitcell,1))

    t_xy=-1
    #t_z=0.1
    lambda=10
    
    n_k=100
    e_max=0.2
    n_e=100

    n_bands=40
    t0=now()
    H_inter=get_H_inter_twisted_bilayer(R_unitcell,inter_vector,t_xy,t_z,d,lambda,mode,theta,r0)
    P_inter=get_valley_operator_inter_twisted_bilayer(R_unitcell,inter_vector,mode,theta,r0)
    #H_inter=H_inter.+0.01*P_inter
    #P=get_valley_operator_twisted_bilayer(R_unitcell,[0,0,0],mode,theta,r0)
    #plot_P(P,R_unitcell)

    #plot_H_2(H_inter,R_unitcell,inter_vector)
    #plot_H(H_inter[1],R_unitcell)
    #get_band_twisted(R_unitcell,inter_vector,H_inter,b1,b2,n_k,n_bands,P_inter)
    
    t1=now()
    println("time to construct H: ", t1-t0)
    energies=get_eigvals_twisted(R_unitcell,inter_vector,H_inter,b1,b2,n_k,n_bands)
    t3=now()
    println("time to compute eigs: ", t3-t1)

    eta=0.005
    omegas=[]
    DOSs=[]
    for i=0:n_e
        omega=e_max*(-1+2*i/n_e)
        DOS=get_DOS(energies,omega,eta)
        push!(omegas,omega)
        push!(DOSs,DOS)
    end
    return omegas, DOSs

    #plt.figure(figsize=(6,6),dpi=80)
    #plt.scatter(omegas, DOSs)
    #plt.axis([-0.21,0.21,0,1000])
    #plt.xlabel("k")
    #plt.ylabel("E")
    #plt.show()
end

function get_band_twisted(R_unitcell,inter_vector,H_inter,b1,b2,n_k,n_bands,P_inter)
    momenta=[]
    energies=[]
    valley_polarization=[]
    for i=0:n_k
        println("summing over k point: ",i)
        k=b2*i/n_k+(b1-b2/2)/2#i/n_k*(2*b1+b2)####-b1/2+i/n_k*b1###
        #if i<=n_k/2
       #     k=i/(n_k/2)*(4*b2/3+2*b1/3)
       # else
        #    k=(4*b2/3+2*b1/3)+(i-n_k/2)/(n_k/2)*(-2*b1/3+2*b2/3)
        #end
        
        H=get_Hk(R_unitcell,inter_vector,k,H_inter)
        P=get_Pk(R_unitcell,inter_vector,k,P_inter)#get_layer_operator_twisted_bilayer(R_unitcell)#
        #H=H+0.05*P
        eigvals,eigvecs=eigs(H,nev=n_bands,sigma=0.01,which=:LM,maxiter=10000)
        #F=eigen(Matrix(H))
        #eigvals,eigvecs=F.values,F.vectors
        #eigvals=real(eigvals)
        #n_bands=size(eigvals,1)
        for i_band=1:n_bands
            push!(momenta,i)
            push!(energies,eigvals[i_band])
            psi=eigvecs[:,i_band]
            push!(valley_polarization,real(adjoint(psi)*P*psi))
        end
    end
    plt.figure(figsize=(6,6),dpi=80)
    colors=[(1,0,0),(0,1,0),(0,0,1)]
    cm=colormap.LinearSegmentedColormap.from_list("my_list", colors)
    sc=plt.scatter(momenta, energies,c=valley_polarization,cmap=cm,vmin=-1,vmax=1)
    plt.colorbar(sc)
    plt.axis([-1,n_k+1,-0.2,0.2])
    plt.xlabel("k")
    plt.ylabel("E")
    plt.show()
    return momenta, energies
end

function get_eigvals_twisted(R_unitcell,inter_vector,H_inter,b1,b2,n_k,n_bands)
    energies=[]
    for i=1:n_k
        for j=1:n_k
            println("summing over k point: ",i," ",j)
            k=i/n_k*b1+j/n_k*b2
        
            H=get_Hk(R_unitcell,inter_vector,k,H_inter)
  
            eigvals=eigs(H,nev=n_bands,sigma=0.01,which=:LM,maxiter=10000)[1]

            for i_band=1:n_bands
                push!(energies,eigvals[i_band])
            end
        end
    end
    return energies
end

function get_DOS(energies,omega,eta)
    DOS=0
    for E in energies
        DOS=DOS+1/(E-omega+im*eta)
    end
    DOS=-1/pi*imag(DOS)
    return DOS
end

omegas=[]
DOSs=[]
tzs=[]
m=1
for i=0:0
    t_z=0.1+0.1*i
    println("computing t_z=", t_z)
    omega,DOS=main(m,t_z)
    n=size(omega,1)
    for ii=1:n
        push!(omegas,omega[ii])
        push!(DOSs,DOS[ii])
        push!(tzs,t_z)
    end
    #plt.figure(figsize=(6,6),dpi=80)
    #sc=plt.scatter(omegas, DOSs)
    #plt.xlabel("k")
    #plt.ylabel("E")
    #plt.show()
end

workbook = xlsxwriter.Workbook(string("DOS_vs_tz_m=35.xlsx"))
worksheet = workbook.add_worksheet()
N=size(omegas,1)
for row=0:N-1
    worksheet.write(row, 0, omegas[row+1])
    worksheet.write(row, 1, tzs[row+1])
    worksheet.write(row, 2, DOSs[row+1])
end      
workbook.close()
