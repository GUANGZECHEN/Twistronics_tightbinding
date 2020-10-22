using LinearAlgebra
using DelimitedFiles
using Statistics
#using Plots
using DelimitedFiles
using Random
using PyCall
using SparseArrays

plt = pyimport("matplotlib.pyplot")

include("Geometry_twisted.jl")


### for triangular lattice, use pi-flux hopping!!!

struct get_hamiltonian
    intra
    a1
    a2
    a12
    a1_2
end

function get_H_twisted_bilayer(R,r,t_xy,t_z,d,lambda,mode,theta,r0)  # d: inter-layer distance, lambda: parameter, mode="uniform_hop" or "pi-flux"
    N=size(R,1)
    R2=[x+r for x in R]

    H=spzeros(Complex,N,N)                        # H_QSL is spin-degen
    #I,J,value=[],[],[]
    for i=1:N
        for j=1:N
            dR=R[i]-R2[j]
            z=dR[3]
            L=norm(dR)
                               
            if abs(z)>0.01
                if abs(t_z*(z^2/L^2)*exp(-lambda*(L-d)))>0.0001
                    H[i,j]=H[i,j]+t_z*(z^2/L^2)*exp(-lambda*(L-d))         # from j to i
                end
            elseif 0.9<L<1.1
                t=get_hopping(R[i],R2[j],t_xy,mode,theta,r0) 
                H[i,j]=H[i,j]+t                           
            end
        end
    end

    return H
end

function get_H_inter_twisted_bilayer(R,inter_vector,t_xy,t_z,d,lambda,mode,theta,r0)
    n_inter=size(inter_vector,1)
    H_inter=[]
    for nn=1:n_inter
        r=inter_vector[nn]
        H=get_H_twisted_bilayer(R,r,t_xy,t_z,d,lambda,mode,theta,r0)
        push!(H_inter,H)
    end
    return H_inter    
end

function get_Hk(R,inter_vector,k,H_inter)
    N=size(R,1)
    n_inter=size(inter_vector,1)
    H=spzeros(Complex,N,N)
    
    for nn=1:n_inter
        r=inter_vector[nn]
        phase=k[1]*r[1]+k[2]*r[2]
        H=H+exp(im*phase)*H_inter[nn]   # R2=R+r, hop=R^dagR2, phase=kr
    end
    if norm(H-adjoint(H))>0.0001
        println("H is non-Hermitian ", norm(H-adjoint(H)))
    end
    return H 
end

function get_hopping(R,R2,t_xy,mode,theta,r0)
    factor=1
    if mode=="uniform_hop"
        if R[3]==1
            factor=1
        end
        return t_xy*factor
    elseif mode=="pi_flux"
        U=rotation_along_z_3D(-(theta))
        if R[3]==1
            R=U*(R-r0)
            R2=U*(R2-r0)
            factor=1
        end
        x=round(R[1]-R2[1],digits=2)
        y=R[2]-R2[2]
        if x*y<-0.01                             # with this method, it is better that the original lattice passes (0,0)
            hop=t_xy
        elseif x>0.1&&y>0.1
            if mod(round(-R[2]*2/sqrt(3)),2)==0  # from B to A, R==A
                hop=-t_xy
            else
                hop=t_xy
            end
        elseif x<-0.1&&y<-0.1
            if mod(round(-R[2]*2/sqrt(3)),2)==0
                hop=t_xy
            else
                hop=-t_xy
            end
        else
            if mod(round(-R[2]*2/sqrt(3)),2)==0
                hop=t_xy
            else
                hop=-t_xy
            end
        end
        return hop*factor                              
    else
        println("invalid mode, default to uniform_hop")
        return t_xy
    end
end

function plot_H(H,R)
    plt.figure(figsize=(6,6),dpi=80)
    N=size(H,1)
    for i=1:N
        for j=1:N
            if real(H[i,j])>0.1
                plt.plot((R[i][1],R[j][1]), (R[i][2],R[j][2]),color="red",linewidth=norm(H[i,j])*1)
            elseif real(H[i,j])<-0.1
                plt.plot((R[i][1],R[j][1]), (R[i][2],R[j][2]),color="blue",linewidth=norm(H[i,j])*1)
            end
        end
    end
    plt.axis([-5,5,-5,5])
    plt.xlabel("")
    plt.ylabel("")
    plt.show()
end

function plot_H_2(H_inter,R,inter_vector)
    plt.figure(figsize=(6,6),dpi=80)
    n=size(H_inter,1)
    
    for i=1:n
        H=H_inter[i]
        r=inter_vector[i]
        R2=[x+r for x in R]

        N=size(H,1)
        for i=1:N
            for j=1:N
                if real(H[i,j])>0.1
                    plt.plot((R[i][1],R2[j][1]), (R[i][2],R2[j][2]),color="red",linewidth=norm(H[i,j])*1)
                elseif real(H[i,j])<-0.1
                    plt.plot((R[i][1],R2[j][1]), (R[i][2],R2[j][2]),color="blue",linewidth=norm(H[i,j])*1)
                end
            end
        end
    end
    plt.axis([-5,5,-5,5])
    plt.xlabel("")
    plt.ylabel("")
    plt.show()
end

