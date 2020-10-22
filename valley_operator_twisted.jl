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

function get_valley_operator_twisted_bilayer(R,r,mode,theta,r0)  # d: inter-layer distance, lambda: parameter, mode="uniform_hop" or "pi-flux"
    N=size(R,1)
    R2=[x+r for x in R]

    P=spzeros(Complex,N,N)                        # H_QSL is spin-degen
    #I,J,value=[],[],[]
    for i=1:N
        for j=1:N
            dR=R[i]-R2[j]
            z=dR[3]
                               
            if abs(z)<0.01
                t=get_valley_operator_hopping(R[i],R2[j],mode,theta,r0) 
                P[i,j]=P[i,j]+t                           
            end
        end
    end
    return P
end

function get_valley_operator_hopping(R,R2,mode,theta,r0)
    layer_factor=1
    U=rotation_along_z_3D(-(theta))
    if R[3]==1
        R=U*(R-r0)
        R2=U*(R2-r0)
        layer_factor=1
    end
    x=R[1]-R2[1]
    y=R[2]-R2[2]
    if mode=="uniform_hop"    #valley operator for honeycomb
        hop=0
        if round(x/sqrt(3),digits=3)==1 && round(y,digits=3)==0
            hop=im
        elseif round(x/sqrt(3),digits=3)==-1 && round(y,digits=3)==0
            hop=-im
        elseif round(x*2/sqrt(3),digits=3)==-1 && round(y,digits=3)==3/2
            hop=im
        elseif round(x*2/sqrt(3),digits=3)==1 && round(y,digits=3)==-3/2
            hop=-im
        elseif round(x*2/sqrt(3),digits=3)==-1 && round(y,digits=3)==-3/2
            hop=im
        elseif round(x*2/sqrt(3),digits=3)==1 && round(y,digits=3)==3/2
            hop=-im
        end
        return hop*layer_factor/(3*sqrt(3))
    elseif mode=="pi_flux"
        hop=0
        if round(x,digits=3)==1 && round(y,digits=3)==0
            hop=-im
        elseif round(x,digits=3)==-1 && round(y,digits=3)==0
            hop=im
        elseif round(x,digits=3)==0 && round(y/sqrt(3),digits=3)==-1
            hop=-im
        elseif round(x,digits=3)==0 && round(y/sqrt(3),digits=3)==1
            hop=im
        end
        return hop*layer_factor/4                              
    else
        println("invalid mode for valley operator")
    end
end

function get_valley_operator_inter_twisted_bilayer(R,inter_vector,mode,theta,r0)
    n_inter=size(inter_vector,1)
    P_inter=[]
    for nn=1:n_inter
        r=inter_vector[nn]
        P=get_valley_operator_twisted_bilayer(R,r,mode,theta,r0)
        push!(P_inter,P)
    end
    return P_inter    
end

function get_Pk(R,inter_vector,k,P_inter)
    N=size(R,1)
    n_inter=size(inter_vector,1)
    P=spzeros(Complex,N,N)
    
    for nn=1:n_inter
        r=inter_vector[nn]
        phase=k[1]*r[1]+k[2]*r[2]
        P=P+exp(im*phase)*P_inter[nn]   # R2=R+r, hop=R^dagR2, phase=kr
    end

    return P 
end

function plot_H(H,R)
    plt.figure(figsize=(6,6),dpi=80)
    N=size(R,1)
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

function plot_P(P,R)
    plt.figure(figsize=(6,6),dpi=80)
    N=size(R,1)
    for i=1:N
        for j=i:N
            if imag(P[i,j])>0.1
                plt.plot((R[i][1],R[j][1]), (R[i][2],R[j][2]),color="red",linewidth=norm(P[i,j])*1)
            elseif imag(P[i,j])<-0.1
                plt.plot((R[i][1],R[j][1]), (R[i][2],R[j][2]),color="blue",linewidth=norm(P[i,j])*1)
            end
        end
    end
    plt.axis([-5,5,-5,5])
    plt.xlabel("")
    plt.ylabel("")
    plt.show()
end


function get_layer_operator_twisted_bilayer(R)  # d: inter-layer distance, lambda: parameter, mode="uniform_hop" or "pi-flux"
    N=size(R,1)

    L=spzeros(Complex,N,N)   
    for i=1:N
        if R[i][3]==0
            L[i,i]=1
        else
            L[i,i]=-1
        end
    end
    return L
end