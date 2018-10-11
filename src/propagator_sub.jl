#using PyCall
#@pyimport matplotlib.pyplot as plt

#
function propagator_sub(u,w,I0,beta,h,skip,Ttotal,L)
# Compute propagators

#beta = 1

N=length(u)

# Set L = 1
#L=1

dz = L/N
dz2 =L^2/N/N

#z=[0:N-1]*dz

# Define nu and nu^2

nu2all=4*(I0+u)

# Find positive I_0 + a, suprathreshold region (necessary for bump case)
Lst = findall(!iszero, nu2all.>0)
Nst = length(Lst)

# Find subthreshold region
Lsub = findall(!iszero, nu2all .<= 0)
Nsub = length(Lsub)

print(Nst,",",Nsub,"\n")

nu2 = zeros(N)

nu2[Lst] = nu2all[Lst]
nu=sqrt.(nu2)

#print(nu2," ",nu[Lst],"\n")

# Initialize propagators

r = zeros(N,N)
s = zeros(N,N)
U = zeros(N,N)

rv = zeros(N,N)
sv = zeros(N,N)
Dvv = Matrix{Float64}(I, size(rand(N,N)))/dz
#Dvv = eye(N)

rp = zeros(N,N)
sp = zeros(N,N)
Dvp = zeros(N,N)
Dvp[:,Lst] = .5*beta*w[:,Lst]   # note factor of 1/2 for causality

Cuu = zeros(N)
Cuu1 = zeros(N)
Cuu2 = zeros(N)

#Construct Backwards Euler update matrix

In = Matrix{Float64}(I, size(rand(Nst,Nst)))
Z = zeros(Nst,Nst)

hp = h/pi;
hb = h*beta
hbN = hb*dz

#B11 = Z
#B12 = In
#B13 = Z
#B21 = -nu2[Lst].*In
#B22 = Z
#B23 = In/pi
#B31 = beta*dz*w[Lst,Lst].*nu[Lst]'
#B32 = Z
#B33 = -beta*In

#B = [[B11 B12 B13], [B21 B22 B23], [B31 B32 B33]]
#print(eigvals(B),"\n")

A11 = In
A12 = -h*In
A13 = Z
A21 = h*nu2[Lst].*In
A22 = In
A23 = -hp*In
A31 = -hbN*w[Lst,Lst].*nu[Lst]'
A32 = Z
A33 = In+hb*In

A = [[A11 A12 A13]; [A21 A22 A23]; [A31 A32 A33]]

Ainv = inv(A)

Asub = Matrix{Float64}(I, size(rand(Nsub, Nsub)))/(1+h*beta)
print(det(Asub),"\n")

phase = zeros(N).-pi

count =0 
Ntotal = round(Ttotal/h)
Cout = zeros(Int(Ntotal/skip)*N,7)


#ind = round(N/2)
ind = Int(round(N/2))

#for k in Lst
#    for l in Lst
#        Cuu1 += Dvv[:,k].*Dvp[:,l]*w[k,l]*nu[l]*beta*dz2/pi*h*skip
#    end
#end

#Cuu2 =  sum(U[:,Lst].*U[:,Lst],2)*dz
#Cuu = Cuu1 .- Cuu2
#print("0"," ",Cuu1[ind]," ",Cuu2[ind]," ",Dvv[ind,ind]," ",Dvp[ind,ind]," ",U[ind,ind],"\n")

# time loop
for n = 1:Ntotal

   t = n*h

# add forcing term to Dvp whenever fictitious oscillator fires
    for k in Lst  
	phase[k] +=nu[k]*h
	if phase[k] >= pi
          for m=1:N
           Dvp[m,k] += beta*w[m,k]*exp(-beta*(phase[k]-pi)/nu[k])
          end
#           Dvp[Lst,k] += beta*w[Lst,k]
	   phase[k] -= 2*pi
#	   print(phase[k]," ","*","\n")
	end
    end

    for j in Lst
    # z' loop
    
    	#Construct update vectors

    	v1j = [r[Lst,j];s[Lst,j];U[Lst,j]+h*beta/2/pi*w[Lst,j]*nu[j]]
    	v2j = [rv[Lst,j];sv[Lst,j];Dvv[Lst,j]]    
	v3j = [rp[Lst,j];sp[Lst,j];Dvp[Lst,j]]

	# Backwards Euler step
	v1j = Ainv*v1j
	v2j = Ainv*v2j
	v3j = Ainv*v3j

	
	# Load propagators
	r[Lst,j] = v1j[1:Nst]
	s[Lst,j] = v1j[Nst+1:2*Nst]
    	U[Lst,j] = v1j[2*Nst+1:3*Nst]
    	rv[Lst,j] = v2j[1:Nst]
	sv[Lst,j] = v2j[Nst+1:2*Nst]
    	Dvv[Lst,j] = v2j[2*Nst+1:3*Nst]
    	rp[Lst,j] = v3j[1:Nst]
	sp[Lst,j] = v3j[Nst+1:2*Nst]
    	Dvp[Lst,j] = v3j[2*Nst+1:3*Nst]

	# Subthreshold update step
        U[Lsub,j] =   Asub*(U[Lsub,j]+h*beta*w[Lsub,Lst]*(nu[Lst].*r[Lst,j])*dz+h*beta/2/pi*w[Lsub,j]*nu[j])
        #U[Lsub,j] =   Asub*(U[Lsub,j]+h*beta*w[Lsub,Lst]*(nu[Lst].*r[Lst,j])*dz)
        Dvv[Lsub,j] = Asub*(Dvv[Lsub,j]+h*beta*w[Lsub,Lst]*(nu[Lst].*rv[Lst,j])*dz)
        Dvp[Lsub,j] = Asub*(Dvp[Lsub,j]+h*beta*w[Lsub,Lst]*(nu[Lst].*rp[Lst,j])*dz)

    end

#    Uf = h*w[Lst[1]-1,Lst]*(nu[Lst].*r[Lst,Lst[1]])*dz

Dvv[Lsub,Lsub]=Asub*Dvv[Lsub,Lsub]

for k=1:N
    for l in Lst
            for m = 1:N
            Cuu1[m] += Dvv[m,k]*Dvp[m,l]*w[k,l]*nu[l]*beta*dz2/pi*h
            end
#                   Cuu1 += Dvv[:,k].*Dvp[:,l]*w[k,l]*nu[l]*beta*dz2/pi*h*skip
    end
end

   if (rem(n,skip) == 0)

           # Update integral for first diagram of <uu>
           # Use time translational invariance, integrating tau is he same as integrating t (in reverse)
           # use rho(z,pi) = nu/4/pi

          Cuu2 =  sum(U[:,Lst].*U[:,Lst],dims=2)*dz
          Cuu = Cuu1 .- Cuu2

#     print(size(Cuu1)," ",size(Cuu2),"\n")
for k=1:N
      Cout[count*N+k,1] = t
      Cout[count*N+k,2]=Cuu[k]
      Cout[count*N+k,3]=Cuu1[k]
      Cout[count*N+k,4]=Cuu2[k]
      Cout[count*N+k,5]=Dvv[k,ind]
      Cout[count*N+k,6]=Dvp[k,ind]
      Cout[count*N+k,7]=U[k,ind]
end
      print(t," ",Cuu[ind]," ",Cuu2[ind]," ",Dvv[ind,ind]," ",Dvp[ind,ind]," ",U[ind,ind],"\n")
      count+=1
    end
end

# Cuu1 = 2*2\beta \int dz1 dz2 dtau Dvv(z,t;z1,tau) Dvp(z,t;z2,tau) w(z1-z2) rho(z2)
# Cuu1  *= beta*dz2*h/pi  # need the 1/4/pi factor to convert nu to rho.

Cuu2 =  sum(U[:,Lst].*U[:,Lst],dims=2)*dz
Cuu = Cuu1 .- Cuu2

#return Cout,Cuu,Cuu1,Cuu2,r,U,rv,Dvv,rp,Dvp,Lst,Lsub,nu
return Cout, Dvp
end
