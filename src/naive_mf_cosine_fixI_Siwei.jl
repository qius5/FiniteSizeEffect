#using PyCall
#@pyimport matplotlib.pyplot as plt

# gain function
gain(x) = sqrt.(max(0,x))

#function weights(N,L,A,B,a,b)
function weights(N,L,A,a,b)
        # Construct periodic weight function

        w=zeros(N,N)
        w1=zeros(N)
        dz=L/N
        #w1=(A*exp(-a*z)-B*exp(-b*z) +A*exp(-a*(L-z))-B*exp(-b*(L-z)))*100/L
        #w1=(A*cos(a*z)+A*cos(a*(L-z))+b)*100/L
        #w1= A*cos(2*pi*(z)/L)+A*cos(2*pi/L*(L-(z)))+b
        #w1= A*cos(2*pi*(z-N/2*dz)/L)+A*cos(2*pi/L*(L-(z-N/2*dz)))+b
        #w1= A*cos(2*pi*(z)/L)+A*cos(2*pi/L*(L-(z)))+b
        for i=1:N
        #w1[i]= -A*cos(2*pi*(i-1)*dz/L)+b
        #w1[i]=(A*exp(-a*(i-1)*dz)-exp(-b*(i-1)*dz) + A*exp(-a*(L-(i-1)*dz)) - exp(-b*(L-(i-1)*dz))*100/L
        w1[i]=(A*exp(-a*(i-1)*dz*(i-1)*dz)-exp(-b*(i-1)*dz*(i-1)*dz) + A*exp(-a*(L-(i-1)*dz)*(L-(i-1)*dz)) - exp(-b*(L-(i-1)*dz)*(L-(i-1)*dz)))*100/L
        end
        for i=1:N
                w[i,:]=circshift(w1,i-1)
        end

        return w
end

#
#function naive_mf(N,L,I0,u,scale,shift,total,h,A,B,a,b)
function naive_mf(N,L,I0,u,scale,shift,total,h,A,a,b)
# Euler solver for theta model naive mean field theory
# u_i = (1/N\pi)*\sum w_{ij} gain(I+u_j)

      	beta=1

#        h=.001
        dz=L/N
        z=[0:N-1]*dz

        if u==0
        u=zeros(N)
#        u=20*exp(-(10/L)^2*(z-L/2).^2)
        u=sin(2*pi/L*z)
        end

	# compute weights
#        A = 2.6
#        a = 20.0/L
#	b = 10.0/L
        a/=L
        b/=L
        w = shift/L .+ scale* weights(N,L,A,a,b)
#        w = scale*weights(N,L,A,a,b)

        # set input
#        if I==0
#                I = I0*exp(-(1/sigma*L)^2*(z-L/2).^2)
                I=I0
#        end

        for  i = 1:total

                f = gain.(I0+u)

               u+=-h*u+h/(pi)*w*f*dz


        end



        return u,w,I0,z,N
end
