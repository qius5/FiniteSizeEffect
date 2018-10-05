include("naive_mf_cosine_fixI_Siwei.jl")
include("propagator_sub.jl")
N=200
L=1
u=reshape(readdlm("micro_sim_u_ave.dat"),N)
I0=reshape(readdlm("I0_negative.dat"), N)
#function naive_mf(N,L,I0,u,scale,shift,total,h,A,a,b)
u,w,I,z,N=naive_mf(N,L,I0,u,1,0,10000,0.01,1.5,30,20);
writedlm("u_N200_bump_test.dat", u)
writedlm("w_N200_bump.dat", w)
@time Cout, Dvp=propagator_sub(u,w,I,1,0.0001,1, 15, L);
writedlm("Cout_Cff_beta_1_N200_bump.dat", Cout)
