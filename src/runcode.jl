push!(LOAD_PATH, "/Users/siweiqiu/FiniteSizeEffect")
import FiniteSizeEffect
N=200
L=1
u=reshape(FiniteSizeEffect.readdlm("micro_sim_u_ave.dat"),N)
I0=reshape(FiniteSizeEffect.readdlm("I0_negative.dat"), N)
#function naive_mf(N,L,I0,u,scale,shift,total,h,A,a,b)
u,w,I1,z,N=FiniteSizeEffect.naive_mf(N,L,I0,u,1,0,10000,0.01,1.5,30,20);
FiniteSizeEffect.writedlm("u_N200_bump_test.dat", u)
FiniteSizeEffect.writedlm("w_N200_bump.dat", w)
@time Cout, Dvp=FiniteSizeEffect.propagator_sub(u,w,I1,1,0.001,1, 0.1, L);
FiniteSizeEffect.writedlm("Cout_Cff_beta_1_N200_bump.dat", Cout)
