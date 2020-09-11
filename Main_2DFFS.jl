# Trajectory planning, finite Fourier series (FFS) idea
# T is a constant parameter
##### In this code, I just define the acceleration and safety constraints. Note that, I didn't use norm.
##### So, I should rewrite the constraint using norm.
##### I put objective_func=0. So, I have to implement numerical integration (trapz).

include("./lgwt.jl")
import .tau_mod
include("./Coeff_Matrix_Generation.jl")
import .MatrixGeneration
include("./Generate_F_vectors.jl")
import .F_VectorGeneration
include("./Fourier_Coefficients_Initialization.jl")
import .Initialization
include("./Run_FFS.jl")
import .MyNLP


using LinearAlgebra, Trapz, JuMP, KNITRO
# Juniper, Ipopt, Cbc, Alpine

T = 100;
T2 = T^2

## ---------------------- Initial parameters for vehicles ----------------------- ##

k = 2;       # the number of vehicles
ds = 1.6;     # safe distance between two vehicles
a_max = 0.5;  # max thrust acceleration (m/s^2)
m = 120;      # the number of Discritizations Point (DPs)
na = 600;     # the number of na DPs of the Cubic Polynomial (CP)

# Note: the number of Fourier terms should be even (and not odd) values
Nx = 6;       # the number of Fourier terms
Ny = 6;       # the number of Fourier terms

# the initial positions of Vehicles (m), each row corresponds to one vehicle
# i_p = [0 3
#       -4 0
#       -2 0
#        0 -2
#        0 1
#        5 0
#        0 -4
#        3 0
#        0 5
#        7 0] .+ 30;
i_p = [2.0 0.0 ;4.5 0.0];
f_p = [-1.25 0; 0.5 0];
#

# the final positions of Vehicles (m), each row corresponds to one vehicle
# f_p = [0 -3
#      4 0
#      9 0
#      0 2
#      0 -1
#      -5 0
#      0 4
#      -9 0
#       0 -5
#      -7 0] .+ 30;
#

# the initial and terminal velocitie of vehicles (m/s),
# each row corresponds to one vehicle
i_v = [0 0; 0 0; 0 0; 0 0; 0 0;0 0;0 0;0 0; 0 0;0 0];
f_v = [0 0; 0 0; 0 0; 0 0;0 0;0 0;0 0; 0 0; 0 0;0 0];

# ----- test --------------------
# This part is just to ensure the initial and final distances are greater than ds.

# Initial distance between each two vehicles
S1 = []
for i = 1:k
    global  S1
    for j = i+1:k
            s1 = sqrt( (i_p[i,1] - i_p[j,1])^2 + (i_p[i,2] - i_p[j,2])^2 )
            S1 = [S1;s1]
    end
end

# Final distance between each two vehicles
S2 = [];
for i = 1:k
    global  S2
    for j = i+1:k
        s2 = sqrt( (f_p[i,1] - f_p[j,1])^2 + (f_p[i,2] - f_p[j,2])^2);
        S2 = [S2;s2];
    end
end

println("The initial and final distances between each two vehicles are greater than the safe distance (ds).")

## ------------------------ Time scaling  (0<= \tau = t/T <= 1) ------------------------- ##
tau_points = tau_mod.lgwt(m, 1, 0)
tau_points = [0;tau_points;1];
N_DP = length(tau_points);
println("The number of discretization point is $N_DP.")
println("The discretization in tau_points are obtained and also are fixed.")

## ---------------------- Velocity coordinaes in the new \tau-unit ----------------------- ##
i_v_prime = T .* i_v;
f_v_prime = T .* f_v;

## --------------- Generation of The Required constant Matrices and Vectors -------------- ##
# once discretization in \tau and number of Fourier terms are fixed these
# matrices are fixed! of course these values are problem specific

Coeff_Matrix = MatrixGeneration.Coeff_Matrix_Generation(Nx, Ny, tau_points);
println("The coefficient matrix is successfully built.")

F_Vectors = F_VectorGeneration.Generate_F_vectors(i_p, f_p, i_v_prime,
            f_v_prime,Coeff_Matrix,k);
# println(size(F_Vectors.FY[1]))
println("The F_vectors in x_axis and y-axis are seccessfully created.")


## -------------- Initialization of the Fourier coefficients -------------- ##
x_coeffs,y_coeffs,T_app = Initialization.Fourier_Coefficients_Initialization(i_p,
     f_p,i_v_prime, f_v_prime,na, Nx, Ny, a_max, k);

println("The initialization for decision variables is done!")

## --------------------- Nested Constrained Fmincon---------------------- ##
x0 = [];
for i=1:k
    global x0 = [x0;x_coeffs[i];y_coeffs[i]];
end
println("Finally, the size of initialization vector is $(size(x0)).")

tt = MyNLP.Run_FFS(x0,Coeff_Matrix,F_Vectors,T,T2,tau_points,k,a_max,ds)
# x_sol, Delta_v
