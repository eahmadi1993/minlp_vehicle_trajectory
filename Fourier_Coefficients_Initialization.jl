module Initialization

    include("./Cubic_polynomial_approximation.jl")
    import .Polynomial
    include("./x_and_y_coefficient_calculation.jl")
    import .Coeff_calculation


    function Fourier_Coefficients_Initialization(i_p,f_p,i_v_prime, f_v_prime,na, Nx, Ny, a_max, k)

        ta_app,x_app,y_app = Polynomial.Cubic_polynomial_approximation(i_p, f_p,i_v_prime,f_v_prime,na,k);
        # println(size(x_app[1]))

        x_coeffs = [];  # it has k cell, e.g. for two vehicles: [x_coeffs;x2_coeffs]
        y_coeffs = [];
        for i=1:k
            xx_coeffs, yy_coeffs = Coeff_calculation.x_and_y_coefficient_calculation(ta_app,
            x_app[i], y_app[i], Nx, Ny, i_p[i,:], f_p[i,:], i_v_prime[i,:], f_v_prime[i,:]);
            push!(x_coeffs,xx_coeffs)
            push!(y_coeffs,yy_coeffs)
        end

        # -------- Initialization for time of flight (Tapp) ---------------
        # This is for the case that "T" is as a decision vaiable

        # The initialization of flight time T_App is estimated under the
        # following conditions: randomly selecting a vwhicle to drive
        # directly from the starting point to the end point with maximum
        # thrust acceleration and both start and end speeds are 0.

        # *I think this initialization for T_App is not suitable for traffic*

        S = sqrt( (i_p[1,1] - f_p[1,1])^2 + (i_p[1,2] - f_p[1,2])^2);
        T_app = sqrt(2*S/a_max);
        println("The initial value fot T is $T_app.")

        # println(size(x_coeffs))
        println("The function Fourier_Coefficients_Initialization, properly produces the initial value of x_coeff and y_coeff for all the vehicles.")

        return x_coeffs, y_coeffs,T_app
    end

end
