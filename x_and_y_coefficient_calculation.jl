module Coeff_calculation

    # include("./lgwt.jl")
    # import .tau_mod

    function x_and_y_coefficient_calculation(ta_app,x_app, y_app, Nx, Ny, i_p, f_p,i_v_prime, f_v_prime)

        x_count = 2*Nx - 3;    # total number of unknown coefficients in x approximation
        y_count = 2*Ny - 3;    # total number of unknown coefficients in y approximation

        # ----------- x coefficient approximation -------------------
        # tau_x = tau_mod.lgwt(3*x_count+2, 1, 0);
        # tau_x = [0;tau_x;1];
        # RHS_x = interpolate(ta_app, x_app, tau_x, method="cubic")
        tau_x = ta_app';
        RHS_x = x_app';

        PiTauR = pi*tau_x;
        SPiTauR = sin.(PiTauR);
        CPiTauR = cos.(PiTauR);
        S2PiTauR = sin.(2*PiTauR);
        C2PiTauR = cos.(2*PiTauR);

        A_x = 0.5*(1 .- C2PiTauR);

        for i=3:2:Nx-1
            A_x = [A_x   cos.(i*PiTauR)-CPiTauR    sin.(i*PiTauR) - i*SPiTauR    cos.((i+1)*PiTauR) - C2PiTauR    sin.((i+1)*PiTauR) - (i+1)/2*S2PiTauR];
        end
        F_x =  (i_p[1] - f_p[1])./2*CPiTauR .+ 1/(2*pi)*(i_v_prime[1]-f_v_prime[1])*SPiTauR .+
        (i_p[1] + f_p[1])./2*C2PiTauR .+ 1/(4*pi)*(i_v_prime[1] + f_v_prime[1])*S2PiTauR;

        xx_coeffs = A_x\(RHS_x.-F_x);
        # ----------- y coefficient approximation -------------------
        # tau_y = lgwt(3*y_count+2,1,0);
        # tau_y = [0;tau_y;1];
        # RHS_y = interp1(ta_app, y_app, tau_y, 'spline');
        tau_y = ta_app';
        RHS_y = y_app';

        PiTaut = pi*tau_y;
        SPiTaut = sin.(PiTaut);
        CPiTaut = cos.(PiTaut);
        S2PiTaut = sin.(2*PiTaut);
        C2PiTaut = cos.(2*PiTaut);

        A_y = 0.5*(1 .- C2PiTaut);
        for i=3:2:Ny-1
            A_y = [A_y   cos.(i*PiTaut) - CPiTaut    sin.(i*PiTaut) - i*SPiTaut     cos.((i+1)*PiTaut) - C2PiTaut     sin.((i+1)*PiTaut) - (i+1)/2*S2PiTaut];
        end
        F_y = (i_p[2] - f_p[2])./2*CPiTaut .+ 1/(2*pi)*(i_v_prime[2]-f_v_prime[2])*SPiTaut .+
        (i_p[2] + f_p[2])./2*C2PiTaut .+ 1/(4*pi)*(i_v_prime[2] + f_v_prime[2])*S2PiTaut;
        yy_coeffs = A_y\(RHS_y.-F_y);

        # println(size(xx_coeffs))
        return xx_coeffs, yy_coeffs

    end # function
end  # module
