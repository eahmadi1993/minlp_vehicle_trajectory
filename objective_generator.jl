module CostFunc
# This module produces the objective function (\Delta_v).

    using Trapz

    function objective_generator(x0,A_x,A_xdot,A_xddot,A_y,A_ydot,A_yddot,FX,
        FXd,FXdd,FY,FYd,FYdd,nuNx,nuNy,T,T2,tau_points,k)

            #----------------- * Decision variables *-------------------------
            nuNxpNy = nuNx+nuNy;  # the total number of decision variables

            coeffs_x = []; coeffs_y = [];
            for i = 1:k
                Co_x = x0[1+(i-1)*nuNxpNy : nuNx+(i-1)*nuNxpNy];
                Co_y = x0[(i-1)*nuNxpNy+nuNx+1 : i*nuNxpNy];
                push!(coeffs_x,Co_x)
                push!(coeffs_y,Co_y)
            end

            # ------------------*** The compact form of FFF ***-------------

            xddot = [];
            yddot = [];
            for i=1:k
                X_ddot = (A_xddot * coeffs_x[i] + FXdd[i])/T2;
                Y_ddot = (A_yddot*coeffs_y[i] + FYdd[i])/T2;
                push!(xddot,X_ddot)
                push!(yddot,Y_ddot)
            end

            # -------------*** EOM (based on Huo's paper,2020) ***-------------
            O = [];
            for i=1:k
                fx = xddot[i];
                fy = yddot[i];
                f = sqrt.(fx.^2 + fy.^2);
                o = trapz(tau_points',f')
                push!(O,o)
            end

            Obj = sum(O); # minimize DeltaV
            # println(Obj)  # the value of objective_function
            println("The objective function is calculated successfully.")
            return Obj

    end  # function



end  # modul CostFunc
