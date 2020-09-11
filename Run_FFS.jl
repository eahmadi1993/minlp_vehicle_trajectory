module MyNLP

    include("objective_generator.jl")
    import .CostFunc

    using LinearAlgebra, JuMP, KNITRO
     # Ipopt, Cbc, Alpine

    function Run_FFS(x0,Coeff_Matrix,F_Vectors,T,T2,tau_points,k,a_max,ds)

        n = 2;   # The dimension is 2D.
        nuNx = Coeff_Matrix.nuNx;
        nuNy = Coeff_Matrix.nuNy;
        nuNxpNy = nuNx+nuNy;

        N_DP = length(tau_points);

        #  for x coordinate
        A_x      = Coeff_Matrix.A_x;
        A_xdot   = Coeff_Matrix.A_xdot;
        A_xddot  = Coeff_Matrix.A_xddot;

        FX = F_Vectors.FX;     # FX is 1*k and has k cell
        FXd = F_Vectors.FXd;
        FXdd = F_Vectors.FXdd;

        # for y coordinate
        A_y      = Coeff_Matrix.A_y;
        A_ydot   = Coeff_Matrix.A_ydot;
        A_yddot  = Coeff_Matrix.A_yddot;

        FY = F_Vectors.FY;
        FYd = F_Vectors.FYd;
        FYdd = F_Vectors.FYdd;

        # ------------ Solving the NLP using JuMP ------------------
        # setting for using Knitro
         model = JuMP.Model(solver=KNITRO.KnitroSolver(options_file="knitro.opt"))

        # setting for using Juniper
        # optimizer = Juniper.Optimizer
        # nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
        # model = Model(optimizer_with_attributes(optimizer, "nl_solver"=>nl_solver))

        # setting for using Ipopt
        # model = Model(Ipopt.Optimizer)

        # coeffs_x = zeros(nuNxpNy,k);
        # coeffs_y = [];
        # defining variables
        # for i=1:k
            coeffs_x = @variable(model, coeffs_x[1:nuNx, 1:k])
            coeffs_y = @variable(model, coeffs_y[1:nuNy, 1:k])
            @variable(model, δ, Bin)
            # @variable(model, X_y[1:nuNy])
            # push!(coeffs_x, X_x)
            # push!(coeffs_y, X_y)
        # end

        # The compact form of FFS
        x = [];        y = [];
        xdot =[];      ydot =[];
        xddot = [];    yddot = [];
        for i=1:k
            X      =  A_x    * coeffs_x[:,i] + FX[i];
            X_dot   = (A_xdot * coeffs_x[:,i] + FXd[i]) /T;
            X_ddot = (A_xddot * coeffs_x[:,i] + FXdd[i])/T2;
            Y      =  A_y    *coeffs_y[:,i] + FY[i];
            Y_dot   = (A_ydot *coeffs_y[:,i] + FYd[i]) /T;
            Y_ddot = (A_yddot*coeffs_y[:,i] + FYdd[i])/T2;
            push!(x,X)
            push!(xdot,X_dot)
            push!(xddot,X_ddot)
            push!(y,Y)
            push!(ydot,Y_dot)
            push!(yddot,Y_ddot)
        end
        println(size(x[1]))
        # ------- Objectivr function definition --------
        # Obj = CostFunc.objective_generator(x0,A_x,A_xdot,A_xddot,A_y,A_ydot,
        # A_yddot,FX,FXd,FXdd,FY,FYd,FYdd,nuNx,nuNy,T,T2,tau_points,k)

        # O = [];
        # for i=1:k
        #     fx = xddot[i];
        #     fy = yddot[i];
        #     f = sqrt.(fx.^2 + fy.^2);
        #     o = trapz(tau_points',f')
        #     push!(O,o)
        # end
        # println("OK?")
        #
        # @NLobjective(model, Min,  sum(O));   # minimize DeltaV
        # println("The objective function is calculated successfully.")

        # --------- Constraints definition ---------
        # Constraints on thrust acceleration
        A_total = zeros(1,n);
        # A_total = [[] []];
        for i=1:k
            fx = xddot[i];
            fy = yddot[i];
            Aa = [fx fy];
            A_total = [A_total;Aa];
        end
        A_total = A_total[2:end,:];
        # println(A_total[1,1])

        # here I should use norm or sqrt, but, JuMP don't accept them. So, I used a_max^2 instead of a_max
        for i = 1:k*N_DP
            @constraint(model, A_total[i,1]^2 + A_total[i,2]^2 <= a_max^2); # it works
        end

        # for i = 1:k*N_DP
        #     @NLconstraint(model, @NLexpression(model, sqrt(A_total[i,1]^2 + A_total[i,2]^2)) <= a_max);
        # end
        println("The acceleration constraints are defind.")


        # Constraints on safety
        d = [];
        for i=1:k
            dd = [x[i] y[i]];
            push!(d, dd)
        end
        # println(size(d[1]))

        for i = 1:k
            for   j = i + 1:k
                @constraint(model, (d[i][:,1] - d[j][:,1]).^2 + (d[i][:,2] - d[j][:,2]).^2  .<= ds^2);
            end
        end
        println("The safety constraints are defind.")


        # Constraints for intersection geometry modeling (Binary constraints)
        ux = 5;   uy = 1;   Hx = 5.25;     Hy = 1.25;
        lx = 1;   ly = -1;  hx = 0.75;     hy = -1.25;
        eps = 1e-5;

        for i = 1:k
            @constraint(model, x[i].-ux .<= (Hx-ux)*(1-δ))
            @constraint(model, x[i].-ux .>= eps + δ*(hx-ux-eps))
            @constraint(model, lx.-x[i] .<= (lx-hx)*(lx-δ))
            @constraint(model, lx.-x[i] .>= eps + δ*(lx-Hx-eps))
            @constraint(model, y[i].-uy .<= (Hy-uy)*(1-δ))
            @constraint(model, y[i].-uy .>= eps + δ*(hy-uy-eps))
            @constraint(model, ly.-y[i] .<= (ly-hy)*(ly-δ))
            @constraint(model, ly.-y[i] .>= eps + δ*(ly-Hy-eps))
        end
        println("The binary constraints are built.")

        @objective(model, Min,  0)
        # @time begin
        optimize!(model)
        # end
        # #
        println("The value of objective function is :  ", objective_value(model::Model))
        # println("The size of first set of decision variable (coeff_x) is: ", size(value.(coeffs_x)))
        println("Optimal solution of coeff_x is = \n", value.(coeffs_x))
        # # println("The value of first coulmn of coeff_x is = \n", value.(coeffs_x[:,1]))
        # println("Optimal solution of coeff_y is = \n", value.(coeffs_y))
        # # println(value.(coeffs_x))




        tt = 3;
        return tt
    end # funct

end  # module, MyNLP
