module F_VectorGeneration

    struct F_Vectors
        FX
        FXd
        FXdd
        FY
        FYd
        FYdd
    end

    include("./F_formula.jl")
    import .FixedVector

    function Generate_F_vectors(i_p, f_p, i_v_prime, f_v_prime,Coeff_Matrix,k)

        CPiTau = Coeff_Matrix.CPiTau;
        SPiTau = Coeff_Matrix.SPiTau;
        C2PiTau = Coeff_Matrix.C2PiTau;
        S2PiTau = Coeff_Matrix.S2PiTau;
        Pi2 = Coeff_Matrix.Pi2;

        # each of its cell (FX) is correspond to F_x of each vehicle,...
        # e.g. for two vehicles: [F_x1; F_x2].
        FX = []; FXd = []; FXdd = [];
        for i = 1:k
            F_s,F_sdot,F_sddot = FixedVector.F_formula(i_p[i,1], f_p[i,1],
            i_v_prime[i,1], f_v_prime[i,1], CPiTau,SPiTau,C2PiTau,S2PiTau,Pi2)
            push!(FX,F_s)
            push!(FXd,F_sdot)
            push!(FXdd,F_sddot)
        end

        FY = []; FYd = []; FYdd = [];
        for i = 1:k
            F_s,F_sdot,F_sddot = FixedVector.F_formula(i_p[i,2], f_p[i,2],
             i_v_prime[i,2],f_v_prime[i,2], CPiTau, SPiTau,C2PiTau,S2PiTau,Pi2)
            push!(FY,F_s)
            push!(FYd,F_sdot)
            push!(FYdd,F_sddot)
        end
        # println(size(FY[1]))
        return F_Vectors(FX,FXd,FXdd,FY,FYd,FYdd)
    end

end
