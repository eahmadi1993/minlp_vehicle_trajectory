module MatrixGeneration
    struct Coeff_Matrix
        A_x
        A_xdot
        A_xddot
        A_y
        A_ydot
        A_yddot;
        nuNx
        nuNy
        CPiTau
        SPiTau
        S2PiTau
        C2PiTau
        Pi2
    end

    function Coeff_Matrix_Generation(Nx,Ny,tau_points)

        nuNx = 2*Nx - 3;   # total number of unknown coefficients in x approximation
        nuNy = 2*Ny - 3;   # total number of unknown coefficients in y approximation

        PiTau = pi.*tau_points;
        CPiTau = cos.(PiTau)
        SPiTau = sin.(PiTau);
        S2PiTau = sin.(2*PiTau);
        C2PiTau = cos.(2*PiTau);
        Pi2 = pi^2;

        #-------------------------------- x ---------------------------------------
        A_x = 0.5.*(1 .- C2PiTau);
        A_xdot  = pi.*S2PiTau;
        A_xddot = 2*Pi2.*C2PiTau;
        for i=3:2:Nx-1
            il = i + 1;
            CiPiTau  = cos.(i.*PiTau);
            SiPiTau  = sin.(i.*PiTau);
            CilPiTau = cos.(il.*PiTau);
            SilPiTau = sin.(il.*PiTau);
             A_x     = [A_x   CiPiTau .- CPiTau  SiPiTau .- i.*SPiTau   CilPiTau .- C2PiTau    SilPiTau.-il./2*S2PiTau];
             A_xdot  = [A_xdot  pi.*(SPiTau.-i.*SiPiTau)   i*pi.*(CiPiTau.-CPiTau)  pi.*(2*S2PiTau.-il.*SilPiTau)   il*pi.*(CilPiTau.-C2PiTau)];
             A_xddot = [A_xddot  Pi2.*(CPiTau.-i^2*CiPiTau)  Pi2*i.*(SPiTau.-i.*SiPiTau)  Pi2.*(4*C2PiTau.-il^2*CilPiTau)  Pi2*il.*(2*S2PiTau.-il*SilPiTau)];
        end

        #-------------------------------- y ---------------------------------------
        A_y     = 0.5.*(1 .- C2PiTau);
        A_ydot  = pi.*S2PiTau;
        A_yddot = 2*Pi2.*C2PiTau;
        for i=3:2:Ny-1
            il = i+1;
            CiPiTau  = cos.(i.*PiTau);
            SiPiTau  = sin.(i.*PiTau);
            CilPiTau = cos.(il.*PiTau);
            SilPiTau = sin.(il.*PiTau);
             A_y     = [A_y   CiPiTau .- CPiTau   SiPiTau .- i.*SPiTau   CilPiTau .- C2PiTau   SilPiTau .- il./2*S2PiTau];
             A_ydot  = [A_ydot  pi.*(SPiTau.-i.*SiPiTau)   i*pi.*(CiPiTau.-CPiTau)  pi.*(2*S2PiTau.-il.*SilPiTau)   il*pi.*(CilPiTau.-C2PiTau)];
             A_yddot = [A_yddot  Pi2.*(CPiTau.-i^2*CiPiTau)  Pi2*i.*(SPiTau.-i.*SiPiTau)  Pi2.*(4*C2PiTau.-il^2*CilPiTau)  Pi2*il.*(2*S2PiTau.-il*SilPiTau)];
        end

        return Coeff_Matrix(A_x,A_xdot,A_xddot,A_y,A_ydot,A_yddot,nuNx,nuNy,CPiTau,SPiTau,S2PiTau,C2PiTau,Pi2)
    end

end
