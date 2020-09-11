module FixedVector

function F_formula(s1, s2, sdot1, sdot2,CPiTau, SPiTau, C2PiTau, S2PiTau, Pi2)
    C1 = s1 - s2;
    C2 = s1 + s2;
    C3 = ( sdot1 - sdot2 )/2;
    C4 = sdot1 + sdot2;
    C5 = C1*CPiTau;
    C6 = C2*C2PiTau;
    C7 = C3*SPiTau;
    C8 = C4*S2PiTau;
    F_s = 0.5*(C5.+C6) .+ 1/pi*(C7.+0.25*C8);
    F_sdot  = -pi/2*(C1*SPiTau.+2*C2*S2PiTau).+C3*CPiTau.+0.5*C4*C2PiTau;
    F_sddot = -0.5*Pi2*(C5.+4*C6) .- pi*(C7.+C8);
    return F_s, F_sdot, F_sddot
end

end  # moduleF
