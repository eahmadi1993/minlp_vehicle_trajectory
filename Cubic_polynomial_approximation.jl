module Polynomial


    function Cubic_polynomial_approximation(i_p, f_p,i_v_prime,f_v_prime,na,k)

        tau = range(0, stop=1, length=na)';
        n= 2;      #n is the dimention, here is 2D
        X = [];
        Y = [];
        for i=1:k
            for j =1:n
                a = 2*(i_p[i,j] - f_p[i,j]) + f_v_prime[i,j] + i_v_prime[i,j];
                b = 3*(f_p[i,j] - i_p[i,j]) - f_v_prime[i,j] - 2*i_v_prime[i,j];
                c = i_v_prime[i,j];
                d = i_p[i,j];
                if j==1
                    x = a*tau.^3 .+ b*tau.^2 .+ c*tau .+ d;
                    push!(X,x)
                else
                    y = a*tau.^3 .+ b*tau.^2 .+ c*tau .+ d;
                    push!(Y,y)
                end
            end
        end
        # println(size(Y))
        # println(size(tau))
        return tau, X, Y
    end
end
