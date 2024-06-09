"steady state distribution"
function steady_state(R::Array{Float64})
    V = nullspace(R)
    distn = V[:, 1] / sum(V[:, 1])
    return distn
end

function eigen_analysis(R::Array{Float64})
    ei = eigen(R)
    val = ei.values
    vecs = ei.vectors

    # Compute permutation to sort eigenvalues in descending order
    perm = sortperm(val, rev=true)

    # Sort eigenvalues and eigenvectors
    val_descending = val[perm]
    vecs_descending = vecs[:, perm]
    if findmax(val_descending)[2] != 1
        println("Error in eigen orders")
    end
    time_step = 1/(abs(last(val_descending)) * 500)
    time_inf = 40/(abs(val_descending[2]))
    return val_descending, vecs_descending
end


function coefficients(Pf::Vector{Float64}, v::Array{Float64}, Pi::Vector{Float64})
    vec = v[:, 2:end]  # Exclude the first column
    m = hcat(Pf, vec)
    coeff = inv(m) * Pi  # Solve the system of linear equations
    return coeff, m
end

"calculate alpha"
function update_alpha(c1::Float64, c2::Float64)
    if c2 < 0
        error("c2 must be positive.")
    end
    a = log(c1 / c2) + (c2 - c1) / c1
    return a
end

"calculate accumulated excess information at infinite time"

function aei_inf(R_real::Array{Float64}, R::Array{Float64}, coeff::Vector{Float64}, lambda::Vector{Float64}, v::Array{Float64}, force::Float64, P_ch::Vector{Float64}, P_cl::Vector{Float64}, P_tau::Vector{Float64}, c::Float64)
    accumulated_excess_info = 0.0  # Renaming the variable to avoid conflict with the function name
    log_cutoff = -4.6
    t_eq = []
    
    if norm(P_tau - v[:,1]) < 10^-18
        max_t_eq = -1.0
        return accumulated_excess_info, max_t_eq
    end

    for k = 2:ns
        num_rv = 0.0
        for i = 1:ns
            for j = 1:ns
                if i != j
                    g = R[i, j] * v[j, k]
                    num_rv = num_rv + g
                end
            end
        end
        h = abs(coeff[k] * v[1, k])
        num_rv = num_rv * coeff[k]
        accumulated_excess_info = accumulated_excess_info + num_rv * (-1/ lambda[k])
        t_eq_k = (log_cutoff - log(h)) / lambda[k]
        push!(t_eq, t_eq_k)
        if lambda[k] >= 0.0
            print("Lambda crisis!!")
            return 0.0 
        end
    end

    accumulated_excess_info = force * accumulated_excess_info
    
    max_t_eq = maximum(t_eq)
    if max_t_eq < 0.0
        #println(t_eq)
        max_t_eq = 0.0
        #accumulated_excess_info = 0.0
        # print("Negative equilibration time due to log_cutoff")
        # for k = 2:ns
        #     num_rv = 0.0
        #     for i = 1:ns
        #         for j = 1:ns
        #             if i != j
        #                 g = R[i, j] * v[j, k]
        #                 num_rv += g
        #             end
        #         end
        #     end
        #     num_rv = num_rv * coeff[k]
        #     println(log(abs(num_rv)))
            
        #     accumulated_excess_info += num_rv * (-1/ lambda[k])
        #     t_eq_k = (log_cutoff - log(abs(num_rv))) / lambda[k]
        #     push!(t_eq, t_eq_k)
        #     if lambda[k] >= 0.0
        #         print("Lambda crisis!!")
        #         return 0.0 
        #     end
        # end
        # println("P_tau-P_cl", P_tau-P_cl)
        # println("P_tau-P_ch", P_tau-P_ch)
        # println("concentration", c)
        # println("coeff", coeff)
        # println("lambda", lambda)
        # println("aei", accumulated_excess_info)
        # if c < 4.5
        #     vec_cl = P_ch
        #     for m = 2:ns
        #         vec_cl += coeff[m] * v[:,m]
        #     end
        # end

        # println("v_cl-P_cl", vec_cl-P_cl)

        # return 0.0 
    end
    
    return accumulated_excess_info, max_t_eq
end

#D_dot_NESS
function D_dot_NESS(R::Array{Float64}, P::Vector{Float64}, force::Float64)
    d_dot = 0.0
    for i = 1:ns
        for j = 1:ns
                d_dot += R[i,j] * P[j]
        end
    end
    d_dot = d_dot * force
    return d_dot
end


    #Calculating AEI(tau)
function aei(R::Array{Float64}, coeff::Vector{Float64}, lambda::Vector{Float64}, v::Array{Float64}, force::Float64, tau::Float64)
    aei = 0
    for i = 1:ns
        for j = 1:ns
            if i != j
                for k = 2:ns
                    g = R[i, j] * coeff[k] * v[j, k] * (1 / lambda[k]) * (exp(lambda[k] * tau) - 1)
                    aei += g
                end
            end
        end
    end
    aei = force * aei
    return aei
end