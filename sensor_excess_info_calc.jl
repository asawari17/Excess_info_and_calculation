module distinguish
using ArgParse, DelimitedFiles, LinearAlgebra
export run_all
const ns = 4
const beta = 1.0


mutable struct Results
    aei_infinity::Float64
    aei_t_eq::Float64
    max_t_eq::Float64
    ati::Float64
    ati_t_eq::Float64
    speed_inf::Float64
    speed_t_eq::Float64
    d_dot::Float64
end

function run_all(cl::Float64, ch::Float64, index::Int)
    a = Make_matrices(cl, ch)
    tau_range = range(0.0, stop=6.0, step=0.01) |> collect
    backbone = ["1", "2a", "2b", "2c", "3a", "3b", "3c"]

    for field in backbone
        field_name1 = Symbol("R$field" * "_cl")  # Construct the field name
        matrices_cl = getproperty(a, field_name1)  # Access the field using getproperty
        field_name2 = Symbol("R$field" * "b_cl")
        Rb_cl = getproperty(a, field_name2)
        field_name3 = Symbol("R$field" * "_ch")  # Construct the field name
        matrices_ch = getproperty(a, field_name3)  # Access the field using getproperty
        field_name4 = Symbol("R$field" * "b_ch")
        Rb_ch = getproperty(a, field_name4)
    
        for k in 1:length(matrices_cl)
            R_cl = matrices_cl[k]
            R_ch = matrices_ch[k]
            println(index, field, k)
            # Construct the variable name and retrieve its value
            variable_name = Symbol("R$(field)_$(k)")
            
            # Create folder based on index range
            folder_index = 1 + floor(Int, (index - 1) / 100) * 100
            folder_name = "itr$folder_index"
            if !isdir(folder_name)
                mkdir(folder_name)
            end

            for i =  1:length(tau_range)
                #println(tau)
                if tau_range[i] != 0.0
                    b = run_binding(R_cl, Rb_cl, R_ch, Rb_ch, cl, ch, tau_range[i])
                else
                    b = run_binding_tau_zero(R_cl, Rb_cl, R_ch, Rb_ch, cl, ch, tau_range[i])
                end
                filename = joinpath(folder_name, "aei_inf_$(index)_R$(field)_$(k)_tau_$(tau_range[i]).txt")
                open(filename, "w") do file
                    write(file, "index: $index\t matrix: $(variable_name)\t E1: $(a.E1)\t E2: $(a.E2)\t E3: $(a.E3)\t E4: $(a.E4)\t B12: $(a.B12)\t B13: $(a.B13)\t B14: $(a.B14)\t B23: $(a.B23)\t B24: $(a.B24)\t B34: $(a.B34)\t tau_$(i): $(tau_range[i])\t aei_inf_$(i): $(b.aei_infinity)\t aei_t_eq_$(i): $(b.aei_t_eq)\t max_t_eq_$(i): $(b.max_t_eq)\t ati_$(i):$(b.ati)\t ati_t_eq_$(i):$(b.ati_t_eq)\t speed_inf_$(i):$(b.speed_inf)\t speed_t_eq_$(i):$(b.speed_t_eq)\t d_dot_$(i):$(b.d_dot)\n")
                end
            end
        end
    end
    return nothing
end


function run_binding(R_cl::Array{Float64}, Rb_cl::Array{Float64}, R_ch::Array{Float64}, Rb_ch::Array{Float64}, cl::Float64, ch::Float64, tau::Float64)
    P_ch = steady_state(R_ch)
    P_cl = steady_state(R_cl)
    D_ref = log(ns)
    #Aei for DAB
    lambda_ch, v_ch = eigen_analysis(R_ch)
    force_AB = update_alpha(ch, cl)
    eigenvalue = lambda_ch[2]

    #Aei for DBA
    lambda_cl, v_cl = eigen_analysis(R_cl)
    force_BA = update_alpha(cl, ch)

    # time_step = min(time_step_AB, time_step_BA)
    # time_inf = max(time_inf_AB, time_inf_BA)

    c_h_in_l, normed_v_cl = coefficients(P_cl, v_cl, P_ch) #coefficients for the decomposition of P_ch to the eigen basis of R_cl; to later obtain the probability distribution after time t_r of resetting 
    P_tau = P_cl
        for i = 2:ns
            P_tau += c_h_in_l[i] * exp(lambda_cl[i] * tau) * normed_v_cl[:,i]
        end
        c_tau_in_h, normed_v_ch = coefficients(P_ch, v_ch, P_tau)
        c_tau_in_l, normed_v_cl = coefficients(P_cl, v_cl, P_tau)
        aei_inf_AB, max_t_eq_Rch = aei_inf(R_ch,Rb_ch, c_tau_in_h, lambda_ch, normed_v_ch, force_AB, P_ch, P_cl, P_tau, ch)
        aei_inf_BA, max_t_eq_Rcl = aei_inf(R_cl,Rb_cl, c_tau_in_l, lambda_cl, normed_v_cl, force_BA, P_ch, P_cl, P_tau, cl)
        
        aei_infinity = aei_inf_AB + aei_inf_BA
        
        max_t_eq = max(max_t_eq_Rch, max_t_eq_Rcl)
        
        aei_t_eq_AB = aei(Rb_ch, c_tau_in_h, lambda_ch, normed_v_ch, force_AB, max_t_eq)
        aei_t_eq_BA = aei(Rb_cl, c_tau_in_l, lambda_cl, normed_v_cl, force_BA, max_t_eq)
        
        aei_t_eq = aei_t_eq_AB + aei_t_eq_BA
        
        d_dot_AB = D_dot_NESS(Rb_ch, P_ch, force_AB)
        d_dot_BA = D_dot_NESS(Rb_cl, P_cl, force_BA)
        ati = aei_infinity + (d_dot_AB + d_dot_BA) * max_t_eq

        ati_t_eq = aei_t_eq + (d_dot_AB + d_dot_BA) * max_t_eq

        if max_t_eq == 0.0
            speed_inf = d_dot_AB + d_dot_BA
            speed_t_eq = d_dot_AB + d_dot_BA
        else
            speed_inf = ati / max_t_eq
            speed_t_eq = ati_t_eq / max_t_eq
        end
        d_dot = d_dot_AB + d_dot_BA
    
    return Results(aei_infinity, aei_t_eq, max_t_eq, ati, ati_t_eq, speed_inf, speed_t_eq, d_dot)
end

function run_binding_tau_zero(R_cl::Array{Float64}, Rb_cl::Array{Float64}, R_ch::Array{Float64}, Rb_ch::Array{Float64}, cl::Float64, ch::Float64, tau::Float64)
    P_ch = steady_state(R_ch)
    P_cl = steady_state(R_cl)
    D_ref = log(ns)
    #Aei for DAB
    lambda_ch, v_ch = eigen_analysis(R_ch)
    force_AB = update_alpha(ch, cl)
    eigenvalue = lambda_ch[2]

    #Aei for DBA
    lambda_cl, v_cl = eigen_analysis(R_cl)
    force_BA = update_alpha(cl, ch)

    c_h_in_l, normed_v_cl = coefficients(P_cl, v_cl, P_ch)
    c_h_in_h, normed_v_ch = coefficients(P_ch, v_ch, P_ch)

    aei_inf_AB = 0 #by design
    aei_inf_BA, max_t_eq = aei_inf(R_cl, Rb_cl, c_h_in_l, lambda_cl, normed_v_cl, force_BA, P_ch, P_cl, P_ch, cl)
    aei_infinity = aei_inf_AB + aei_inf_BA

    aei_t_eq_AB = aei(Rb_ch, c_h_in_h, lambda_ch, normed_v_ch, force_AB, max_t_eq)
    aei_t_eq_BA = aei(Rb_cl, c_h_in_l, lambda_cl, normed_v_cl, force_BA, max_t_eq)
        
    aei_t_eq = aei_t_eq_AB + aei_t_eq_BA

    d_dot_AB = D_dot_NESS(Rb_ch, P_ch, force_AB)
    d_dot_BA = D_dot_NESS(Rb_cl, P_cl, force_BA)

    ati = aei_infinity + (d_dot_AB + d_dot_BA) * max_t_eq

    ati_t_eq = aei_t_eq + (d_dot_AB + d_dot_BA) * max_t_eq

    if max_t_eq == 0.0
        speed_inf = d_dot_AB + d_dot_BA
        speed_t_eq = d_dot_AB + d_dot_BA
    else
        speed_inf = ati / max_t_eq
        speed_t_eq = ati_t_eq / max_t_eq
    end
    d_dot = d_dot_AB + d_dot_BA

    return Results(aei_infinity, aei_t_eq, max_t_eq, ati, ati_t_eq, speed_inf, speed_t_eq, d_dot)
end