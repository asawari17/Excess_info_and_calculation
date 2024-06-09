mutable struct Matrices
    E1::Float64
    E2::Float64
    E3::Float64
    E4::Float64
    B12::Float64
    B13::Float64
    B14::Float64
    B23::Float64
    B24::Float64
    B34::Float64
    R1_cl::Vector{Array{Float64, 2}}
    R1b_cl::Array{Float64, 2}
    R2a_cl::Vector{Array{Float64, 2}}
    R2ab_cl::Array{Float64, 2}
    R2b_cl::Vector{Array{Float64, 2}}
    R2bb_cl::Array{Float64, 2}
    R2c_cl::Vector{Array{Float64, 2}}
    R2cb_cl::Array{Float64, 2}
    R3a_cl::Vector{Array{Float64, 2}}
    R3ab_cl::Array{Float64, 2}
    R3b_cl::Vector{Array{Float64, 2}}
    R3bb_cl::Array{Float64, 2}
    R3c_cl::Vector{Array{Float64, 2}}
    R3cb_cl::Array{Float64, 2}
    R1_ch::Vector{Array{Float64, 2}}
    R1b_ch::Array{Float64, 2}
    R2a_ch::Vector{Array{Float64, 2}}
    R2ab_ch::Array{Float64, 2}
    R2b_ch::Vector{Array{Float64, 2}}
    R2bb_ch::Array{Float64, 2}
    R2c_ch::Vector{Array{Float64, 2}}
    R2cb_ch::Array{Float64, 2}
    R3a_ch::Vector{Array{Float64, 2}}
    R3ab_ch::Array{Float64, 2}
    R3b_ch::Vector{Array{Float64, 2}}
    R3bb_ch::Array{Float64, 2}
    R3c_ch::Vector{Array{Float64, 2}}
    R3cb_ch::Array{Float64, 2}
end

function Make_matrices(cl::Float64, ch::Float64)
    E1 = rand()
    E2 = rand()
    E3 = rand()
    E4 = rand()
    B12 = rand()
    B13 = rand()
    B14 = rand()
    B23 = rand()
    B24 = rand()
    B34 = rand()

    r12 = exp(-(B12 - E2) / beta)
    r13 = exp(-(B13 - E3) / beta)
    r14 = exp(-(B14 - E4) / beta)
    r23 = exp(-(B23 - E3) / beta)
    r24 = exp(-(B24 - E4) / beta)
    r34 = exp(-(B34 - E4) / beta)
    r21 = exp(-(B12 - E1) / beta)
    r31 = exp(-(B13 - E1) / beta)
    r41 = exp(-(B14 - E1) / beta)
    r32 = exp(-(B23 - E2) / beta)
    r42 = exp(-(B24 - E2) / beta)
    r43 = exp(-(B34 - E3) / beta)
    
    R = zeros(Float64, 4, 4)
    R[1, 2] = r12
    R[2, 1] = r21
    
    R[1, 3] = r13
    R[3, 1] = r31
    
    R[1, 4] = r14
    R[4, 1] = r41
    
    R[2, 3] = r23
    R[3, 2] = r32
    
    R[2, 4] = r24
    R[4, 2] = r42
    
    R[3, 4] = r34
    R[4, 3] = r43

    # Initialize R1, R2a, R2b, R2c, R3a, R3b as empty vectors
    R1_cl = Array{Array{Float64, 2}}(undef, 0)
    R2a_cl = Array{Array{Float64, 2}}(undef, 0)
    R2b_cl = Array{Array{Float64, 2}}(undef, 0)
    R2c_cl = Array{Array{Float64, 2}}(undef, 0)
    R3a_cl = Array{Array{Float64, 2}}(undef, 0)
    R3b_cl = Array{Array{Float64, 2}}(undef, 0)

    R1_ch = Array{Array{Float64, 2}}(undef, 0)
    R2a_ch = Array{Array{Float64, 2}}(undef, 0)
    R2b_ch = Array{Array{Float64, 2}}(undef, 0)
    R2c_ch = Array{Array{Float64, 2}}(undef, 0)
    R3a_ch = Array{Array{Float64, 2}}(undef, 0)
    R3b_ch = Array{Array{Float64, 2}}(undef, 0)
    

    R1_cl, R1b_cl = matrix_R1(R, cl)
    R2a_cl, R2ab_cl = matrix_R2a(R, cl)
    R2b_cl, R2bb_cl = matrix_R2b(R, cl)
    R2c_cl, R2cb_cl = matrix_R2c(R, cl)
    R3a_cl, R3ab_cl = matrix_R3a(R, cl)
    R3b_cl, R3bb_cl = matrix_R3b(R, cl)
    R3c_cl, R3cb_cl = matrix_R3c(R, cl)

    R1_ch, R1b_ch = matrix_R1(R, ch)
    R2a_ch, R2ab_ch = matrix_R2a(R, ch)
    R2b_ch, R2bb_ch = matrix_R2b(R, ch)
    R2c_ch, R2cb_ch = matrix_R2c(R, ch)
    R3a_ch, R3ab_ch = matrix_R3a(R, ch)
    R3b_ch, R3bb_ch = matrix_R3b(R, ch)
    R3c_ch, R3cb_ch = matrix_R3c(R, ch)

    return Matrices(E1, E2, E3, E4, B12, B13, B14, B23, B24, B34, R1_cl, R1b_cl, R2a_cl, 
    R2ab_cl, R2b_cl, R2bb_cl, R2c_cl, R2cb_cl, R3a_cl, R3ab_cl, R3b_cl, R3bb_cl, R3c_cl, 
    R3cb_cl, R1_ch, R1b_ch, R2a_ch, R2ab_ch, R2b_ch, R2bb_ch, R2c_ch, R2cb_ch, R3a_ch, 
    R3ab_ch, R3b_ch, R3bb_ch, R3c_ch, R3cb_ch)
end

function matrix_R1(R::Array{Float64, 2}, c::Float64)
    R_new = copy(R)
    R1b = zeros(Float64, 4, 4)
    R1b[3, 1] = c * R_new[3, 1]

    R_new[3, 1] = c * R_new[3, 1]
   
    R1_1 = copy(R_new)
    R1_1[2, 3] = R1_1[3, 2] = R1_1[3, 4] = R1_1[4, 3] = 0

    R1_2 = copy(R_new)
    R1_2[1, 2] = R1_2[2, 1] = R1_2[1, 4] = R1_2[4, 1] = 0

    R1_3 = copy(R_new)
    R1_3[1, 4] = R1_3[4, 1] = R1_3[2, 3] = R1_3[3, 2] = R1_3[3, 4] = R1_3[4, 3] = 0

    R1_4 = copy(R_new)
    R1_4[2, 3] = R1_4[3, 2] = R1_4[2, 4] = R1_4[4, 2] = R1_4[3, 4] = R1_4[4, 3] = 0

    R1_5 = copy(R_new)
    R1_5[1, 2] = R1_5[2, 1] = R1_5[2, 3] = R1_5[3, 2] = R1_5[3, 4] = R1_5[4, 3] = 0

    R1_6 = copy(R_new)
    R1_6[1, 2] = R1_6[2, 1] = R1_6[1, 4] = R1_6[4, 1] = R1_6[2, 4] = R1_6[4, 2] = 0
    
    R1_7 = copy(R_new)
    R1_7[1, 2] = R1_7[2, 1] = R1_7[2, 3] = R1_7[3, 2] = R1_7[1, 4] = R1_7[4, 1] = 0

    R1_8 = copy(R_new)
    R1_8[1, 2] = R1_8[2, 1] = R1_8[3, 4] = R1_8[4, 3] = R1_8[1, 4] = R1_8[4, 1] = 0
    
    R1_9 = copy(R_new)
    R1_9[2, 3] = R1_9[3, 2] = R1_9[2, 4] = R1_9[4, 2] = R1_9[1, 4] = R1_9[4, 1] = 0

    R1 = [R1_1, R1_2, R1_3, R1_4, R1_5, R1_6, R1_7, R1_8, R1_9]
    for i = 1:9 
        for j = 1:4
            R1[i][j, j] = -sum(R1[i][:, j])  # set the diagonal element to the negative sum of the other elements in the column
        end
    end
    return R1, R1b
end

function matrix_R2a(R::Array{Float64, 2}, c::Float64)
    R_new = copy(R)
    R2ab = zeros(Float64, 4, 4)
    R2ab[1, 3] = c * R_new[1, 3]
    R2ab[2, 4] = c * R_new[2, 4]
    
    R_new[1, 3] = c * R_new[1, 3]
    R_new[2, 4] = c * R_new[2, 4]
   
    R2a_1 = copy(R_new)
    R2a_1[2, 3] = R2a_1[3, 2] = R2a_1[1, 4] = R2a_1[4, 1] = 0

    R2a_2 = copy(R_new)
    R2a_2[2, 3] = R2a_2[3, 2] = R2a_2[1, 4] = R2a_2[4, 1]  = R2a_2[3, 4] = R2a_2[4, 3] = 0

    R2a_3 = copy(R_new)
    R2a_3[1, 2] = R2a_3[2, 1] = R2a_3[2, 3] = R2a_3[3, 2] = R2a_3[1, 4] = R2a_3[4, 1] = 0
    
    R2a = [R2a_1, R2a_2, R2a_3]
    for i = 1:3
        for j = 1:4
            R2a[i][j, j] = -sum(R2a[i][:, j])  # set the diagonal element to the negative sum of the other elements in the column
        end
    end
    return R2a, R2ab
end

function matrix_R2b(R::Array{Float64, 2}, c::Float64)
    R_new = copy(R)
    R2bb = zeros(Float64, 4, 4)
    R2bb[1, 3] = c * R_new[1, 3]
    R2bb[2, 1] = c * R_new[2, 1]

    R_new[1, 3] = c * R_new[1, 3]
    R_new[2, 1] = c * R_new[2, 1]

    R2b_1 = copy(R_new)
    R2b_1[2, 3] = R2b_1[3, 2] = R2b_1[1, 4] = R2b_1[4, 1]  = R2b_1[3, 4] = R2b_1[4, 3] = 0

    R2b_2 = copy(R_new)
    R2b_2[2, 3] = R2b_2[3, 2] = R2b_2[1, 4] = R2b_2[4, 1]  = R2b_2[2, 4] = R2b_2[4, 2] = 0

    R2b_3 = copy(R_new)
    R2b_3[4, 2] = R2b_3[2, 4] = R2b_3[2, 3] = R2b_3[3, 2] = R2b_3[3, 4] = R2b_3[4, 3] = 0
  
    R2b = [R2b_1, R2b_2, R2b_3]
    for i = 1:3
        for j = 1:4
            R2b[i][j, j] = -sum(R2b[i][:, j])  # set the diagonal element to the negative sum of the other elements in the column
        end
    end
    return R2b, R2bb
end

function matrix_R2c(R::Array{Float64, 2}, c::Float64)
    R_new = copy(R)
    R2cb = zeros(Float64, 4, 4)
    R2cb[3, 1] = c * R_new[3, 1]
    R2cb[2, 1] = c * R_new[2, 1]
    
    R_new[3, 1] = c * R_new[3, 1]
    R_new[2, 1] = c * R_new[2, 1]

    R2c_1 = copy(R_new)
    R2c_1[2, 3] = R2c_1[3, 2] = R2c_1[1, 4] = R2c_1[4, 1] = 0

    R2c_2 = copy(R_new)
    R2c_2[2, 4] = R2c_2[4, 2] = R2c_2[1, 4] = R2c_2[4, 1]  = 0

    R2c_3 = copy(R_new)
    R2c_3[4, 1] = R2c_3[1, 4] = R2c_3[3, 4] = R2c_3[4, 3] = 0
    
    R2c_4 = copy(R_new)
    R2c_4[4, 1] = R2c_4[1, 4] = R2c_4[2, 4] = R2c_4[4, 2] = R2c_4[2, 3] = R2c_4[3, 2]= 0
    
    R2c_5 = copy(R_new)
    R2c_5[4, 1] = R2c_5[1, 4] = R2c_5[3, 4] = R2c_5[4, 3] = R2c_5[2, 3] = R2c_5[3, 2]= 0
    
    R2c_6 = copy(R_new)
    R2c_6[4, 3] = R2c_6[3, 4] = R2c_6[2, 4] = R2c_6[4, 2] = R2c_6[2, 3] = R2c_6[3, 2]= 0
    
    R2c = [R2c_1, R2c_2, R2c_3, R2c_4, R2c_5, R2c_6]
    for i = 1:6
        for j = 1:4
            R2c[i][j, j] = -sum(R2c[i][:, j])  # set the diagonal element to the negative sum of the other elements in the column
        end
    end
    return R2c, R2cb
end

function matrix_R3a(R::Array{Float64, 2}, c::Float64)
    R_new = copy(R)
    R3ab = zeros(Float64, 4, 4)
    R3ab[3, 1] = c * R_new[1, 3]
    R3ab[2, 1] = c * R_new[2, 1]
    R3ab[4, 1] = c * R_new[4, 1]

    R_new[3, 1] = c * R_new[3, 1]
    R_new[2, 1] = c * R_new[2, 1]
    R_new[4, 1] = c * R_new[4, 1]
   
    R3a_1 = copy(R_new)
    R3a_1[2, 3] = R3a_1[3, 2] = R3a_1[3, 4] = R3a_1[4, 3] = 0

    R3a_2 = copy(R_new)
    R3a_2[2, 4] = R3a_2[4, 2] = R3a_2[2, 3] = R3a_2[3, 2]  = 0

    R3a_3 = copy(R_new)
    R3a_3[2, 4] = R3a_3[4, 2] = R3a_3[3, 4] = R3a_3[4, 3] = 0
    
    R3a_4 = copy(R_new)
    R3a_4[2, 3] = R3a_4[3, 2]= 0
    
    R3a_5 = copy(R_new)
    R3a_5[3, 4] = R3a_5[4, 3] = 0
    
    R3a_6 = copy(R_new)
    R3a_6[2, 4] = R3a_6[4, 2] = 0

    R3a_7 = copy(R_new)
    R3a_7[2, 4] = R3a_7[4, 2] = R3a_7[4, 3] = R3a_7[3, 4] = R3a_7[3, 2] = R3a_7[2, 3] = 0
    
    R3a = [R3a_1, R3a_2, R3a_3, R3a_4, R3a_5, R3a_6, R3a_7]

    for i = 1:7
        for j = 1:4
            R3a[i][j, j] = -sum(R3a[i][:, j])  # set the diagonal element to the negative sum of the other elements in the column
        end
    end
    return R3a, R3ab
end

function matrix_R3b(R::Array{Float64, 2}, c::Float64)
    R_new = copy(R)
    R3bb = zeros(Float64, 4, 4)
    R3bb[1, 3] = c * R_new[1, 3]
    R3bb[1, 2] = c * R_new[1, 2]
    R3bb[1, 4] = c * R_new[1, 4]

    R_new[1, 3] = c * R_new[1, 3]
    R_new[1, 2] = c * R_new[1, 2]
    R_new[1, 4] = c * R_new[1, 4]

    R3b_1 = copy(R_new)
    R3b_1[2, 3] = R3b_1[3, 2] = R3b_1[2, 4] = R3b_1[4, 2] = 0

    R3b_2 = copy(R_new)
    R3b_2[3, 4] = R3b_2[4, 3] = R3b_2[2, 3] = R3b_2[3, 2]  = 0

    R3b_3 = copy(R_new)
    R3b_3[2, 4] = R3b_3[4, 2] = R3b_3[3, 4] = R3b_3[4, 3] = 0
    
    R3b_4 = copy(R_new)
    R3b_4[2, 3] = R3b_4[3, 2]= 0
    
    R3b_5 = copy(R_new)
    R3b_5[3, 4] = R3b_5[4, 3] = 0
    
    R3b_6 = copy(R_new)
    R3b_6[2, 4] = R3b_6[4, 2] = 0

    R3b_7 = copy(R_new)
    R3b_7[2, 4] = R3b_7[4, 2] = R3b_7[4, 3] = R3b_7[3, 4] = R3b_7[2, 3] = R3b_7[3, 2] = 0
    
    R3b = [R3b_1, R3b_2, R3b_3, R3b_4, R3b_5, R3b_6, R3b_7]

    for i = 1:7
        for j = 1:4
            R3b[i][j, j] = -sum(R3b[i][:, j])  # set the diagonal element to the negative sum of the other elements in the column
        end
    end
    return R3b, R3bb
end

function matrix_R3c(R::Array{Float64, 2}, c::Float64)
    R_new = copy(R)
    R3cb = zeros(Float64, 4, 4)
    R3cb[1, 3] = c * R_new[1, 3]
    R3cb[2, 1] = c * R_new[2, 1]
    R3cb[4, 2] = c * R_new[4, 2]

    R_new[1, 3] = c * R_new[1, 3]
    R_new[2, 1] = c * R_new[2, 1]
    R_new[4, 2] = c * R_new[4, 2]

    
    R_new[1, 4] = R_new[4, 1] = R_new[2, 3] = R_new[3, 2] = R_new[3, 4] = R_new[4, 3] = 0

    R3c = [R_new]
    for j = 1:4
        R3c[1][j, j] = -sum(R3c[1][:, j])  # set the diagonal element to the negative sum of the other elements in the column
    end
    return R3c, R3cb
end