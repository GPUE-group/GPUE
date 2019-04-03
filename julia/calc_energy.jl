using DelimitedFiles
using FFTW

function derive(data, dim, dx)
    out = Array{Complex,2}(undef,size(data,1), size(data,2))
    dim2 = dim+1
    if dim2 == 3
        dim2 = 1
    end
    for i = 1:size(data,dim)
        for j = 1:size(data,dim2)
            #println(j,'\t',i)
            if (dim == 1)
                if (i == size(data,dim))
                    out[i,j] = (data[1, j] - data[i, j])/dx
                else
                    out[i,j] = (data[i + 1, j] - data[i, j])/dx
                end
            else
                if (i == size(data,dim))
                    out[j,i] = (data[j,1] - data[j,i])/dx
                else 
                    out[j,i] = (data[j,i + 1] - data[j,i])/dx
                end
            end
        end
    end
    return out
end

# We are calculating the energy to check <Psi|H|Psi>
function calculate_energy(wfc, H_k, H_r, Ax, Ay, g, xDim, yDim, dx, dy)
    hbar = 1.05457148e-34
    omega = 0.6
    omegaX = 6.283

    # Creating momentum and conjugate wavefunctions
    wfc_k = fft(wfc)
    wfc_c = conj(wfc)

    # Finding the momentum and real-space energy terms
    energy_k = wfc_c.*ifft((H_k) .* wfc_k)
    energy_r = wfc_c.* (H_r).* wfc
    energy_i = wfc_c.* (0.5*g*abs2.(wfc)).* wfc

    energy_l = wfc_c.*(im*hbar*(Ax.*derive(wfc,1,dx) + Ay.*derive(wfc,2,dy)))

    # Integrating over all space
    energy_if = 0
    energy_lf = 0
    energy_kf = 0
    energy_rf = 0
    for i = 1:xDim*yDim
        energy_if += real(energy_i[i])*dx*dy
        energy_rf += real(energy_r[i])*dx*dy
        energy_lf += real(energy_l[i])*dx*dy
        energy_kf += real(energy_k[i])*dx*dy
    end 

    println("Kinetic energy:", "\t\t\t", energy_kf)
    println("Potential energy:", "\t\t", energy_rf)
    println("Internal energy:", "\t\t", energy_if)
    println("Angular Momentum energy:", '\t', energy_lf)
    println("Total energy:", "\t\t\t", energy_kf+energy_if+energy_rf+energy_lf)

end

# Read in param.cfg file
function calculate(param_file::String, data_dir::String)
    parameters = Dict()

    for line in readlines(data_dir*param_file)
        if line != "[Params]"
            tmp = split(line,"=")
            parameters[tmp[1]] = tmp[2]
        end
    end

    xDim = parse(Int64, parameters["xDim"])
    yDim = parse(Int64, parameters["yDim"])
    dx = parse(Float64, parameters["dx"])
    dy = parse(Float64, parameters["dy"])

    omega = parse(Float64, parameters["omega"])
    omegaX = parse(Float64, parameters["omegaX"])

    g = parse(Float64, parameters["gDenConst"])

    Ax = readdlm(data_dir*"Ax_0")
    Ay = readdlm(data_dir*"Ay_0")
    K = readdlm(data_dir*"K_0")
    V = readdlm(data_dir*"V_0")

    Ax = reshape(Ax, xDim, yDim)
    Ay = reshape(Ay, xDim, yDim)
    K = reshape(K, xDim, yDim)
    V = reshape(V, xDim, yDim)

    start = 0
    ende = parse(Int64, parameters["esteps"])
    endg = parse(Int64, parameters["gsteps"])
    incr = parse(Int64, parameters["printSteps"])

    # Ground State Evolution
    println("Starting imaginary time energy calculation")
    for i = start:incr:endg
        wfc = readdlm(data_dir*"wfc_0_const_"*string(i)) +
              readdlm(data_dir*"wfc_0_consti_"*string(i))*im
        println(data_dir*"wfc_0_const_"*string(i), '\t',
                data_dir*"wfc_0_consti_"*string(i))
        wfc = reshape(wfc, xDim, yDim)
        calculate_energy(wfc, K, V, Ax, Ay, g, xDim, yDim, dx, dy)
        println()
    end

    println()

    # Ground State Evolution
    println("Starting real time energy calculation")
    for i = start:incr:ende
        wfc = readdlm(data_dir*"wfc_ev_"*string(i)) +
              readdlm(data_dir*"wfc_evi_"*string(i))*im
        println(data_dir*"wfc_0_const_"*string(i), '\t', 
                data_dir*"wfc_0_consti_"*string(i))
        wfc = reshape(wfc, xDim, yDim)
        calculate_energy(wfc, K, V, Ax, Ay, g, xDim, yDim, dx, dy)
        println()
    end


end

calculate("Params.dat", "../data/")

