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

function calculate()
    xDim = 256
    yDim = 256
    dx = 1.72817e-06
    dy = 1.72817e-06
    omega = 0.6
    omegaX = 6.283
    data_dir = "../data/"

    Ax = readdlm(data_dir*"Ax_0")
    Ay = readdlm(data_dir*"Ay_0")
    K = readdlm(data_dir*"K_0")
    V = readdlm(data_dir*"V_0")

    #wfc = readdlm("wfc_load") + readdlm("wfci_load")*im
    wfc = readdlm("../data/wfc_0_const_50000") + readdlm("../data/wfc_0_consti_50000")*im
    println("../data/wfc_0_const_50000", '\t', "../data/wfc_0_consti_50000")

    Ax = reshape(Ax, xDim, yDim)
    Ay = reshape(Ay, xDim, yDim)
    #Ax /= omega*omegaX
    #Ay /= omega*omegaX
    K = reshape(K, xDim, yDim)
    V = reshape(V, xDim, yDim)
    wfc = reshape(wfc, xDim, yDim)
    g = 6.80274e-41

    calculate_energy(wfc, K, V, Ax, Ay, g, xDim, yDim, dx, dy)

end

function calculate_range(start::Int64, incr::Int64, final::Int64)

    xDim = 256
    yDim = 256
    dx = 1.72817e-06
    dy = 1.72817e-06
    omega = 0.6
    omegaX = 6.283
    data_dir = "../data/"
    Ax = readdlm(data_dir*"Ax_0")
    Ay = readdlm(data_dir*"Ay_0")
    K = readdlm(data_dir*"K_0")
    V = readdlm(data_dir*"V_0")

    Ax = reshape(Ax, xDim, yDim)
    Ay = reshape(Ay, xDim, yDim)
    K = reshape(K, xDim, yDim)
    V = reshape(V, xDim, yDim)
    g = 6.80274e-41

    for i = start:incr:final
        wfc = readdlm(data_dir*"wfc_0_const_"*string(i)) +
              readdlm(data_dir*"wfc_0_consti_"*string(i))*im
        println(data_dir*"wfc_0_const_"*string(i), '\t', data_dir*"wfc_0_consti_"*string(i))
        wfc = reshape(wfc, xDim, yDim)
        calculate_energy(wfc, K, V, Ax, Ay, g, xDim, yDim, dx, dy)
        println()
    end
end

#calculate()
calculate_range(0, 5000, 50000)
