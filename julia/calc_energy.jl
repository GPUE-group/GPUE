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
    energy_r = wfc_c.* (H_r .+ 0.5*g*abs2.(wfc)).* wfc

    #energy_l = wfc_c.*(-omega*omegaX*(ifft(-Ay.*fft(wfc,2),2) + ifft(Ax.*fft(wfc,1),1)))
    #energy_l = wfc_c.*(-(Ax.*ifft(Ay.*fft(wfc,2),2) + Ay.*ifft(Ax.*fft(wfc,1),1)))
    #energy_l = wfc_c.*(ifft(Ax.*fft(wfc,2),2) + ifft(Ay.*fft(wfc,1),1))
    #energy_l = wfc_c.*ifft(Ax.*fft(ifft(Ay.*fft(wfc,1),1),2),2)
    energy_l = wfc_c.*(im*hbar*(Ax.*derive(wfc,1,dx) + Ay.*derive(wfc,2,dy)))

    # Integrating over all space
    energy_final = 0
    for i = 1:xDim*yDim
        #energy_final += real(energy_r[i] + energy_k[i] + energy_l[i])
        energy_final += energy_l[i]
    end 

    return energy_final*dx*dy
end

function calculate()
    xDim = 1024
    yDim = 1024
    dx = 5.11523e-07
    dy = 5.11523e-07
    omega = 0.3
    omegaX = 6.283
    data_dir = "./"
    Ax = readdlm(data_dir*"Ax_0")
    Ay = readdlm(data_dir*"Ay_0")
    K = readdlm(data_dir*"K_0")
    V = readdlm(data_dir*"V_0")

    wfc = readdlm("wfc_load") + readdlm("wfci_load")*im

    Ax = reshape(Ax, xDim, yDim)
    Ay = reshape(Ay, xDim, yDim)
    #Ax /= omega*omegaX
    #Ay /= omega*omegaX
    K = reshape(K, xDim, yDim)
    V = reshape(V, xDim, yDim)
    wfc = reshape(wfc, xDim, yDim)
    g = 6.80274e-40

    energy = calculate_energy(wfc, K, V, Ax, Ay, g, xDim, yDim, dx, dy)

    println(energy)

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
    #Ax /= omega*omegaX
    #Ay /= omega*omegaX
    K = reshape(K, xDim, yDim)
    V = reshape(V, xDim, yDim)
    g = 6.80274e-41

    for i = start:incr:final
        wfc = readdlm(data_dir*"wfc_0_const_"*string(i))
              + readdlm(data_dir*"wfc_0_consti_"*string(i))*im
        wfc = reshape(wfc, xDim, yDim)
        energy = calculate_energy(wfc, K, V, Ax, Ay, g, xDim, yDim, dx, dy)
        println(energy)
    end
end

#calculate()
#calculate_range(0, 5000, 50000)
