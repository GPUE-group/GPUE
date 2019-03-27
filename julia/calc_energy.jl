using DelimitedFiles
using FFTW

# We are calculating the energy to check <Psi|H|Psi>
function calculate_energy(wfc, H_k, H_r, Ax, Ay, g, xDim, yDim, dx, dy)
    # Creating momentum and conjugate wavefunctions
    wfc_k = fft(wfc)
    wfc_c = conj(wfc)

    # Finding the momentum and real-space energy terms
    energy_k = wfc_c.*ifft((H_k) .* wfc_k)
    energy_r = wfc_c.* (H_r .+ 0.5*g*abs2.(wfc)).* wfc

    #energy_l = wfc_c.*(-(ifft(Ay.*fft(wfc,2),2) - ifft(Ax.*fft(wfc,1),1)))
    #energy_l = wfc_c.*(ifft(Ax.*fft(wfc,2),2) + ifft(Ay.*fft(wfc,1),1))
    #energy_l = wfc_c.*ifft(Ax.*fft(ifft(Ay.*fft(wfc,1),1),2),2)

    # Integrating over all space
    energy_final = 0
    for i = 1:xDim*yDim
        energy_final += real(energy_r[i] + energy_r[i] + energy_l[i])
    end 

    return energy_final*dx*dy
end

xDim = 1024
yDim = 1024
dx = 5.11523e-07
dy = 5.11523e-07
wfc = readdlm("wfc_ev_0") + readdlm("wfc_evi_0")*im
#wfc = readdlm("wfc_load") + readdlm("wfci_load")*im
Ax = readdlm("Ax_0")
Ay = readdlm("Ay_0")
K = readdlm("K_0")
V = readdlm("V_0")

wfc = reshape(wfc, xDim, yDim)
Ax = reshape(Ax, xDim, yDim)
Ay = reshape(Ay, xDim, yDim)
K = reshape(K, xDim, yDim)
V = reshape(V, xDim, yDim)
g = 6.80274e-40

energy = calculate_energy(wfc, K, V, Ax, Ay, g, xDim, yDim, dx, dy)
println(energy)
