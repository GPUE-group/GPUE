import h5py
import numpy as np

h5py.get_config().complex_names = ("re", "im")

def write_var(shape, h5path, filename_real, filename_imag=None, filename_output="../data/input.h5"):
    var = np.loadtxt(filename_real)
    var = np.reshape(var, shape)
    if (filename_imag is not None):
        var_im = np.loadtxt(filename_imag)
        var_im = np.reshape(var_im, shape)
        var = var.astype("complex128") + var_im.astype("complex128") * 1j

    f = h5py.File(filename_output, "a")
    f[h5path] = var
    f.close()

def write_wfc(shape, gstate, i, data_real, data_im, filename_output="../data/input.h5"):
    dset_name = "/WFC/{}/{}".format("CONST" if gstate else "EV", i)
    write_var(shape, dset_name, data_real, data_im, filename_output)

def write_gauge(shape, data_ax, data_ay, data_az, filename_output="../data/input.h5"):
    write_var(shape, "/A/AX/0", data_ax, None, filename_output)
    write_var(shape, "/A/AY/0", data_ay, None, filename_output)
    write_var(shape, "/A/AZ/0", data_az, None, filename_output)

if __name__ == "__main__":
    # write_gauge((2,2,2), "../data/ax_re.txt", "../data/ay_re.txt", "../data/az_re.txt")
    write_wfc((4,4), True, 1234, "../data/wfc_re.txt", "../data/wfc_im.txt")

