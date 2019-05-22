# Simple script to plot variables. Will grow with time.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse

parser = argparse.ArgumentParser(description='reading strings for plotting')
parser.add_argument('strings', metavar='ID', nargs='+', help='string to plot')

args = parser.parse_args()
print(args.strings)

#class for all variables with initial definitions
class params:

    # Defining static / default values
    xDim = 512
    yDim = 512
    
    # data_dir is assumed to be in the previous directory
    data_dir = "data"
    
    # Defaulting to first element after imaginary time evolution
    start = 0
    end = 1
    incr = 1

    # specifying component to use for plotting
    comp = 0

    # item to work with
    item = "wfc"

# Function to plot specific variable
def plot_var(xDim, yDim, data_dir, pltval):
    if data_dir[0] != "/":
        data_dir = "../" + data_dir
    data = data_dir + "/" + pltval
    lines = np.loadtxt(data)
    val = np.reshape(lines, (xDim,yDim))
    '''
    val = -np.log(val) * 1E4 * 1.0545718E-34
    data_V = "../data/K_0"
    lines_V = np.loadtxt(data_V)
    V_val = np.reshape(lines_V, (xDim, yDim))
    final_val = V_val - val
    '''
    plt.imshow(val, extent=(1,xDim,1,yDim), interpolation='nearest',
                   cmap = cm.jet)
    plt.colorbar()
    fig = plt.gcf()
    #plt.clim(0,1)
    plt.show()

# function to plot a variable with a range
def plot_var_range(xDim, yDim, data_dir, pltval, start, end, incr, comp):
    if data_dir[0] != "/":
        data_dir = "../" + data_dir
    for i in range(start, end, incr):
        print(i)
        output = pltval + "_%s_%s" % (comp, i)
        data = data_dir + "/" + output
        lines = np.loadtxt(data)
        val = np.reshape(lines, (xDim,yDim))
        plt.imshow(val, extent=(1,xDim,1,yDim), interpolation='nearest',
                       cmap = cm.jet)
        plt.colorbar()
        fig = plt.gcf()
        #plt.show()
        plt.draw()
        num_str = "%s" % i
        output_str = pltval + num_str.rjust(5,'0') + ".png"
        fig.savefig(output_str)
        plt.clf()

    
# Function to plot wfc with pltvar as a variable to modify the type of plot
def plot_wfc(xDim, yDim, data_dir, pltval, start, end, incr, comp):
    if data_dir[0] != "/":
        data_dir = "../" + data_dir
    for i in range(start,end,incr):
        print(i)
        data_real = data_dir + "/wfc_%s_const_%s" % (comp, i)
        data_im = data_dir + "/wfc_%s_consti_%s" % (comp, i)
        if pltval == "wfc_ev":
            data_real = data_dir + "/wfc_%s_ev_%s" % (comp, i)
            data_im = data_dir + "/wfc_%s_evi_%s" % (comp, i)

        #data_x = data_dir + "x_0" % i
        #data_y = data_dir + "y_0" % i
        

        #print(i)
    
        lines_real = np.loadtxt(data_real)
        lines_im = np.loadtxt(data_im)
        wfc_real = np.reshape(lines_real, (xDim,yDim));
        wfc_im = np.reshape(lines_im, (xDim,yDim));
    
        wfc = abs(wfc_real + 1j * wfc_im)
        wfc = wfc * wfc
    
        #wfc_k = np.fft.fft2(wfc) 
        #wfc_k_plot = np.abs(np.fft.fftshift(wfc_k))
        #wfc_k_plot = wfc_k_plot**2

        plt.imshow(wfc, extent=(-6.9804018707623236e-04,6.9804018707623236e-04,-6.9804018707623236e-04,6.9804018707623236e-04), interpolation='nearest',
                   cmap = cm.jet)
        plt.colorbar()
        #plt.clim(0,1)
        plt.show()
        #fig = plt.figure()
        #fig.savefig('wfc.png')

# Function to plot complex vals with pltvar as the variable
def plot_complex(xDim, yDim, data_dir, pltval, start, end, incr, comp):
    if data_dir[0] != "/":
        data_dir = "../" + data_dir

    for i in range(start, end, incr):
        data_real = data_dir + "/" + pltval + "%s_%s" % (comp, i)
        data_im = data_dir + "/" + pltval + "i_%s_%s" % (comp, i)
    
        lines_real = np.loadtxt(data_real)
        lines_im = np.loadtxt(data_im)
        wfc_real = np.reshape(lines_real, (xDim,yDim));
        wfc_im = np.reshape(lines_im, (xDim,yDim));
    
        wfc = abs(wfc_real + 1j * wfc_im)
        wfc = wfc * wfc
        
        plt.imshow(wfc, extent=(1,xDim,1,yDim), interpolation='nearest',
                   cmap = cm.jet)
        plt.colorbar()
        plt.show()
        #fig = plt.figure()
        #fig.savefig('wfc.png')

# Function to plot wfc with pltvar as a variable to modify the type of plot
def plot_wfc_k(xDim, yDim, data_dir, pltval, start, end, incr, comp):
    if data_dir[0] != "/":
        data_dir = "../" + data_dir

    for i in range(start,end,incr):
        print(i)
        data_real = data_dir + "/wfc_%s_const_%s" % (comp, i)
        data_im = data_dir + "/wfc_%s_consti_%s" % (comp, i)
        if pltval == "wfc_k_ev":
            data_real = data_dir + "/wfc_%s_ev_%s" % (comp, i)
            data_im = data_dir + "/wfc_%s_evi_%s" % (comp, i)

        lines_real = np.loadtxt(data_real)
        lines_im = np.loadtxt(data_im)
        wfc_real = np.reshape(lines_real, (xDim,yDim));
        wfc_im = np.reshape(lines_im, (xDim,yDim));

        wfc = (wfc_real + 1j * wfc_im)

        wfc_k = np.fft.fft2(wfc)
        wfc_k_plot = np.abs(np.fft.fftshift(wfc_k))
        wfc_k_plot = wfc_k_plot**2

        plt.imshow(wfc_k_plot, extent=(1,xDim,1,yDim), interpolation='nearest',
                   cmap = cm.jet)
        plt.colorbar()
        plt.show()
        #fig = plt.figure()
        #fig.savefig('wfc.png')

# Function to plot wfc with pltvar as a variable to modify the type of plot
def plot_wfc_phase(xDim, yDim, data_dir, pltval, start, end, incr, comp):
    if data_dir[0] != "/":
        data_dir = "../" + data_dir
    for i in range(start,end,incr):
        print(i)
        data_real = data_dir + "/wfc_%s_const_%s" % (comp, i)
        data_im = data_dir + "/wfc_%s_consti_%s" % (comp, i)
        if pltval == "wfc_phase_ev":
            data_real = data_dir + "/wfc_%s_ev_%s" % (comp, i)
            data_im = data_dir + "/wfc_%s_evi_%s" % (comp, i)

        lines_real = np.loadtxt(data_real)
        lines_im = np.loadtxt(data_im)
        wfc_real = np.reshape(lines_real, (xDim,yDim));
        wfc_im = np.reshape(lines_im, (xDim,yDim));

        wfc = (wfc_real + 1j * wfc_im)

        wfc = np.angle(wfc)
        plt.imshow(wfc, extent=(1,xDim,1,yDim), interpolation='nearest',
                   cmap = cm.jet)
        plt.colorbar()
        plt.show()
        #fig = plt.figure()
        #fig.savefig('wfc.png')

# Function to plot wfc cut
def plot_wfc_cut(xDim, yDim, data_dir, pltval, start, end, incr, comp):
    if data_dir[0] != "/":
        data_dir = "../" + data_dir
    for i in range(start,end,incr):
        print(i)
        data_real = data_dir + "/wfc_%s_const_%s" % (comp, i)
        data_im = data_dir + "/wfc_%s_consti_%s" % (comp, i)
        if pltval == "wfc_cut_ev":
            data_real = data_dir + "/wfc_%s_ev_%s" % (comp, i)
            data_im = data_dir + "/wfc_%s_evi_%s" % (comp, i)

        lines_real = np.loadtxt(data_real)
        lines_im = np.loadtxt(data_im)
        wfc_real = np.reshape(lines_real, (xDim,yDim));
        wfc_im = np.reshape(lines_im, (xDim,yDim));

        wfc = abs(wfc_real + 1j * wfc_im)
        wfc = wfc*wfc

        max = 0
        for j in range(xDim):
            for k in range(yDim):
                if (wfc[j][k] > max):
                    max = wfc[j][k]

        print("Max value is: ",max)
        for j in range(xDim):
            for k in range(yDim):
                if (wfc[j][k] > max*0.4):
                    wfc[j][k] = 1.0
                else:
                    wfc[j][k] = 0.0

        plt.imshow(wfc, extent=(-6.9804018707623236e-04,6.9804018707623236e-04,-6.9804018707623236e-04,6.9804018707623236e-04), interpolation='nearest',
                   cmap = cm.jet)
        plt.colorbar()
        plt.show()
        #fig = plt.figure()
        #fig.savefig('wfc.png')



# Function to parse arguments for plotting
# Note: We assume that the parameters come in sets
def parse_args(string_list):
    i = 0
    par = params()
    print(string_list[i])
    while i < len(string_list):
        # -d for "data_dir"
        if (string_list[i] == "d"):
            par.data_dir = string_list[i+1]
            i += 2
        # -i for "item" -- The thing to plot
        elif (string_list[i] == "i"):
            par.item = string_list[i+1]
            i += 2
        # -c for "component" -- The component of the thing to plot
        elif (string_list[i] == "c"):
            par.comp = string_list[i+1]
            i += 2
        # -r for "range"
        elif (string_list[i] == "r"):
            par.start = int(string_list[i+1])
            par.end = int(string_list[i+2]) + 1
            par.incr = int(string_list[i+3])
            par.range = True
            i+= 4
        # -g for "grid"
        elif (string_list[i] == "g"):
            par.xDim = int(string_list[i+1])
            par.yDim = int(string_list[i+2])
            i += 3
    return par

def plot(par):
    if (par.item == "wfc" or par.item == "wfc_ev"):
        plot_wfc(par.xDim, par.yDim, par.data_dir, par.item, 
                 par.start, par.end, par.incr, par.comp)
    elif (par.item == "wfc_k" or par.item == "wfc_k_ev"):
        plot_wfc_k(par.xDim, par.yDim, par.data_dir, par.item,
                   par.start, par.end, par.incr, par.comp)
    elif (par.item == "wfc_phase" or par.item == "wfc_phase_ev"):
        plot_wfc_phase(par.xDim, par.yDim, par.data_dir, par.item,
                       par.start, par.end, par.incr, par.comp)
    elif (par.item == "GK" or par.item == "GV"):
        plot_complex(par.xDim, par.yDim, par.data_dir, par.item,
                     par.start, par.end, par.incr, par.comp)
    elif (par.item == "wfc_cut" or par.item == "wfc_cut_ev"):
        plot_wfc_cut(par.xDim, par.yDim, par.data_dir, par.item, 
                     par.start, par.end, par.incr, par.comp)
    elif (par.end != 1):
        plot_var_range(par.xDim, par.yDim, par.data_dir, par.item,
                       par.start, par.end, par.incr, par.comp)
    else:
        plot_var(par.xDim, par.yDim, par.data_dir, par.item)

par = parse_args(args.strings)
plot(par)
