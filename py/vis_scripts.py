from gen_data import *
'''
proj_var2d(xDim, yDim, zDim, "data", "fieldz")
proj_2d(xDim, yDim, zDim,"data","wfc",50000)
#proj_2d(xDim, yDim, zDim,"data","wfc",1000000)
#proj_k2d(xDim, yDim, zDim,"data","wfc",500000)

item_wfc = wfc_density(xDim, yDim, zDim,"data","wfc",50000)
#item_wfc = wfc_density(xDim, yDim, zDim,"data","wfc",1000000)
item_phase = wfc_phase(xDim, yDim, zDim,"data","wfc",50000)
#item_phase = wfc_phase(xDim, yDim, zDim,"data","wfc",1000000)
#item_var = var(xDim, yDim, zDim,"data","By_0")
#item2_var = var(xDim, yDim, zDim,"data","Edges_0")

#to_bvox(item_wfc, xDim, yDim, zDim, 1, "test_wfc.bvox")
#to_bvox(item_phase, xDim, yDim, zDim, 1, "test_phase.bvox")
#to_bvox(item_var, xDim, yDim, zDim, 1, "test_var.bvox")
#to_bvox(item2_var, xDim, yDim, zDim, 1, "test_edges.bvox")
'''
'''

print("iterating through data...")
# Writing out the potential
#item_pot = var(xDim, yDim, zDim, "data", "V_0")
#to_bvox(item_pot, xDim, yDim, zDim, 1, "test_pot.bvox")
#item_gauge = var(xDim, yDim, zDim, "data", "Az_0")
#to_bvox(item_gauge, xDim, yDim, zDim, 1, "test_gauge.bvox")
#item_edges = var(xDim, yDim, zDim,"data","Edges_0")
#to_bvox(item_edges, xDim, yDim, zDim, 1, "test_edges.bvox")
#to_vtk(item_edges, xDim, yDim, zDim, 1, "test_edges.vtk")
for i in range(0,101):
#for i in [86, 143, 195, 95, 152, 204]:
#for i in [4, 40]:
    num = i * 100
    proj_2d(xDim, yDim, zDim, "data", "wfc_ev", num)
    #item = wfc_density(xDim, yDim, zDim, "data", "wfc_ev", num)
    #item_ph = wfc_phase(xDim, yDim, zDim, "data", "wfc", num)
    #to_bvox(item, xDim, yDim, zDim, 1, "wfc_%s.bvox" %num)
    #to_vtk(item, xDim, yDim, zDim, 1, "wfc_%s.vtk" %num)
    #to_bvox(item_ph, xDim, yDim, zDim, 1, "wfc_ph_%s.bvox" %num)

'''
'''

for i in range(0,21):
    num = i * 1
    print(num)
    item = var(xDim, yDim, zDim, "data", "Edges_%s"%num)
    to_bvox(item, xDim, yDim, zDim, 1, "EDGES_%s.bvox" %num)

'''
'''
proj_var2d(xDim, yDim, zDim, "data", "fieldz")
item = var(xDim, yDim, zDim, "data", "fieldz")
to_bvox(item, xDim, yDim, zDim, 1, "HE11check.bvox")
'''

'''
proj_2d(xDim, yDim, zDim, "data", "wfc_ramp", 100000)
#item_edges = var(xDim, yDim, zDim,"data","Edges_0")
#to_bvox(item_edges, xDim, yDim, zDim, 1, "test_edges.bvox")
'''

'''
print("start")
item_edges = var(xDim, yDim, zDim,"data","Edges_std")
to_bvox(item_edges, xDim, yDim, zDim, 1, "edges_std.bvox")
to_vtk(item_edges, xDim, yDim, zDim, 1, "edges_std.vtk")
print("starting linpol")
item_edges = var(xDim, yDim, zDim,"data","Edges_linpol")
to_bvox(item_edges, xDim, yDim, zDim, 1, "edges_linpol.bvox")
to_vtk(item_edges, xDim, yDim, zDim, 1, "edges_linpol.vtk")
'''

'''
item = wfc_density(xDim, yDim, zDim, "data", "wfc", 100000)
to_vtk(item, xDim, yDim, zDim, 1, "wfc.vtk")
to_bvox(item, xDim, yDim, zDim, 1, "wfc_check.bvox")
'''

#proj_phase_2d(xDim, yDim, zDim, "data", "wfc", 660)
#proj_phase_2d(xDim, yDim, zDim, "data", "wfc", 320)
#proj_phase_2d(xDim, yDim, zDim, "data", "wfc", 1000)
#item_pot = var(xDim, yDim, zDim, "data", "V_0")
#to_bvox(item_pot, xDim, yDim, zDim, 1, "test_pot.bvox")
#item_edges = var(xDim, yDim, zDim,"data","Edges_0")
#to_bvox(item_edges, xDim, yDim, zDim, 1, "test_edges.bvox")
#item = wfc_density(xDim, yDim, zDim,"data","wfc", 0)
#to_bvox(item, xDim, yDim, zDim, 1, "test_wfc.bvox")

'''
item_gx = var(xDim, yDim, zDim,"data","Ax_0")
to_bvox(item_gx, xDim, yDim, zDim, 1, "Ax_0.bvox")
item_gy = var(xDim, yDim, zDim,"data","Ay_0")
to_bvox(item_gy, xDim, yDim, zDim, 1, "Ay_0.bvox")
item_gz = var(xDim, yDim, zDim,"data","Az_0")
to_bvox(item_gz, xDim, yDim, zDim, 1, "Az_0.bvox")
item_V = var(xDim, yDim, zDim,"data","V_0")
to_bvox(item_V, xDim, yDim, zDim, 1, "V_0.bvox")
item_K = var(xDim, yDim, zDim,"data","K_0")
to_bvox(item_K, xDim, yDim, zDim, 1, "K_0.bvox")
'''
'''
for i in range(0,11):
    num = i*100000
    print(num)
    item = var(xDim, yDim, zDim, "data", "Edges_%s" %num)
    to_bvox(item, xDim, yDim, zDim, 1, "EDGES_%s.bvox" %num)

'''
'''
num = 500000
proj_2d(xDim, yDim, zDim, "data", "wfc", num)
item = wfc_density(xDim, yDim, zDim, "data", "wfc", num)
item_edges = var(xDim, yDim, zDim, "data", "Edges_0")
to_bvox(item, xDim, yDim, zDim, 1, "wfc_%s.bvox" %num)
to_bvox(item_edges, xDim, yDim, zDim, 1, "Edges_0.bvox")
to_vtk(item, xDim, yDim, zDim, 1, "wfc_%s.vtk" %num)

'''
'''
proj_var2d(xDim, yDim, zDim, "data", "V_0")
#proj_2d(xDim, yDim, zDim, "data", "wfc_ev", 0)
'''
'''
item = var(xDim, yDim, zDim,"data","Bx_0")
to_vtk(item, xDim, yDim, zDim, 1, "test_Bx.vtk")
item = var(xDim, yDim, zDim,"data","By_0")
to_vtk(item, xDim, yDim, zDim, 1, "test_By.vtk")
item = var(xDim, yDim, zDim,"data","Bz_0")
to_vtk(item, xDim, yDim, zDim, 1, "test_Bz.vtk")
item = var(xDim, yDim, zDim,"data","Br_0")
to_vtk(item, xDim, yDim, zDim, 1, "test_Br.vtk")
item = var(xDim, yDim, zDim,"data","Bphi_0")
to_vtk(item, xDim, yDim, zDim, 1, "test_Bphi.vtk")
#item = var(xDim, yDim, zDim,"data","Edges_0")
#to_bvox(item, xDim, yDim, zDim, 1, "test_wfc.bvox")

#comx, comy = wfc_com(xDim, yDim, zDim, "data", "wfc", 0)
#print(comx)
#print(comy)
'''

'''
for i in range(0,3):
    if (i == 2):
        data_dir = "data_HE21lin2_w50"
    if (i == 1):
        data_dir = "data_lin2_w50_real_time"
    if (i == 0):
        data_dir = "data_broad2_w50_realtime1"

    item = proj_var2d(xDim, yDim, zDim, data_dir, "Bx_0", "BX2D.dat")
    item = proj_var2d(xDim, yDim, zDim, data_dir, "By_0", "BY2D.dat")
    item = proj_var2d(xDim, yDim, zDim, data_dir, "Bz_0", "BZ2D.dat")
    item = proj_var2d(xDim, yDim, zDim, data_dir, "Br_0", "BR2D.dat")
    item = proj_var2d(xDim, yDim, zDim, data_dir, "Bphi_0", "BPHI2D.dat")
    print(data_dir + " is done")
'''
#item = proj_var1d(xDim, yDim, zDim, "data_check", "Br_0", "BR1D.dat")
#item = proj_2d(xDim, yDim, zDim, "data_check", "wfc", 0)
#item = var(xDim, yDim, zDim,"data_broad2_w50_realtime1","Edges_0")
'''
for i in range(0,10):
    print(i)
    item = proj_var2d(xDim, yDim, zDim,"data_check","Edges_%s"%(i*100), "Edges_%s.dat"%(i*100))
    #to_vtk(item, xDim, yDim, zDim, 1, "test_Edges%s.vtk"%(i*100))
'''

'''
xDim = 128
yDim = 128
zDim = 128

item = var(xDim, yDim, zDim,"data_lin","Az_0")
to_vtk(item, xDim, yDim, zDim, 1, "test_Az.vtk")

item = var(xDim, yDim, zDim,"data_lin","Ax_0")
to_vtk(item, xDim, yDim, zDim, 1, "test_Ax.vtk")

item = var(xDim, yDim, zDim,"data_lin","Ay_0")
to_vtk(item, xDim, yDim, zDim, 1, "test_Ay.vtk")
'''

#item = var_r2(xDim, yDim, zDim,"/media/james/ExtraDrive1/GPUE/data_lin")
#to_vtk(item, xDim, yDim, zDim, 1, "test_Ar.vtk")

#item = var(xDim, yDim, zDim,"/media/james/ExtraDrive1/GPUE/data_elong_w0.5n",
#           "Edges_0")
#to_vtk(item, xDim, yDim, zDim, 1, "test_Edges.vtk")
#item = var(xDim, yDim, zDim,"/media/james/ExtraDrive1/GPUE/data_lin","Edges_0")

item = var(xDim, yDim, zDim,"/media/james/ExtraDrive1/GPUE/data_elong_w0.5n","Edges_0")
to_vtk(item, xDim, yDim, zDim, 1, "test_Edges.vtk")
item = var(xDim, yDim, zDim,"/media/james/ExtraDrive1/GPUE/data_elong_w0.5n","V_0")
to_vtk(item, xDim, yDim, zDim, 1, "test_V.vtk")
item = var(xDim, yDim, zDim,"/media/james/ExtraDrive1/GPUE/data_elong_w0.5n","Az_0")
to_vtk(item, xDim, yDim, zDim, 1, "test_Az.vtk")
#to_bvox(item, xDim, yDim, zDim, 1, "test_Edges.bvox")

#thresh = find_thresh(xDim, yDim, "/media/james/ExtraDrive1/GPUE/data_elong",
#                     "Pwfc_", 0, 0.5)
#print(thresh)
#angle = find_angle(xDim, yDim, "/media/james/ExtraDrive1/GPUE/data_elong",
#                   "Pwfc_", 0, thresh)
#print(angle)
