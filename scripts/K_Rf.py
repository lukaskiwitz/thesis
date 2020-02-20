#!/home/kiwitz/anaconda3/envs/fenicsproject/bin/python

from parameter_scan_test import get_parameter_space, p_c

name = "K_Rf"
# fraction of receptors on secreting cells
pList = get_parameter_space("fraction", 1)#quick fix

def get_Rf_from_K(K):
    f = p_c["fraction"]
    Rs = p_c["R_l"] #quick fix
    K = [0.99 if i == 1 else i for i in K]
    return  (-(-1+f)**2*K*Rs)/(f**2*(-1+K))

# f_list = get_Rf_from_K(pList)

# plt.plot(pList,f_list)
# plt.show()
# path = path_prefix+name+"/"
# scan = [{"R_il2_f": i[0]} for i in f_list]
# execute_scan(scan, p_c, T, path, ext_cache)