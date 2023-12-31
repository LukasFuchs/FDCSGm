# 
#         lapazS
#                   www.fabiocrameri.ch/colourmaps
from matplotlib.colors import LinearSegmentedColormap      
      
cm_data = [[0.10352, 0.047787, 0.39353],      
           [0.99706, 0.94979, 0.95121],      
           [0.36082, 0.54989, 0.64017],      
           [0.17465, 0.32419, 0.57638],      
           [0.70348, 0.67492, 0.58747],      
           [0.23932, 0.44474, 0.6286],      
           [0.91657, 0.8067, 0.72519],      
           [0.13975, 0.19504, 0.49367],      
           [0.51987, 0.61974, 0.60988],      
           [0.12359, 0.12648, 0.4447],      
           [0.60733, 0.64453, 0.59016],      
           [0.20114, 0.3857, 0.60689],      
           [0.97468, 0.88337, 0.84014],      
           [0.81459, 0.72884, 0.62938],      
           [0.29243, 0.5002, 0.64018],      
           [0.15554, 0.26061, 0.53818],      
           [0.43529, 0.58847, 0.6295],      
           [0.13199, 0.16132, 0.46961],      
           [0.1867, 0.35524, 0.59266],      
           [0.56309, 0.63241, 0.59923],      
           [0.47722, 0.60531, 0.62039],      
           [0.98876, 0.9175, 0.89625],      
           [0.65365, 0.65791, 0.58514],      
           [0.2639, 0.47302, 0.63575],      
           [0.75762, 0.69806, 0.60115],      
           [0.11424, 0.089643, 0.41924],      
           [0.14749, 0.22811, 0.51661],      
           [0.86964, 0.76622, 0.67213],      
           [0.21855, 0.41558, 0.61892],      
           [0.95153, 0.84642, 0.78261],      
           [0.16443, 0.29264, 0.55815],      
           [0.32486, 0.52597, 0.64168],      
           [0.39466, 0.569, 0.6364],      
           [0.11911, 0.10838, 0.43203],      
           [0.18039, 0.33979, 0.58476],      
           [0.14363, 0.21165, 0.50529],      
           [0.58503, 0.63844, 0.59436],      
           [0.34246, 0.53819, 0.6413],      
           [0.13589, 0.17829, 0.48176],      
           [0.30818, 0.51329, 0.6413],      
           [0.96444, 0.86528, 0.81152],      
           [0.54139, 0.62623, 0.60448],      
           [0.78598, 0.71249, 0.61334],      
           [0.99344, 0.93379, 0.92382],      
           [0.45613, 0.59723, 0.62518],      
           [0.27766, 0.48677, 0.63831],      
           [0.41476, 0.57907, 0.63327],      
           [0.63014, 0.65093, 0.58696],      
           [0.84275, 0.74689, 0.64912],      
           [0.49849, 0.61279, 0.61523],      
           [0.73003, 0.68559, 0.59265],      
           [0.19361, 0.37054, 0.60004],      
           [0.25112, 0.459, 0.6325],      
           [0.12791, 0.14407, 0.45724],      
           [0.93562, 0.82683, 0.75368],      
           [0.15981, 0.2767, 0.54838],      
           [0.16935, 0.30849, 0.5675],      
           [0.20945, 0.40074, 0.61318],      
           [0.98265, 0.90074, 0.86839],      
           [0.1091, 0.06976, 0.4064],      
           [0.15142, 0.2444, 0.52758],      
           [0.2285, 0.43025, 0.62407],      
           [0.89445, 0.78634, 0.69774],      
           [0.67805, 0.6658, 0.58515],      
           [0.37993, 0.56104, 0.63828],      
           [0.90588, 0.79651, 0.71129],      
           [0.57402, 0.63543, 0.59673],      
           [0.74371, 0.69159, 0.59646],      
           [0.37029, 0.55553, 0.63932],      
           [0.10638, 0.059148, 0.39996],      
           [0.15765, 0.26868, 0.54333],      
           [0.28491, 0.49353, 0.63933],      
           [0.44567, 0.59294, 0.62741],      
           [0.23377, 0.43753, 0.62641],      
           [0.55222, 0.62935, 0.60182],      
           [0.48784, 0.60914, 0.61784],      
           [0.77173, 0.70503, 0.60677],      
           [0.14558, 0.21989, 0.51099],      
           [0.33356, 0.53214, 0.64158],      
           [0.2139, 0.40817, 0.61612],      
           [0.99126, 0.92569, 0.91007],      
           [0.96988, 0.87442, 0.82587],      
           [0.94396, 0.83671, 0.76812],      
           [0.24508, 0.4519, 0.63064],      
           [0.69062, 0.67019, 0.58599],      
           [0.95835, 0.85595, 0.79708],      
           [0.19007, 0.3629, 0.59641],      
           [0.12999, 0.1527, 0.46344],      
           [0.35154, 0.5441, 0.64083],      
           [0.18348, 0.34753, 0.58876],      
           [0.8564, 0.75642, 0.66025],      
           [0.14171, 0.20336, 0.49951],      
           [0.15344, 0.25253, 0.53293],      
           [0.40464, 0.57411, 0.63491],      
           [0.12136, 0.11751, 0.43837],      
           [0.16684, 0.30057, 0.56289],      
           [0.59612, 0.64146, 0.59217],      
           [0.13398, 0.16983, 0.47572],      
           [0.53062, 0.62303, 0.60718],      
           [0.92649, 0.81681, 0.73935]]      
      
lapazS_map = LinearSegmentedColormap.from_list('lapazS', cm_data)      
# For use of "viscm view"      
test_cm = lapazS_map      
      
if __name__ == "__main__":      
    import matplotlib.pyplot as plt      
    import numpy as np      
      
    try:      
        from viscm import viscm      
        viscm(lapazS_map)      
    except ImportError:      
        print("viscm not found, falling back on simple display")      
        plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',      
                   cmap=lapazS_map)      
    plt.show()      
