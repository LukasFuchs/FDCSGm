# 
#         hawaiiS
#                   www.fabiocrameri.ch/colourmaps
from matplotlib.colors import LinearSegmentedColormap      
      
cm_data = [[0.55054, 0.006842, 0.45198],      
           [0.70378, 0.94898, 0.99377],      
           [0.61115, 0.58982, 0.11132],      
           [0.59051, 0.3043, 0.24172],      
           [0.41971, 0.83304, 0.55824],      
           [0.57371, 0.18114, 0.3346],      
           [0.53707, 0.73968, 0.29114],      
           [0.60713, 0.4336, 0.1582],      
           [0.40326, 0.91353, 0.83451],      
           [0.47814, 0.78987, 0.42206],      
           [0.58235, 0.24327, 0.28591],      
           [0.37617, 0.87491, 0.69869],      
           [0.58608, 0.67243, 0.17528],      
           [0.56366, 0.11287, 0.39002],      
           [0.59872, 0.36683, 0.19953],      
           [0.61384, 0.50759, 0.12051],      
           [0.52972, 0.93842, 0.93259],      
           [0.57814, 0.21254, 0.30954],      
           [0.37572, 0.8951, 0.76855],      
           [0.61664, 0.94539, 0.96756],      
           [0.61447, 0.54785, 0.10911],      
           [0.55762, 0.06997, 0.42038],      
           [0.61099, 0.46954, 0.13827],      
           [0.58646, 0.27377, 0.26344],      
           [0.39395, 0.85406, 0.6281],      
           [0.59458, 0.33523, 0.22051],      
           [0.44844, 0.8118, 0.48955],      
           [0.50799, 0.76629, 0.35578],      
           [0.56381, 0.7086, 0.22974],      
           [0.56894, 0.14842, 0.36138],      
           [0.60293, 0.39949, 0.17878],      
           [0.60203, 0.63211, 0.13369],      
           [0.45254, 0.92703, 0.88545],      
           [0.59492, 0.65266, 0.15242],      
           [0.59255, 0.3197, 0.23106],      
           [0.57138, 0.16501, 0.34776],      
           [0.58027, 0.22798, 0.29757],      
           [0.60082, 0.38301, 0.18915],      
           [0.61342, 0.56868, 0.10794],      
           [0.6126, 0.48829, 0.12892],      
           [0.60741, 0.61106, 0.11981],      
           [0.46322, 0.80096, 0.45567],      
           [0.61453, 0.52746, 0.11366],      
           [0.55423, 0.04117, 0.43606],      
           [0.56075, 0.092811, 0.40501],      
           [0.57565, 0.69112, 0.20131],      
           [0.66055, 0.94751, 0.98145],      
           [0.49309, 0.77836, 0.38876],      
           [0.4891, 0.9333, 0.91062],      
           [0.59664, 0.35092, 0.20999],      
           [0.60504, 0.41635, 0.16844],      
           [0.38553, 0.90462, 0.80228],      
           [0.38359, 0.86453, 0.66337],      
           [0.60913, 0.45132, 0.14807],      
           [0.37304, 0.88514, 0.73387],      
           [0.58848, 0.28901, 0.25249],      
           [0.55085, 0.72477, 0.25987],      
           [0.58442, 0.25853, 0.27455],      
           [0.5227, 0.75346, 0.32319],      
           [0.40624, 0.84356, 0.59304],      
           [0.56638, 0.13117, 0.37547],      
           [0.57597, 0.19698, 0.32185],      
           [0.43386, 0.82246, 0.52375],      
           [0.57269, 0.94241, 0.95146],      
           [0.4217, 0.91968, 0.85731],      
           [0.56224, 0.10313, 0.39747],      
           [0.58339, 0.25092, 0.2802],      
           [0.52994, 0.7467, 0.30708],      
           [0.59768, 0.35883, 0.20478],      
           [0.57706, 0.20478, 0.31565],      
           [0.58544, 0.26616, 0.26898],      
           [0.59461, 0.94402, 0.95982],      
           [0.60187, 0.39122, 0.18395],      
           [0.59356, 0.32743, 0.22577],      
           [0.48561, 0.78418, 0.40538],      
           [0.55243, 0.023795, 0.444],      
           [0.57017, 0.15681, 0.35452],      
           [0.60814, 0.44241, 0.15309],      
           [0.59561, 0.34305, 0.21523],      
           [0.42672, 0.82776, 0.54096],      
           [0.44111, 0.81714, 0.50662],      
           [0.39338, 0.90916, 0.81862],      
           [0.58105, 0.68192, 0.18799],      
           [0.56988, 0.70002, 0.21528],      
           [0.59869, 0.64247, 0.1425],      
           [0.37967, 0.89992, 0.78557],      
           [0.61461, 0.53759, 0.11103],      
           [0.61183, 0.47884, 0.13351],      
           [0.54404, 0.73238, 0.27539],      
           [0.59977, 0.37488, 0.19437],      
           [0.56768, 0.13987, 0.36838],      
           [0.61327, 0.49788, 0.12457],      
           [0.58747, 0.28137, 0.25793],      
           [0.45581, 0.80641, 0.47258],      
           [0.37398, 0.88005, 0.71632],      
           [0.61245, 0.57923, 0.10903],      
           [0.50055, 0.77241, 0.37222],      
           [0.59153, 0.31197, 0.23638],      
           [0.60493, 0.62162, 0.12612],      
           [0.60608, 0.42493, 0.16332]]      
      
hawaiiS_map = LinearSegmentedColormap.from_list('hawaiiS', cm_data)      
# For use of "viscm view"      
test_cm = hawaiiS_map      
      
if __name__ == "__main__":      
    import matplotlib.pyplot as plt      
    import numpy as np      
      
    try:      
        from viscm import viscm      
        viscm(hawaiiS_map)      
    except ImportError:      
        print("viscm not found, falling back on simple display")      
        plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',      
                   cmap=hawaiiS_map)      
    plt.show()      
