# 
#         glasgowS
#                   www.fabiocrameri.ch/colourmaps
from matplotlib.colors import LinearSegmentedColormap      
      
cm_data = [[0.21181, 0.073933, 0.22061],      
           [0.85897, 0.82865, 1],      
           [0.42748, 0.44107, 0.17541],      
           [0.44013, 0.17727, 0.021779],      
           [0.41302, 0.63858, 0.64524],      
           [0.32208, 0.10509, 0.1252],      
           [0.62636, 0.73453, 0.82848],      
           [0.45474, 0.31761, 0.0044146],      
           [0.38673, 0.54066, 0.41656],      
           [0.51666, 0.68948, 0.74447],      
           [0.40689, 0.49213, 0.29625],      
           [0.26968, 0.089618, 0.17115],      
           [0.38192, 0.12619, 0.07636],      
           [0.44576, 0.38265, 0.062145],      
           [0.73189, 0.77731, 0.90777],      
           [0.45688, 0.25015, 0.0021644],      
           [0.37458, 0.58916, 0.53562],      
           [0.43741, 0.41306, 0.1165],      
           [0.39673, 0.51653, 0.35654],      
           [0.79113, 0.80128, 0.95154],      
           [0.4515, 0.35047, 0.020328],      
           [0.41438, 0.14644, 0.048631],      
           [0.67852, 0.75568, 0.86778],      
           [0.57249, 0.71259, 0.78766],      
           [0.46115, 0.66486, 0.69762],      
           [0.41715, 0.46718, 0.2358],      
           [0.45334, 0.21392, 0.0074732],      
           [0.45647, 0.28442, 0.0009539],      
           [0.29581, 0.097347, 0.14796],      
           [0.24219, 0.081714, 0.19539],      
           [0.37796, 0.56466, 0.47616],      
           [0.35037, 0.1139, 0.1017],      
           [0.38308, 0.61174, 0.58769],      
           [0.36575, 0.11935, 0.089337],      
           [0.42231, 0.45432, 0.20554],      
           [0.28284, 0.093477, 0.15943],      
           [0.22749, 0.07767, 0.20791],      
           [0.33585, 0.10931, 0.11364],      
           [0.38204, 0.55265, 0.44641],      
           [0.30883, 0.10118, 0.13657],      
           [0.39525, 0.62511, 0.61695],      
           [0.39169, 0.52863, 0.38659],      
           [0.39842, 0.13503, 0.062973],      
           [0.45573, 0.30108, 0.0018503],      
           [0.44187, 0.39814, 0.08848],      
           [0.44831, 0.19533, 0.013206],      
           [0.59971, 0.72369, 0.80834],      
           [0.41201, 0.47975, 0.26605],      
           [0.25617, 0.085604, 0.18312],      
           [0.54475, 0.7012, 0.76641],      
           [0.7605, 0.7889, 0.92902],      
           [0.45592, 0.2323, 0.0041464],      
           [0.65254, 0.74515, 0.84822],      
           [0.45337, 0.33407, 0.0098688],      
           [0.42864, 0.16066, 0.033747],      
           [0.37513, 0.57679, 0.50586],      
           [0.82403, 0.81457, 0.97536],      
           [0.48858, 0.67739, 0.72163],      
           [0.43255, 0.42736, 0.14564],      
           [0.43548, 0.6519, 0.6722],      
           [0.449, 0.3667, 0.038112],      
           [0.70476, 0.76632, 0.88748],      
           [0.4018, 0.50438, 0.3264],      
           [0.45692, 0.26749, 0.0011036],      
           [0.37646, 0.59869, 0.55798],      
           [0.4541, 0.32583, 0.0066794],      
           [0.31541, 0.10316, 0.13093],      
           [0.37989, 0.55864, 0.46129],      
           [0.43003, 0.43428, 0.16043],      
           [0.26298, 0.087669, 0.17711],      
           [0.41457, 0.47348, 0.25092],      
           [0.45655, 0.2413, 0.003037],      
           [0.6131, 0.72914, 0.81847],      
           [0.2198, 0.075756, 0.21424],      
           [0.40654, 0.14038, 0.055873],      
           [0.66553, 0.75042, 0.858],      
           [0.43479, 0.16873, 0.027345],      
           [0.45613, 0.29278, 0.0012453],      
           [0.45251, 0.34228, 0.014455],      
           [0.53073, 0.69538, 0.75554],      
           [0.3942, 0.52258, 0.37156],      
           [0.55867, 0.70693, 0.77712],      
           [0.28934, 0.095451, 0.15367],      
           [0.39017, 0.13034, 0.069771],      
           [0.37634, 0.5707, 0.49102],      
           [0.27629, 0.091604, 0.16524],      
           [0.448, 0.65844, 0.6851],      
           [0.3892, 0.53464, 0.40159],      
           [0.50258, 0.68348, 0.73318],      
           [0.37376, 0.12251, 0.082919],      
           [0.42374, 0.64527, 0.6589],      
           [0.47473, 0.67118, 0.7098],      
           [0.71818, 0.77176, 0.89753],      
           [0.4439, 0.39046, 0.075044],      
           [0.42179, 0.1532, 0.041226],      
           [0.45674, 0.27603, 0.0009164],      
           [0.45527, 0.30938, 0.0028627],      
           [0.43502, 0.42027, 0.13095],      
           [0.45033, 0.3586, 0.028028],      
           [0.3745, 0.58293, 0.52074]]      
      
glasgowS_map = LinearSegmentedColormap.from_list('glasgowS', cm_data)      
# For use of "viscm view"      
test_cm = glasgowS_map      
      
if __name__ == "__main__":      
    import matplotlib.pyplot as plt      
    import numpy as np      
      
    try:      
        from viscm import viscm      
        viscm(glasgowS_map)      
    except ImportError:      
        print("viscm not found, falling back on simple display")      
        plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',      
                   cmap=glasgowS_map)      
    plt.show()      
