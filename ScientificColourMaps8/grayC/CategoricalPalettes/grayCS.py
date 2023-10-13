# 
#         grayCS
#                   www.fabiocrameri.ch/colourmaps
from matplotlib.colors import LinearSegmentedColormap      
      
cm_data = [[0.028493, 0.028493, 0.028493],      
           [0.97733, 0.97733, 0.97733],      
           [0.46781, 0.46781, 0.46781],      
           [0.27163, 0.27163, 0.27163],      
           [0.67708, 0.67708, 0.67708],      
           [0.81433, 0.81433, 0.81433],      
           [0.16127, 0.16127, 0.16127],      
           [0.37425, 0.37425, 0.37425],      
           [0.56183, 0.56183, 0.56183],      
           [0.42047, 0.42047, 0.42047],      
           [0.61505, 0.61505, 0.61505],      
           [0.10622, 0.10622, 0.10622],      
           [0.74052, 0.74052, 0.74052],      
           [0.21544, 0.21544, 0.21544],      
           [0.32242, 0.32242, 0.32242],      
           [0.88976, 0.88976, 0.88976],      
           [0.51303, 0.51303, 0.51303],      
           [0.93262, 0.93262, 0.93262],      
           [0.44428, 0.44428, 0.44428],      
           [0.34864, 0.34864, 0.34864],      
           [0.24377, 0.24377, 0.24377],      
           [0.77661, 0.77661, 0.77661],      
           [0.075266, 0.075268, 0.075267],      
           [0.64536, 0.64536, 0.64536],      
           [0.53527, 0.53527, 0.53527],      
           [0.3961, 0.3961, 0.3961],      
           [0.84875, 0.84875, 0.84875],      
           [0.29897, 0.29897, 0.29897],      
           [0.71025, 0.71025, 0.71025],      
           [0.1903, 0.1903, 0.1903],      
           [0.58964, 0.58964, 0.58964],      
           [0.48859, 0.48859, 0.48859],      
           [0.13195, 0.13195, 0.13195],      
           [0.40837, 0.40837, 0.40837],      
           [0.57558, 0.57558, 0.57558],      
           [0.50069, 0.50069, 0.50069],      
           [0.22966, 0.22966, 0.22966],      
           [0.5484, 0.5484, 0.5484],      
           [0.86903, 0.86903, 0.86903],      
           [0.055232, 0.055233, 0.055233],      
           [0.17579, 0.17579, 0.17579],      
           [0.75837, 0.75837, 0.75837],      
           [0.09151, 0.091513, 0.091513],      
           [0.25777, 0.25777, 0.25777],      
           [0.36152, 0.36152, 0.36152],      
           [0.91096, 0.91096, 0.91096],      
           [0.66104, 0.66104, 0.66104],      
           [0.1466, 0.1466, 0.1466],      
           [0.43244, 0.43244, 0.43244],      
           [0.33561, 0.33561, 0.33561],      
           [0.45605, 0.45605, 0.45605],      
           [0.28536, 0.28536, 0.28536],      
           [0.95476, 0.95476, 0.95476],      
           [0.79526, 0.79526, 0.79526],      
           [0.69348, 0.69348, 0.69348],      
           [0.63004, 0.63004, 0.63004],      
           [0.38367, 0.38367, 0.38367],      
           [0.11721, 0.11721, 0.11721],      
           [0.20109, 0.20109, 0.20109],      
           [0.30909, 0.30909, 0.30909],      
           [0.47964, 0.47964, 0.47964],      
           [0.52563, 0.52563, 0.52563],      
           [0.7274, 0.7274, 0.7274],      
           [0.60404, 0.60404, 0.60404],      
           [0.82892, 0.82892, 0.82892],      
           [0.54179, 0.54179, 0.54179],      
           [0.15393, 0.15393, 0.15393],      
           [0.12457, 0.12457, 0.12457],      
           [0.065886, 0.065887, 0.065887],      
           [0.23676, 0.23676, 0.23676],      
           [0.50682, 0.50682, 0.50682],      
           [0.80474, 0.80474, 0.80474],      
           [0.34215, 0.34215, 0.34215],      
           [0.2508, 0.2508, 0.2508],      
           [0.098827, 0.098829, 0.098828],      
           [0.92173, 0.92173, 0.92173],      
           [0.31577, 0.31577, 0.31577],      
           [0.5968, 0.5968, 0.5968],      
           [0.27853, 0.27853, 0.27853],      
           [0.32904, 0.32904, 0.32904],      
           [0.70181, 0.70181, 0.70181],      
           [0.22258, 0.22258, 0.22258],      
           [0.78588, 0.78588, 0.78588],      
           [0.16854, 0.16854, 0.16854],      
           [0.96599, 0.96599, 0.96599],      
           [0.6225, 0.6225, 0.6225],      
           [0.3899, 0.3899, 0.3899],      
           [0.83877, 0.83877, 0.83877],      
           [0.042788, 0.042788, 0.042788],      
           [0.68523, 0.68523, 0.68523],      
           [0.85883, 0.85883, 0.85883],      
           [0.55507, 0.55507, 0.55507],      
           [0.65315, 0.65315, 0.65315],      
           [0.45017, 0.45017, 0.45017],      
           [0.66902, 0.66902, 0.66902],      
           [0.083764, 0.083767, 0.083766],      
           [0.29219, 0.29219, 0.29219],      
           [0.20829, 0.20829, 0.20829],      
           [0.35511, 0.35511, 0.35511],      
           [0.9003, 0.9003, 0.9003]]      
      
grayCS_map = LinearSegmentedColormap.from_list('grayCS', cm_data)      
# For use of "viscm view"      
test_cm = grayCS_map      
      
if __name__ == "__main__":      
    import matplotlib.pyplot as plt      
    import numpy as np      
      
    try:      
        from viscm import viscm      
        viscm(grayCS_map)      
    except ImportError:      
        print("viscm not found, falling back on simple display")      
        plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',      
                   cmap=grayCS_map)      
    plt.show()      