import numpy as np

def source_name_finder(source_name):
    source_type_str = ["plane_IC", "square_IC", "square_s", 
                        "gaussian_IC", "MMS", "gaussian_s",
                        "gaussian_IC_2D", "line_source", "P1_su_olson_rad",
                        "P1_su_olson_mat", "P1_gaussian_rad", "P1_gaussian_mat", 
                        "P1_gaussian_rad_thick", "P1_gaussian_mat_thick",
                        "P1_su_olson_rad_thick", "P1_su_olson_mat_thick"]
    if source_name not in source_type_str:
        print('source name not here')
        print(source_name)
        assert(0)
    source_array = np.zeros(len(source_type_str))
    index = source_type_str.index(source_name)
    source_array[index] = 1
    return source_array