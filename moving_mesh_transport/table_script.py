from .solution_plotter import plot, make_tables_su_olson, make_tables_gaussian_thin, make_tables_thick_problems

def make_all_tables():
    # Square sources

    # su olson
    make_tables_su_olson(Ms=[3,4,5,6,6,10,10], N_spaces = [256, 256, 256, 256, 128, 32, 32], problem_name = 'su_olson', rad_or_transport = 'rad', 
                                c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['su_olson_phi.csv','su_olson_e.csv'], source_name_list = ['square_s'], uncollided = True, moving = True)

    # su olson s2
    make_tables_su_olson(Ms=[4,4,4,4,6,8,8], N_spaces = [128, 128, 128, 128, 128, 32, 32], problem_name = 'su_olson_s2', rad_or_transport = 'rad', 
                                c = 0.0, s2 = True, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['su_olson_s2_phi.csv','su_olson_s2_e.csv'], source_name_list = ['square_s'], uncollided = True, moving = True)

    # # const cv
    # make_tables_su_olson(Ms=[3,4,5,6,6,10,10], N_spaces = [256, 256, 256, 256, 128, 32, 32], problem_name = 'transfer_const_cv=0.03', rad_or_transport = 'rad', 
    #                             c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['cv=0.03_phi.csv','cv=0.03_e.csv'], source_name_list = ['square_s'], uncollided = True, moving = True)

    # # const cv s2
    # make_tables_su_olson(Ms=[4,4,4,4,8,8,8], N_spaces = [128, 128, 128, 128, 128, 64, 64], problem_name = 'transfer_const_cv=0.03_s2', rad_or_transport = 'rad', 
    #                             c = 0.0, s2 = True, cv0=0.03, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['cv=0.03_s2_phi.csv','cv=0.03_s2_e.csv'], source_name_list = ['square_s'], uncollided = True, moving = True)

    # # su olson gaussian
    # make_tables_gaussian_thin(Ms=[12,12,12,12,12,12,12], N_spaces = [64, 64, 64, 64, 64, 64, 64], problem_name = 'su_olson', rad_or_transport = 'rad', 
    #                             c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['su_olson_g_phi.csv','su_olson_g_e.csv'], source_name_list = ['gaussian_s'], uncollided = True, moving = False)

    # # su olson gaussian s2
    # make_tables_gaussian_thin(Ms=[12,12,12,12,12,8,12], N_spaces = [64, 64, 64, 64, 64, 64, 32], problem_name = 'su_olson_s2', rad_or_transport = 'rad', 
    #                             c = 0.0, s2 = True, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['su_olson_s2_g_phi.csv','su_olson_s2_g_e.csv'], source_name_list = ['gaussian_s'], uncollided = True, moving = False)

    # # const cv gauss

    # make_tables_gaussian_thin(Ms=[12,12,12,12,10,10,10], N_spaces =  [64,64,64,64,128,128,128], problem_name = 'transfer_const_cv=0.03', rad_or_transport = 'rad', 
    #                             c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['cv=0.03_g_phi.csv','cv=0.03_g_e.csv'], source_name_list = ['gaussian_s'], uncollided = True, moving = False)
    # # const cv gauss s2

    # make_tables_gaussian_thin(Ms=[12,12,12,12,12,12,12], N_spaces = [64,64,64,64,128,128,128], problem_name = 'transfer_const_cv=0.03_s2', rad_or_transport = 'rad', 
    #                         c = 0.0, s2 = True, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['cv=0.03_s2_g_phi.csv','cv=0.03_s2_g_e.csv'], source_name_list = ['gaussian_s'], uncollided = True, moving = False)

    # # thick problems

    # # thick su olson

    # make_tables_thick_problems(Ms=[10,10,10], N_spaces = [128,128,128], problem_name = 'su_olson_thick', rad_or_transport = 'rad', 
    #                         c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['su_olson_thick_phi.csv','su_olson_thick_e.csv'], source_name_list = ['square_s'], uncollided = False, moving = False)

    # make_tables_thick_problems(Ms=[10,10,10], N_spaces = [128,128,128], problem_name = 'su_olson_thick_s2', rad_or_transport = 'rad', 
    #                         c = 0.0, s2 = True, cv0=0.0, x0_or_sigma = 0.5, mat_or_rad ='rad', filenames = ['su_olson_thick_s2_phi.csv','su_olson_thick_s2_e.csv'], source_name_list = ['square_s'], uncollided = False, moving = False)

    # make_tables_thick_problems(Ms=[10,10,10], N_spaces = [128,128,128], problem_name = 'su_olson_thick', rad_or_transport = 'rad', 
    #                         c = 0.0, s2 = False, cv0=0.0, x0_or_sigma = 0.375, mat_or_rad ='rad', filenames = ['su_olson_thick_g_phi.csv','su_olson_thick_g_e.csv'], source_name_list = ['gaussian_s'], uncollided = False, moving = False)
    
    # make_tables_thick_problems(Ms=[10,10,10], N_spaces = [128,128,128], problem_name = 'su_olson_thick_s2', rad_or_transport = 'rad', 
    #                         c = 0.0, s2 = True, cv0=0.0, x0_or_sigma = 0.375, mat_or_rad ='rad', filenames = ['su_olson_thick_s2_g_phi.csv','su_olson_thick_s2_g_e.csv'], source_name_list = ['gaussian_s'], uncollided = False, moving = False)

    # make_tables_thick_problems(Ms=[10,10,10], N_spaces = [128,128,128], problem_name = 'transfer_const_cv=0.03_thick', rad_or_transport = 'rad', 
    #                         c = 0.0, s2 = False, cv0=0.03, x0_or_sigma = 0.375, mat_or_rad ='rad', filenames = ['cv=0.03_thick_g_phi.csv','cv=0.03_thick_g_e.csv'], source_name_list = ['gaussian_s'], uncollided = False, moving = False)
    
    # make_tables_thick_problems(Ms=[10,10,10], N_spaces = [128,128,128], problem_name = 'transfer_const_cv=0.03_thick_s2', rad_or_transport = 'rad', 
    #                         c = 0.0, s2 = True, cv0=0.03, x0_or_sigma = 0.375, mat_or_rad ='rad', filenames = ['cv=0.03_thick_s2_g_phi.csv','cv=0.03_thick_s2_g_e.csv'], source_name_list = ['gaussian_s'], uncollided = False, moving = False)