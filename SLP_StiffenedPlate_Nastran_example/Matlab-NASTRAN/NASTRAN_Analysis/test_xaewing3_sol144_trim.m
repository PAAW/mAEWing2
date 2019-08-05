
sol144_bdf_base_fname = 'sol144_multiple_subcases_analysis';
sol144_bdf_fname = 'sol144_multiple';


dynamic_pressure_all = [0.160957  0.160957  0.058803 ];


flaps=[min_DRAG_Opt_flap';Optimal_Flap_WRBM';liftmax_Flap_rotations'];
load_factors = [1 3.8 max_load_factor];

[max_strain,max_disp]=xAEWing3_sol144_trim(base_bdf_folder,design_folder,...
    flaps,load_factors,dynamic_pressure_all,...
    sol144_bdf_base_fname,sol144_bdf_fname,...
    ref_chord,ref_area,...
    TE_flaps_num,Mach,...
    nastran,main_program_folder)