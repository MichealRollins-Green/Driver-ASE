function save_variables(variable_path)

cd(variable_path);

matlab_vars = 'MatLab_Variables';
somatic_calls = 'somatic_calls';
mutations = 'mutations';
matrix = 'matrix';
annotations = 'annotations';
mut_ase_auto = 'mut_ase_auto';
hits = 'hits';
driver_beds = 'Driver_Beds';

save('directory_variables','matlab_vars','somatic_calls','matrix','mutations','annotations','mut_ase_auto','hits','driver_beds');
