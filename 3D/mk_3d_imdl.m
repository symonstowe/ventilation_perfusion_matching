function imdl = mk_3d_imdl(shape_str)
	% Load the mat file if it exists - otherwise make an inverse model
	% input is the shape string required to get a shape from:
	% shape_library function in eidors
	if exist('imdl.mat','file')
        ld = load('imdl.mat','imdl');
        imdl = ld.imdl;
    else
        skip4 = {32,1,[0,5],[0,5],{'no_meas_current_next1'},1};

        elec_rotation = 0.03;
        n_elecs_per_plane = 16;
        spacing = 1;
        z_planes = [0.8,1.2];
        elec_pos = [n_elecs_per_plane,spacing+elec_rotation,z_planes];
        elec_shape = [0.05, [0, 0.02]];
        shape = shape_library('get',shape_str);
        fmdl = ng_mk_extruded_model({2,shape.boundary,1,0.15},elec_pos, ...
            elec_shape);
		% TODO - I have code to automate this somewhere in eidors...
        % Renumber the electrodes
        elec_order = [32,16,15,31,30,14,13,29, ...
                      28,12,11,27,26,10, 9,25, ...
                      24, 8, 7,23,22, 6, 5,21, ...
                      20, 4, 3,19,18, 2, 1,17];
        fmdl.electrode(:) = fmdl.electrode(elec_order);  

        [fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(skip4{:});

        vopt.imgsz = [32 32];
        vopt.square_pixels = true;
        vopt.save_memory = 1;
        opt.noise_figure = 1.0;
        vopt.zvec = linspace(0.6,1.4,10);
        % GREIT 3D with 2x16 electrode belt
        [imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);
        imdl3= mk_GREIT_model(imdl, 0.20, [], opt);
        imdl = imdl3;
        % Save the model to save time
        save('imdl.mat','imdl');
    end