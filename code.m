
function [] = test()
    clc
    clear
    warning('off');
    
    time_step = 0.01;

    globalSys.x = [1 0 0];
    globalSys.y = [0 1 0];
    globalSys.z = [0 0 1];
    %% get angle
    joint_angle = get_joint_angle();
    plot_joint_angle(joint_angle, get_toeoff_percentage(262,316,351), get_toeoff_percentage(218,270,308));

    %% get omega / angular velocity
    joints = {'Lankle', 'Lknee', 'Lhip', 'Rankle', 'Rknee', 'Rhip'};
    for i = 1:numel(joints)
        angle_tmp = joint_angle.([joints{i} '_degrees']);
        joint_angle.([joints{i} '_omega']) = get_derivative(deg2rad(angle_tmp), time_step);
        joint_angle.([joints{i} '_acc']) = get_2th_derivative(deg2rad(angle_tmp), time_step);
    end
    
    %% get accerlation
    get_coordinates();      % generate Bone_coordinates and Bone_center
    coordinates = load('Bone_coordinates.mat');
    center = load('Bone_center.mat');
	ends = load('Bone_end.mat');
    
    prefixes = {'LTO', 'LTI', 'LFE', 'RFE', 'RTI', 'RTO'};
    cycle = {'l', 'l', 'l', 'r', 'r', 'r'};
    absolute_acc = struct();
    % 
    for i = 1:numel(prefixes)
        prefixe = prefixes{i};
        
        bone_cor = coordinates.(prefixe);
        angle = get_absolute_angle(bone_cor, cycle{i});
        absolute_acc.(prefixe).('omega') = get_derivative(deg2rad(angle), time_step);
        absolute_acc.(prefixe).('alphaacc') = get_2th_derivative(deg2rad(angle), time_step);

        bone_center = center.(prefixe);
        absolute_acc.(prefixe).('acc') = get_2th_derivative(bone_center, time_step);
    end

    %% get ground force
    [GroundForce_pos, GroundForce] = get_groundforce();
    
    %% get Moment
    segments = {'LTO', 'RTO', 'LTI', 'RTI', 'LFE', 'RFE'};  % bone
    segement_ends = {'LTOO', 'LAO', 'RTOO', 'RAO', 'LAO', 'LKO', 'RAO', 'RKO', 'LKO', 'LHO', 'RKO', 'RHO'};                       % end points of bone
    Joints = {'Lankle', 'Rankle', 'Lknee', 'Rknee', 'Lhip', 'Rhip'};

	% force value, force position
	GroundForces = {'LF', 'RF', 'LF', 'RF', 'LF', 'RF'};        % from GroundForce.LF
	GroundForce_poses = {'LF', 'RF', 'LF', 'RF', 'LF', 'RF'};   % from GroundForce_pos.LF
	
	WeightRatio = [0.0145 0.0145 0.0465 0.0465 0.1 0.1];
	weight = 71.5;  % Kg
    g = 9.81;       % m/s^2
	C_of_G = [0.475 0.475 0.302 0.302 0.323 0.323];

    Force_joint = struct();
	
	for i = 1:numel(Joints)
        segment = segments{i};
        joint = Joints{i};

        % local segment variables
        % Omega_segment = joint_angle.([Joints{i} '_omega']);
        % Acc_segment = joint_angle.([Joints{i} '_acc']);
        Omega_segment = absolute_acc.(segment).omega;
        Acc_segment = absolute_acc.(segment).alphaacc;

        coordinate_segment.x = coordinates.(segment).x_cycle;
        coordinate_segment.y = coordinates.(segment).y_cycle;
        coordinate_segment.z = coordinates.(segment).z_cycle;

        % global variables
	    CoM 		= center.(segment);
		Distal_point 	= ends.(segement_ends{2*i-1});
		Proximal_point 	= ends.(segement_ends{2*i});
		Ground_pos = GroundForce_pos.(GroundForce_poses{i});
		
		length_bone = norm(Proximal_point(1,:) - Distal_point(1,:));
		% moment arm vector
        V_GroundF 	= CoM - GroundForce_pos.(GroundForces{i});
		V_Proximal  = CoM - Proximal_point;
        V_Distal 	= CoM - Distal_point;
		
		% force vector
		F_Ground = GroundForce.(GroundForces{i});
		Sw_weight = WeightRatio(i)*weight;
		F_Weight = [0 0 -Sw_weight*g];
		
		R_of_gyration = C_of_G(i)*length_bone;
        MoI_x = Sw_weight * R_of_gyration^2;
        MoI_y = MoI_x;
        MoI_z = 0;

		Moment_inertia_x = - (MoI_x.*Acc_segment(:,1) + (MoI_z-MoI_y).*Omega_segment(:,2).*Omega_segment(:,3));
		Moment_inertia_y = - (MoI_y.*Acc_segment(:,2) + (MoI_x-MoI_z).*Omega_segment(:,1).*Omega_segment(:,3));
        Moment_inertia_z = - (MoI_z.*Acc_segment(:,3) + (MoI_y-MoI_x).*Omega_segment(:,1).*Omega_segment(:,2));

        Moment_inertia = [Moment_inertia_x Moment_inertia_y Moment_inertia_z];

        

		switch(joint)
			case({'Lankle', 'Rankle'})
			% ReactionForce + F_Ground + F_Weight = m*a
				Lankle_RF =	-F_Ground -F_Weight + Sw_weight*absolute_acc.(segment).acc;
                Lankle_RF_localSegment = transformVector(Lankle_RF, globalSys, coordinate_segment);
                Force_joint.(joint).force = Lankle_RF;
                Force_joint.(joint).localforce = Lankle_RF_localSegment;

			% Moment at XYZ axis  direction is the same as global
			% Moment_ground + Moment_Joint + Moment_joint_force + Moment_inertia = 0
% 			Lankle_M = -( cross(V_Proximal,Lankle_RF) + cross(V_GroundF,F_Ground) + Moment_inertia(:,1));
%                 Lankle_Mx = Lankle_M(:,1)
                M1 = cross(V_Proximal,Lankle_RF);
                M1_localSegment = transformVector(M1, globalSys, coordinate_segment);
                M2 = cross(V_GroundF,F_Ground);
                M2_localSegment = transformVector(M2, globalSys, coordinate_segment);

                Lankle_M_segment = -( M1_localSegment + M2_localSegment + Moment_inertia);
				Lankle_M = transformVectorBack(Lankle_M_segment, globalSys, coordinate_segment);
				Force_joint.(joint).M = Lankle_M;
                Force_joint.(joint).localM = Lankle_M_segment;

			case('Lknee')
			% knee_RF + -ankle_RF + F_Weight = m*a
				Lknee_RF =	 -F_Weight --Force_joint.('Lankle').force + Sw_weight*absolute_acc.(segment).acc;
                Lknee_RF_localSegment = transformVector(Lknee_RF, globalSys, coordinate_segment)
                Force_joint.(joint).force = Lknee_RF;
                Force_joint.(joint).localforce = Lknee_RF_localSegment;
				
				% Moment at XYZ axis  direction is the same as global
				% Moment_ankle+Moment_ankleRF+Moment_kneeRF + Moment_knee + Moment_inertia = 0
                M1 = cross(V_Distal	, - Force_joint.('Lankle').force);
                M1_local = transformVector(M1, globalSys, coordinate_segment);
                M2 = cross(V_Proximal ,Lknee_RF);
                M2_local = transformVector(M2, globalSys, coordinate_segment);
				Lankle_M_local = transformVector(Force_joint.('Lankle').M, globalSys, coordinate_segment);

                Lknee_M_segment = -( M1_local + M2_local + Moment_inertia-Lankle_M_local);
				Lknee_M = transformVectorBack(Lknee_M_segment, globalSys, coordinate_segment);
				Force_joint.(joint).M = Lknee_M;
                Force_joint.(joint).localM = Lknee_M_segment;

			case('Rknee')
			% knee_RF + -ankle_RF + F_Weight = m*a
				Rknee_RF = -F_Weight --Force_joint.('Rankle').force + Sw_weight*absolute_acc.(segment).acc;
                Rknee_RF_localSegment = transformVector(Rknee_RF, globalSys, coordinate_segment)
                Force_joint.(joint).force = Rknee_RF;
                Force_joint.(joint).localforce = Rknee_RF_localSegment;
				
				% Moment at XYZ axis  direction is the same as global
				% Moment_ankle+Moment_ankleRF+Moment_kneeRF + Moment_knee + Moment_inertia = 0
                M1 = cross(V_Distal	, - Force_joint.('Rankle').force);
                M1_local = transformVector(M1, globalSys, coordinate_segment);
                M2 = cross(V_Proximal ,Rknee_RF);
                M2_local = transformVector(M2, globalSys, coordinate_segment);
				Rankle_M_local = transformVector(Force_joint.('Rankle').M, globalSys, coordinate_segment);

                Rknee_M_segment = -( M1_local + M2_local + Moment_inertia-Rankle_M_local);
				Rknee_M = transformVectorBack(Rknee_M_segment, globalSys, coordinate_segment);
				Force_joint.(joint).M = Rknee_M;
                Force_joint.(joint).localM = Rknee_M_segment;
				
			case('Lhip')
			% hip_RF + -knee_RF + F_Weight = m*a
				Lhip_RF = -F_Weight	--Force_joint.('Lknee').force + Sw_weight*absolute_acc.(segment).acc;
                Lhip_RF_localSegment = transformVector(Lhip_RF, globalSys, coordinate_segment)
                Force_joint.(joint).force = Lhip_RF;
                Force_joint.(joint).localforce = Lhip_RF_localSegment;
				
				% Moment at XYZ axis  direction is the same as global
				% Moment_knee + Moment_kneeRF + Moment_hipRF + Moment_hip + Moment_inertia = 0
                M1 = cross(V_Distal	, - Force_joint.('Lknee').force);
                M1_local = transformVector(M1, globalSys, coordinate_segment);
                M2 = cross(V_Proximal ,Lhip_RF);
                M2_local = transformVector(M2, globalSys, coordinate_segment);
				Lknee_M_local = transformVector(Force_joint.('Lknee').M, globalSys, coordinate_segment);

                Lhip_M_segment = -( M1_local + M2_local + Moment_inertia-Lknee_M_local);
				Lhip_M = transformVectorBack(Lhip_M_segment, globalSys, coordinate_segment);
				Force_joint.(joint).M = Lhip_M;
                Force_joint.(joint).localM = Lhip_M_segment;
			
			case('Rhip')
			% hip_RF + -knee_RF + F_Weight = m*a
				Rhip_RF = -F_Weight	--Force_joint.('Rknee').force + Sw_weight*absolute_acc.(segment).acc;
                Rhip_RF_localSegment = transformVector(Rhip_RF, globalSys, coordinate_segment)
                Force_joint.(joint).force = Rhip_RF;
                Force_joint.(joint).localforce = Rhip_RF_localSegment;
				
				% Moment at XYZ axis  direction is the same as global
				% Moment_knee + Moment_kneeRF + Moment_hipRF + Moment_hip + Moment_inertia = 0
                M1 = cross(V_Distal	, - Force_joint.('Rknee').force);
                M1_local = transformVector(M1, globalSys, coordinate_segment);
                M2 = cross(V_Proximal ,Rhip_RF);
                M2_local = transformVector(M2, globalSys, coordinate_segment);
				Rknee_M_local = transformVector(Force_joint.('Rknee').M, globalSys, coordinate_segment);

                Rhip_M_segment = -( M1_local + M2_local + Moment_inertia-Rknee_M_local);
				Rhip_M = transformVectorBack(Rhip_M_segment, globalSys, coordinate_segment);
				Force_joint.(joint).M = Rhip_M;
                Force_joint.(joint).localM = Rhip_M_segment;
				
		end
 

%         
%         plot(normalized_gaitcycle,R_M_Ankle/w,'DisplayName','Ankle Moment')
    end

    plot_joint_moment(Force_joint, get_toeoff_percentage(262,316,351), get_toeoff_percentage(218,270,308),weight)
    save('Force_joint.mat', 'Force_joint');

	%% get power
    Power_joint = struct();

	for i = 1:numel(Joints)
		joint = Joints{i};
		Power_joint.(joint).power = -Force_joint.(joint).localM .*  joint_angle.([joint '_omega']);
		Power_joint.(joint).sumpower = sum( Power_joint.(joint).power, 2 );
		
    end
    
    plot_joint_power(Power_joint, get_toeoff_percentage(262,316,351), get_toeoff_percentage(218,270,308),weight)
    save('Power_joint.mat', 'Power_joint');
end


function angle = get_joint_angle()

    get_coordinates();
    load('Bone_coordinates.mat');

    % Cal angle of trunk and pelvis. Note: they are refer to global coordinate

    joint_center = {'PEL', 'TRX'};
    joint_upper = {'pelvis', 'thorax'};

    for prefix_idx = 1:(numel(joint_upper))
        jointName = joint_upper{prefix_idx};
        prefix = joint_center{prefix_idx};

        % loop all time step
        eval( [ jointName '_degrees = [];' ] )
        for i = 1:length(eval([prefix '.x']))
            stri = num2str(i);
            x1 = eval([prefix '.x(' stri ',:);']);
            y1 = eval([prefix '.y(' stri ',:);']);
            z1 = eval([prefix '.z(' stri ',:);']);

            x2 = [1,0,0];
            y2 = [0,1,0];
            z2 = [0,0,1];

            tmp = compute_euler_angles([x2;y2;z2], [x1;y1;z1]);
            eval( [ jointName '_degrees = [' jointName '_degrees;tmp];' ] )
        end
    end


    % Cal angle of other joints
    prefixes = {'LTO', 'LTI', 'LFE', 'PEL', 'RFE', 'RTI', 'RTO', 'TRX'};
    SIGN_direction = [-1, -1, -1, +1, +1, +1];
    joint = {'Lankle', 'Lknee', 'Lhip', 'Rhip', 'Rknee', 'Rankle'};

    for prefix_idx = 1:(numel(prefixes)-2)
        jointName = joint{prefix_idx};
        prefix = prefixes{prefix_idx};
        next_prefix = prefixes{prefix_idx + 1};

        SIGN = SIGN_direction(prefix_idx);

        % loop all time step
        eval( [ jointName '_degrees = [];' ] )
        for i = 1:length(eval([prefix '.x']))
            stri = num2str(i);
            x1 = eval([prefix '.x(' stri ',:);']);
            y1 = eval([prefix '.y(' stri ',:);']);
            z1 = eval([prefix '.z(' stri ',:);']);

            x2 = eval([next_prefix '.x(' stri ',:);']);
            y2 = eval([next_prefix '.y(' stri ',:);']);
            z2 = eval([next_prefix '.z(' stri ',:);']);

            globalx = [1,0,0];
            globaly = [0,1,0];
            globalz = [0,0,1];

            tmp1 = compute_euler_angles([globalx; globaly; globalz], [x1;y1;z1]);
            tmp2 = compute_euler_angles([globalx; globaly; globalz], [x2;y2;z2]);
            eval( [ jointName '_degrees = [' jointName '_degrees; (tmp2 - tmp1) .* SIGN];' ] )
        end
    end

    % left leg cycle
    L_index = 262:351;
    Lankle_degrees = Lankle_degrees(L_index,:);
    Lhip_degrees = Lhip_degrees(L_index,:);
    Lknee_degrees = Lknee_degrees(L_index,:);
    Lpelvis_degrees = pelvis_degrees(L_index,:);
    Lthorax_degrees = thorax_degrees(L_index,:);

    % right leg cycle
    R_index = 218:308;
    Rankle_degrees = Rankle_degrees(R_index,:);
    Rhip_degrees = Rhip_degrees(R_index,:);
    Rknee_degrees = Rknee_degrees(R_index,:);
    Rpelvis_degrees = pelvis_degrees(R_index,:);
    Rthorax_degrees = thorax_degrees(R_index,:);

    % save degree to mat
    % Define the list of variables to save
    variables_to_save = {'Lankle_degrees', 'Lhip_degrees', 'Lknee_degrees', ...
        'Lpelvis_degrees', 'Lthorax_degrees', ...
        'Rankle_degrees', 'Rhip_degrees', 'Rknee_degrees', ...
        'Rpelvis_degrees', 'Rthorax_degrees'};

    % Save the specified variables to a .mat file
    save('joint_variables.mat', variables_to_save{:});
    
    angle = load('joint_variables.mat');
end


% output unit rad
function angle = get_absolute_angle(coordinates, cycle)

    L_index = 262:351;
    R_index = 218:308;

    angle = [];

    if cycle == 'r'
        index = R_index;
    elseif cycle == 'l'
        index = L_index;
    end


    % Cal angle refer to global coordinate
    for i = index
        x1 = coordinates.x(i,:);
        y1 = coordinates.y(i,:);
        z1 = coordinates.z(i,:);
        
        x2 = [1,0,0];
        y2 = [0,1,0];
        z2 = [0,0,1];

        angle(end+1,:) = compute_euler_angles([x2;y2;z2], [x1;y1;z1]);
    end
%     save('absolute_angle.mat', variable_names{:});

end

function [ax, ay, az, alpha_x, alpha_y, alpha_z] = get_bone_acceleration()

    Bones  = load('Bone_coordinates.mat');
    
    % get trasnlation accelaration

    
    
    
    % get rotation acceleration
    
end

function [] = get_coordinates()

    % Load data from a text file into a variable
    data = readtable('Samarth_gaitnorm18.txt');
    variable_names = data.Properties.VariableNames;

    % x y z coordinate system please refer to description.pdf Figure 1.

    % Define the strings 'LTI' and 'LTO'
    prefixes = {'LTO', 'LTI', 'LFE', 'PEL', 'RFE', 'RTI', 'RTO', 'TRX'};
    cycle = {'l','l','l','x','r','r','r','x'};

    for prefix_idx = 1:numel(prefixes)
        prefix = prefixes{prefix_idx};

        % find index of variable_name in structure data
        indices_with_prefix = find(~cellfun('isempty', strfind(variable_names, prefix)));

        % Extract and convert mm into m
        eval([prefix 'O' '= data{:, indices_with_prefix(1:3)}./1000;']);
        eval([prefix 'P' '= data{:, indices_with_prefix(4:6)}./1000;']);
        eval([prefix 'L' '= data{:, indices_with_prefix(7:9)}./1000;']);
        eval([prefix 'A' '= data{:, indices_with_prefix(10:12)}./1000;']);

        % Perform calculations for each dimension (x, y, z)
        eval([prefix '_x = ' prefix 'L - ' prefix 'O;']);
        eval([prefix '_y = ' prefix 'A - ' prefix 'O;']);
        eval([prefix '_z = ' prefix 'P - ' prefix 'O;']);

    end

%     % foot case is special
    LTO_x = - (LTOL - LTOO);
    LTO_y = - (LTOP - LTOO);
    LTO_z = LTOA - LTOO;
% 
    RTO_x = - (RTOL - RTOO);
    RTO_y = - (RTOP - RTOO);
    RTO_z = RTOA - RTOO;
% 
%     % adjust direction
    RTI_x = -RTI_x;
    LTI_x = -LTI_x;

    RFE_x = -RFE_x;
    LFE_x = -LFE_x;

    TRX_z = -TRX_z;
    PEL_x = -PEL_x;
%     

%     timeI = 200
%     scatter3(LTOO(timeI,1),LTOO(timeI,2),LTOO(timeI,3)); hold on;
%     scatter3(LTOP(timeI,1),LTOP(timeI,2),LTOP(timeI,3))
%     scatter3(LTOL(timeI,1),LTOL(timeI,2),LTOL(timeI,3))
%     scatter3(LTOA(timeI,1),LTOA(timeI,2),LTOA(timeI,3))
% 
%     scatter3(RTOO(timeI,1),RTOO(timeI,2),RTOO(timeI,3)); hold on;
%     scatter3(RTOP(timeI,1),RTOP(timeI,2),RTOP(timeI,3))
%     scatter3(RTOL(timeI,1),RTOL(timeI,2),RTOL(timeI,3))
%     scatter3(RTOA(timeI,1),RTOA(timeI,2),RTOA(timeI,3))

   
    % save degree to mat
    for i = 1:numel(prefixes)
        prefix = prefixes{i};
        if cycle{i} == 'l'
            index = '(262:351, :)';
            eval([prefix '.x_cycle = ' prefix '_x' index ' ;']);
            eval([prefix '.y_cycle = ' prefix '_y' index ' ;']);
            eval([prefix '.z_cycle = ' prefix '_z' index ' ;']);
        elseif cycle{i} == 'r'
            index = '(218:308, :)';
            eval([prefix '.x_cycle = ' prefix '_x' index ' ;']);
            eval([prefix '.y_cycle = ' prefix '_y' index ' ;']);
            eval([prefix '.z_cycle = ' prefix '_z' index ' ;']);
        end

        eval([prefix '.x = ' prefix '_x;']);
        eval([prefix '.y = ' prefix '_y;']);
        eval([prefix '.z = ' prefix '_z;']);
    end
    
    % Save the specified variables to a .mat file
    save('Bone_coordinates.mat', prefixes{:});
    
    % left leg cycle
    L_index = 262:351;
    LTOO = LTOO(L_index, :);
    LTIO = LTIO(L_index, :);
    LFEO = LFEO(L_index, :);
    LFEP = LFEP(L_index, :);

    % right leg cycle
    R_index = 218:308;
    RTOO = RTOO(R_index, :);
    RTIO = RTIO(R_index, :);
    RFEO = RFEO(R_index, :);
    RFEP = RFEP(R_index, :);

    % bone center of mass
    LTO = (LTOO+LTIO)./2;       % foot
    LTI = 0.433*LFEO+0.567*LTIO;
    LFE = 0.567*LFEO+0.433*LFEP;
    RTO = (RTOO+RTIO)./2;
    RTI = 0.433*RFEO+0.567*RTIO;
    RFE = 0.567*RFEO+0.433*RFEP;

    variable_names = {'LTO', 'LTI', 'LFE', 'RTO', 'RTI', 'RFE'};
    % Save the specified variables to a .mat file
    save('Bone_center.mat', variable_names{:});

    LAO = LTIO;
    LKO = LFEO;
    LHO = LFEP;
    RAO = RTIO;
    RKO = RFEO;
    RHO = RFEP;
    variable_names = {'LTOO', 'LAO', 'LKO', 'LHO', 'RTOO', 'RAO', 'RKO', 'RHO'};
    save('Bone_end.mat', variable_names{:});

end

function [pos, force] = get_groundforce()

    data = readtable('walking_FP.csv');
    variable_names = data.Properties.VariableNames;

    % x y z coordinate system please refer to description.pdf Figure 1.
    GroundForce = struct(); % LF xyz
	GroundForce_pos = struct();

    %% LEFT strike plate
    L_index = (262:351) * 10;

    prefixes_L2_pos = 'FP1_COP';
    prefixes_L2_force = 'FP1_Force';

    % find index of variable_name in structure data
%     indices_with_prefix = find(~cellfun('isempty', strfind(variable_names, prefixes_L2_pos)));
    indices_with_prefix = contains(variable_names, prefixes_L2_pos);
    GroundForce_pos.LF = data{L_index, indices_with_prefix}./1000;

    % find index of variable_name in structure data
    indices_with_prefix = contains(variable_names, prefixes_L2_force);
    GroundForce.LF = data{L_index, indices_with_prefix};


    %% RIGHT strike plate
    % right leg cycle
    R_index = (218:308) * 10;

    prefixes_R2_pos = 'FP2_COP';
    prefixes_R2_force = 'FP2_Force';

    % find index of variable_name in structure data
    indices_with_prefix = contains(variable_names, prefixes_R2_pos);
    GroundForce_pos.RF = data{R_index, indices_with_prefix}./1000;

    % find index of variable_name in structure data
    indices_with_prefix = contains(variable_names, prefixes_R2_force);
    GroundForce.RF = data{R_index, indices_with_prefix};

    pos = GroundForce_pos;
    force = GroundForce;

end

function [derivation] = get_2th_derivative(data, time_step)
    % Calculate the second derivative of the given data
    % with respect to the provided time step.
    % d^2(x+1) = {data(x) - 2data(x+1) - data(x+2)} / (dt)^2
    % data can be nxm

    % Compute the second derivative using central differences
    derivation = diff(data, 2) / time_step^2;
    
    insert = nan(1,size(derivation,2));
    
    derivation = [insert; derivation; insert];
end


function [derivation] = get_derivative(data, time_step)
    % Calculate the derivative of the given data
    % with respect to the provided time step.
    % d(x) = {data(x+1) - data(x)} / dt

    % Compute the second derivative using central differences
    derivation = diff(data, 1) / time_step;
    
    insert = nan(1,size(derivation,2));
    
    derivation = [derivation; insert];
end

function eul_degrees = compute_euler_angles(sys1, sys2)
    % Angle from sys1 to sys2
    
    % Define the direction vectors for each coordinate system
    % Let's say you have two sets of direction vectors: dir1 and dir2
    % Each set of direction vectors should be a 3x3 matrix where each row represents a unit vector (x, y, z)
    
    % Example direction vectors for the first coordinate system
    % sys1 (u,v,w) = [u1_x u1_y u1_z;
    %         v1_x v1_y v1_z;
    %         w1_x w1_y w1_z];
    
    % Example direction vectors for the second coordinate system
    % sys2 (u,v,w) = [u2_x u2_y u2_z;
    %         v2_x v2_y v2_z;
    %         w2_x w2_y w2_z];
    
    % sys1 = [1 0 0;
    %         0 1 0;
    %         0 0 1];
    % 
    % sys2 = [1 0 0;
    %         0 -1 0;
    %         0 0 -1];
    % 
    % compute_euler_angles(sys1, sys2)
    
    %  each row represents a unit vector (x, y, z)
    for i = 1:3
       sys1(i,:) =  sys1(i,:)./norm(sys1(i,:) );
       sys2(i,:) =  sys2(i,:)./norm(sys2(i,:) );
    end
    
    % Compute the rotation matrix from sys1 to sys2
    R = sys2' * (sys1')^-1;
    
    % Extract Euler angles
    eul = rotm2eul(R, 'XYZ');
    
    % Convert Euler angles from radians to degrees
    eul_degrees = rad2deg(eul);

end

function newVector = transformVector(v, sys1, sys2)
    % v is a vector in sys1 ，output vector v2 in a series of sys2.
	% sys1 and sys2 define the direction vectors for each coordinate system
	% Each set of direction vectors should be a 3x3 matrix where each row represents a unit vector (x, y, z)

	% Example direction vectors for the first coordinate system
	% sys1 u = sys1.x; v = sys1.y; w = sys1.z

	% Example direction vectors for the second coordinate system
	% sys2 u = sys2.x; v = sys2.y; w = sys2.z
    
    newVector = [];

    for j = 1: length(sys2.x)
        sys11 = [sys1.x; sys1.y; sys1.z];
        sys22 = [sys2.x(j,:);sys2.y(j,:);sys2.z(j,:)];
        
        %  each row represents a unit vector (x, y, z)
        for i = 1:3
           sys11(i,:) =  sys11(i,:)./norm(sys11(i,:) );
           sys22(i,:) =  sys22(i,:)./norm(sys22(i,:) );
        end
    
	    % Compute the rotation matrix from sys1 to sys2
 	    R = sys22' * (sys11')^-1;
    
        % 计算变换后的向量
        v2 = R' * v(j,:)';
        
        % each row is x,y,z coordinate
        v2 = v2';

        newVector(j,:) = v2;
    end
end

function newVector = transformVectorBack(v, sys1, sys2)

    % v is a vector in sys1 ，each row is a force
	% output newVector in a series of coordinates based on sys2.x, sys2.y, sys2.z
	% sys1 and sys2 define the direction vectors for each coordinate system

	%  direction vectors for the first coordinate system
	% sys1 u = sys1.x; v = sys1.y; w = sys1.z

	%  direction vectors for the second coordinate system in each row
	% sys2 u = sys2.x; v = sys2.y; w = sys2.z
    
    newVector = [];

    for j = 1: length(sys2.x)
        sys11 = [sys1.x; sys1.y; sys1.z];
        sys22 = [sys2.x(j,:);sys2.y(j,:);sys2.z(j,:)];
        
        %  each row represents a unit vector (x, y, z)
        for i = 1:3
           sys11(i,:) =  sys11(i,:)./norm(sys11(i,:) );
           sys22(i,:) =  sys22(i,:)./norm(sys22(i,:) );
        end
    
	    % Compute the rotation matrix from sys2 to sys1
 	    R = sys11' * (sys22')^-1;
    
        % 计算变换后的向量
        v2 = R' * v(j,:)';
        
        % each row is x,y,z coordinate
        v2 = v2';

        newVector(j,:) = v2;
    end
end

function percentage = get_percentage(One_cycle_sample_index)

    samples = length(One_cycle_sample_index);
    percentage = linspace(0, 100, samples);
end

function plot_joint_angle(data_structure, left_toeoff_event, right_toeoff_event)

    font = 14;

    % Plotting all data in the structure
    fields = fieldnames(data_structure); % Get the field names of the structure
    
    for i = 1:numel(fields)
        field = strrep(fields{i}, '_', ' ');
        switch(field)
            case('Rhip degrees')
                figure(1)
                subplot(3,3,7)
                hold on
                data = data_structure.(fields{i})(:,1);
                plot(get_percentage(data), data, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ flexion/- extension', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');

                subplot(3,3,8)
                hold on
                data = data_structure.(fields{i})(:,2);
                plot(get_percentage(data), data, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ addcution/- abdcution', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');

                title('Hip Joint Degree in 3D', FontSize=font)
                
                subplot(3,3,9)
                hold on
                data = data_structure.(fields{i})(:,3);
                plot(get_percentage(data), data, '-g', LineWidth=2, DisplayName='right');
                ylabel('- external/+ internal rotation', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                

%                 legend({ '+ flexion/- extension', '+ addcution/- abdcution', '+ internal/- external rotation'}, FontSize=font)

            case('Lhip degrees')
                figure(1)
                subplot(3,3,7)
                hold on
                data = data_structure.(fields{i})(:,1);
                plot(get_percentage(data), data, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ flexion/- extension', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-20 70])

                subplot(3,3,8)
                hold on
                data = data_structure.(fields{i})(:,2);
                plot(get_percentage(data), -data, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ addcution/- abdcution', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-20 70])
                
                subplot(3,3,9)
                hold on
                data = data_structure.(fields{i})(:,3);
                plot(get_percentage(data), -data, '-r', LineWidth=2, DisplayName='left');
                ylabel('- external/+ internal rotation', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-20 70])
%                 legend({ '+ flexion/- extension', '+ abduction/- addcution', '+ external/- internal rotation'}, FontSize=font)

            case('Rankle degrees')
                figure(1)
                subplot(3,3,1)
                data = data_structure.(fields{i})(:,1);
                plot(get_percentage(data), data, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ dorsi/- plantarflexion', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                

                subplot(3,3,2)
                data = data_structure.(fields{i})(:,2);
                plot(get_percentage(data), -data, '-g', LineWidth=2, DisplayName='right');
                ylabel('- inversion/+ eversion', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');

                title('Ankle Joint Degree in 3D', FontSize=font)

                subplot(3,3,3)
                data = data_structure.(fields{i})(:,3);
                plot(get_percentage(data), data, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ internal/- external rotation', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                

%                 legend(FontSize=font)

            case('Lankle degrees')
                figure(1)
                subplot(3,3,1)
                hold on
                data = data_structure.(fields{i})(:,1);
                plot(get_percentage(data), data, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ dorsi/- plantarflexion', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-30 30])

                subplot(3,3,2)
                hold on
                data = data_structure.(fields{i})(:,2);
                plot(get_percentage(data), data, '-r', LineWidth=2, DisplayName='left');
                ylabel('- inversion/+ eversion', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-30 30])
                
                subplot(3,3,3)
                hold on
                data = data_structure.(fields{i})(:,3);
                plot(get_percentage(data), -data, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ internal/- external rotation', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-30 30])
%                 legend(FontSize=font)

            case('Rknee degrees')
                figure(1)
                subplot(3,3,4)
                data = data_structure.(fields{i})(:,1);
                plot(get_percentage(data), -data, '-g', LineWidth=2, DisplayName='right');
                ylabel('- extension/+ flexion', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                

                subplot(3,3,5)
                data = data_structure.(fields{i})(:,2);
                plot(get_percentage(data), data, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ varus/- valgus', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                

                title('Knee Joint Degree in 3D', FontSize=font)

                subplot(3,3,6)
                data = data_structure.(fields{i})(:,3);
                plot(get_percentage(data), data, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ internal/- external rotation', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                

%                 legend({ '+ extension/- flexion', '+ addcution/- abdcution', '+ internal/- external rotation'}, FontSize=font)
            case('Lknee degrees')
                figure(1)
                subplot(3,3,4)
                hold on
                data = data_structure.(fields{i})(:,1);
                plot(get_percentage(data), -data, '-r', LineWidth=2, DisplayName='left');
                ylabel('- extension/+ flexion', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-25 60])

                subplot(3,3,5)
                hold on
                data = data_structure.(fields{i})(:,2);
                plot(get_percentage(data), -data, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ varus/- valgus', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-25 60])
                
                subplot(3,3,6)
                hold on
                data = data_structure.(fields{i})(:,3);
                plot(get_percentage(data), -data, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ internal/- external rotation', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-25 60])
%                 legend({ '+ extension/- flexion', '+ abduction/- addcution', '+ external/- internal rotation'}, FontSize=font)
            
            case{'Lpelvis degrees'}
                figure(2)
                subplot(2,3,1)
                hold on
                data = data_structure.(fields{i})(:,1);
                plot(get_percentage(data), -data, '-r', LineWidth=2, DisplayName='left');
                ylabel('- posterior/+ anterior', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-30 30])

                subplot(2,3,2)
                hold on
                data = data_structure.(fields{i})(:,2);
                plot(get_percentage(data), data, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ bend right/- left', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-30 30])

                title('Pelvis Joint Degree in 3D', FontSize=font)
                
                subplot(2,3,3)
                hold on
                data = data_structure.(fields{i})(:,3);
                plot(get_percentage(data), -data, '-r', LineWidth=2, DisplayName='left');
                ylabel('- left/+ right rotation', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-30 30])

%                 legend({ '+ posterior/- anterior', '+ bend right/- left', '+ left/- right rotation'}, FontSize=font)
            case{'Rpelvis degrees'}
                figure(2)
                subplot(2,3,1)
                hold on
                data = data_structure.(fields{i})(:,1);
                plot(get_percentage(data), -data, '-g', LineWidth=2, DisplayName='right');
                ylabel('- posterior/+ anterior', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');

                subplot(2,3,2)
                hold on
                data = data_structure.(fields{i})(:,2);
                plot(get_percentage(data), -data, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ bend right/- left', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                
                subplot(2,3,3)
                hold on
                data = data_structure.(fields{i})(:,3);
                plot(get_percentage(data), data, '-g', LineWidth=2, DisplayName='right');
                ylabel('- left/+ right rotation', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');

            case{'Lthorax degrees'}
                figure(2)
                subplot(2,3,4)
                hold on
                data = data_structure.(fields{i})(:,1);
                plot(get_percentage(data), -data, '-r', LineWidth=2, DisplayName='left');
                ylabel('- posterior/+ anterior', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-10 10])

                subplot(2,3,5)
                hold on
                data = data_structure.(fields{i})(:,2);
                plot(get_percentage(data), -data, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ bend right/- left', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-10 10])

                title('Thorax Joint Degree in 3D', FontSize=font)
                
                subplot(2,3,6)
                hold on
                data = data_structure.(fields{i})(:,3);
                plot(get_percentage(data), -data, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ left/- right rotation', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-10 10])

%                 legend({ '+ posterior/- anterior', '+ bend right/- left', '+ left/- right rotation'}, FontSize=font)
            case{'Rthorax degrees'}
                figure(2)
                subplot(2,3,4)
                hold on
                data = data_structure.(fields{i})(:,1);
                plot(get_percentage(data), -data, '-g', LineWidth=2, DisplayName='right');
                ylabel('- posterior/+ anterior', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                

                subplot(2,3,5)
                hold on
                data = data_structure.(fields{i})(:,2);
                plot(get_percentage(data), data, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ bend right/- left', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');

                title('Thorax Joint Degree in 3D', FontSize=font)
                
                subplot(2,3,6)
                hold on
                data = data_structure.(fields{i})(:,3);
                plot(get_percentage(data), data, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ left/- right rotation', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
        end

        % Add legend outside the loop
        legend('FontSize', font);
    end
end

function percentage = get_toeoff_percentage(foot_strike_frame_1, toe_off_frame, foot_strike_frame_2)
    
    cycle = (foot_strike_frame_2-foot_strike_frame_1);
    percentage = (toe_off_frame-foot_strike_frame_1)/cycle * 100;
end

function plot_joint_moment(data_structure, left_toeoff_event, right_toeoff_event, weight)

    font = 14;

    % Plotting all data in the structure
    fields = fieldnames(data_structure); % Get the field names of the structure
    
    for i = 1:numel(fields)
        field = strrep(fields{i}, '_', ' ');
        switch(field)
            case('Rhip')
                figure(3)
                subplot(3,3,7)
                hold on
                data = data_structure.(fields{i}).localM(:,1);
                plot(get_percentage(data), data/weight, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ extension/- flexion (Nm/kg)', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');

                subplot(3,3,8)
                hold on
                data = data_structure.(fields{i}).localM(:,2);
                plot(get_percentage(data), data/weight, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ addcution/- abdcution (Nm/kg)', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');

                title('Hip Joint Moment in 3D', FontSize=font)
                
                subplot(3,3,9)
                hold on
                data = data_structure.(fields{i}).localM(:,3);
                plot(get_percentage(data), data/weight, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ int/- ext rotation (Nm/kg)', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                
%                 legend({ '+ flexion/- extension', '+ addcution/- abdcution', '+ internal/- external rotation'}, FontSize=font)

            case('Lhip')
                figure(3)
                subplot(3,3,7)
                hold on
                data = data_structure.(fields{i}).localM(:,1);
                plot(get_percentage(data), data/weight, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ extension/- flexion (Nm/kg)', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-1 1.5])

                subplot(3,3,8)
                hold on
                data = data_structure.(fields{i}).localM(:,2);
                plot(get_percentage(data), -data/weight, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ addcution/- abdcution (Nm/kg)', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-1 1.5])
                
                subplot(3,3,9)
                hold on
                data = data_structure.(fields{i}).localM(:,3);
                plot(get_percentage(data), -data/weight, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ int/- ext rotation (Nm/kg)', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-0.5 0.5])
%                 legend({ '+ flexion/- extension', '+ abduction/- addcution', '+ external/- internal rotation'}, FontSize=font)

            case('Rankle')
                figure(3)
                subplot(3,3,1)
                data = data_structure.(fields{i}).localM(:,1);
                plot(get_percentage(data), data/weight, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ plantar/- dorsiflexion (Nm/kg)', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                

                subplot(3,3,2)
                data = data_structure.(fields{i}).localM(:,2);
                plot(get_percentage(data), data/weight, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ inversion/- eversion (Nm/kg)', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');

                title('Ankle Joint Moment in 3D', FontSize=font)

                subplot(3,3,3)
                data = data_structure.(fields{i}).localM(:,3);
                plot(get_percentage(data), data/weight, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ int/- ext rotation (Nm/kg)', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                

            case('Lankle')
                figure(3)
                subplot(3,3,1)
                hold on
                data = data_structure.(fields{i}).localM(:,1);
                plot(get_percentage(data), data/weight, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ plantar/- dorsiflexion (Nm/kg)', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-1 1.5])

                subplot(3,3,2)
                hold on
                data = data_structure.(fields{i}).localM(:,2);
                plot(get_percentage(data), -data/weight, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ inversion/- eversion (Nm/kg)', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-0.5 0.5])
                
                subplot(3,3,3)
                hold on
                data = data_structure.(fields{i}).localM(:,3);
                plot(get_percentage(data), -data/weight, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ int/- ext rotation (Nm/kg)', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-0.5 0.5])

            case('Rknee')
                figure(3)
                subplot(3,3,4)
                data = data_structure.(fields{i}).localM(:,1);
                plot(get_percentage(data), -data/weight, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ extension/- flexion (Nm/kg)', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                

                subplot(3,3,5)
                data = data_structure.(fields{i}).localM(:,2);
                plot(get_percentage(data), data/weight, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ varus/- valgus (Nm/kg)', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                

                title('Knee Joint Moment in 3D', FontSize=font)

                subplot(3,3,6)
                data = data_structure.(fields{i}).localM(:,3);
                plot(get_percentage(data), data/weight, '-g', LineWidth=2, DisplayName='right');
                ylabel('+ int/- ext rotation (Nm/kg)', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                

%                 legend({ '+ extension/- flexion', '+ addcution/- abdcution', '+ internal/- external rotation'}, FontSize=font)
            case('Lknee')
                figure(3)
                subplot(3,3,4)
                hold on
                data = data_structure.(fields{i}).localM(:,1);
                plot(get_percentage(data), -data/weight, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ extension/- flexion (Nm/kg)', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-0.5 1])

                subplot(3,3,5)
                hold on
                data = data_structure.(fields{i}).localM(:,2);
                plot(get_percentage(data), -data/weight, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ varus/- valgus (Nm/kg)', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-0.5 1])
                
                subplot(3,3,6)
                hold on
                data = data_structure.(fields{i}).localM(:,3);
                plot(get_percentage(data), -data/weight, '-r', LineWidth=2, DisplayName='left');
                ylabel('+ int/- ext rotation (Nm/kg)', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-0.5 1])
%                 legend({ '+ extension/- flexion', '+ abduction/- addcution', '+ external/- internal rotation'}, FontSize=font)
            
        end

        % Add legend outside the loop
        legend('FontSize', font);
    end

end

function plot_joint_power(data_structure, left_toeoff_event, right_toeoff_event,weight)

    font = 14;

    % Plotting all data in the structure
    fields = fieldnames(data_structure); % Get the field names of the structure
    
    for i = 1:numel(fields)
        field = strrep(fields{i}, '_', ' ');
        switch(field)
            case('Rhip')
                figure(5)
                subplot(2,3,1); hold on
                data = data_structure.(fields{i}).sumpower;
                plot(get_percentage(data), data/weight, '-g', LineWidth=2, DisplayName='right');
                ylabel('- Absorptions/+ Generation (W/kg)', FontSize=font);
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');

            case('Lhip')
                figure(5)
                subplot(2,3,1); hold on
                data = data_structure.(fields{i}).sumpower;
                plot(get_percentage(data), data/weight, '-r', LineWidth=2, DisplayName='left');
                ylabel('- Absorptions/+ Generation (W/kg)', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                ylim([-2 2])
                title('Hip Joint power', FontSize=font)

            case('Rankle')
                figure(5)
                subplot(2,3,2); hold on
                data = data_structure.(fields{i}).sumpower;
                plot(get_percentage(data), data/weight, '-g', LineWidth=2, DisplayName='right');
                ylabel('- Absorptions/+ Generation (W/kg)', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');

            case('Lankle')
                figure(5)
                subplot(2,3,2); hold on
                hold on
                data = data_structure.(fields{i}).sumpower;
                plot(get_percentage(data), data/weight, '-r', LineWidth=2, DisplayName='left');
                ylabel('- Absorptions/+ Generation (W/kg)', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-2 5])
                title('Ankle Joint power', FontSize=font)

            case('Rknee')
                figure(5)
                subplot(2,3,3)
                data = data_structure.(fields{i}).sumpower;
                plot(get_percentage(data), data/weight, '-g', LineWidth=2, DisplayName='right');
                ylabel('- Absorptions/+ Generation (W/kg)', FontSize=font);
                hold on
                xline(right_toeoff_event, '--', 'Color', 'g', 'LineWidth', 1, DisplayName='right toeoff');yline(0, '--', 'Color', 'k', 'LineWidth', 1, DisplayName='zero');
                
            case('Lknee')
                figure(5)
                subplot(2,3,3)
                hold on
                data = data_structure.(fields{i}).sumpower;
                plot(get_percentage(data), data/weight, '-r', LineWidth=2, DisplayName='left');
                ylabel('- Absorptions/+ Generation (W/kg)', FontSize=font);
                xline(left_toeoff_event, '--', 'Color', 'r', 'LineWidth', 1, DisplayName='left toeoff');
                ylim([-2 2])
                title('Knee Joint power', FontSize=font)
        end

        % Add legend outside the loop
        legend('FontSize', font);
    end
end