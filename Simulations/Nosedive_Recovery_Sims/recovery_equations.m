function states = recovery_equations(delta_time, max_iterations, mass, drag_coefficient, lift_coefficient, drag_offset, lift_offset, lift_direction, inertia_tensor, stable_bounds) 
	%RECOVERY EQUATIONS: 

	arguments
	delta_time (1, 1) {mustBeNumeric}
	max_iterations (1, 1) {mustBeNumeric} = 10000
	mass (1, 1) {mustBeNumeric} = 1
	drag_coefficient (1, 1) {mustBeNumeric} = 0.2
	lift_coefficient (1, 1) {mustBeNumeric} = 0.5;
	drag_offset (3, 1) {mustBeNumeric} = [0 0 0]'
	lift_offset (3, 1) {mustBeNumeric} = [0 0 -0.24]'
	lift_direction (3, 1) {mustBeNumeric} = [1 0 0]'
	inertia_tensor (:, :) {mustBeNumeric} = []
	stable_bounds (1, 1) {mustBeNumeric} = 0.005;
	end

	if isempty(inertia_tensor)
		radius = 0.25;
		height = 2;
		fprintf("Approximating the Inertia Tensor using a cylinder of radius: %fm and height: %fm\n", radius, height);
		inertia_tensor = [0.0833*mass*height*height+0.25*mass*radius*radius 0 0;
			0 0.0833*mass*height*height+0.25*mass*radius*radius 0;
			0 0 0.5*mass*radius];
	end
    	
	stable = [sqrt(2)/2 0 -sqrt(2)/2 0];

	inv_inertia_tensor = inv(inertia_tensor);

	states = struct();
	states.pos = zeros(max_iterations, 3);
	states.vel = zeros(max_iterations, 3);
	states.avel = zeros(max_iterations, 3);
	states.orient = zeros(max_iterations, 4);
    states.pitch = zeros(max_iterations, 1);
    states.pitch(1) = -pi/2;
	states.orient(1) = 1;

	inv_mass = 1 / mass;

	weight = [0 0 -1]' * mass * 9.81;

	for ii = 2:max_iterations
		R = generate_rotation_matrix(states.orient(ii-1, :));
		v = states.vel(ii-1,:)';
		I = R * inertia_tensor * R';
		I_inv = R * inv_inertia_tensor * R';
		
        if norm(v) ~= 0
		    drag = -v / norm(v) * dot(v, v) * drag_coefficient;
        else
            drag = [0 0 0]';
        end
		lift = R * lift_direction * dot(v, v) * lift_coefficient;
		
		force = drag + lift + weight;
		torque = cross(R*drag_offset, drag) + cross(R*lift_offset, lift);
		am = I * states.avel(ii-1, :)';

		delta_vel = force' * inv_mass * delta_time;
		delta_avel = I_inv * (cross(am, states.avel(ii, :)') + torque) * delta_time;
		delta_pos = states.vel(ii-1, :) * delta_time;
		delta_rot = 0.5 * quatxquat([0 states.avel(ii-1, :)], states.orient(ii-1, :)) * delta_time;

		states.pos(ii, :) = states.pos(ii-1, :) + delta_pos;
		states.vel(ii, :) = states.vel(ii-1, :) + delta_vel;
		states.orient(ii, :) = states.orient(ii-1, :) + delta_rot;
		states.orient(ii, :) = states.orient(ii, :) / norm(states.orient(ii, :));
		states.avel(ii, :) = states.avel(ii-1, :) + delta_avel';

        p = R * [0 0 -1]';
        p = norm(cross(p, [0 0 -1]'));

        states.pitch(ii)=p;

		if stable_bounds >= norm(stable - states.orient(ii, :)) && false
			states.pos(ii+1:end, :) = [];
			states.vel(ii+1:end, :) = [];
			states.orient(ii+1:end, :) = [];
			states.avel(ii+1:end, :) = [];
            states.pitch(ii+1:end, :) = [];
			break;
		end
    end

    d = min(states.pos(:, 3));
    fprintf("The payload will stablelize after falling a distance of %fm\n", d);
end

function out = quatxquat(q1, q2) 
	q10 = q1(1);
	q11 = q1(2);
	q12 = q1(3);
	q13 = q1(4);

	q20 = q2(1);
	q21 = q2(2);
	q22 = q2(3);
	q23 = q2(4);

	o0 = q10*q20-q11*q21-q12*q22-q13*q23;
	o1 = q10*q21+q20*q11+q12*q23-q13*q22;
	o2 = q10*q22+q20*q12+q13*q21-q11*q23;
	o3 = q10*q23+q20*q13+q11*q22-q12*q21;

	out = [o0 o1 o2 o3];
end

function out = generate_rotation_matrix(q)
	q0 = q(1);
	q1 = q(2);
	q2 = q(3);
	q3 = q(4);

	% First row of the rotation matrix
	r00 = 2 * (q0 * q0 + q1 * q1) - 1;
	r01 = 2 * (q1 * q2 - q0 * q3);
	r02 = 2 * (q1 * q3 + q0 * q2);

	% Second row of the rotation matrix
	r10 = 2 * (q1 * q2 + q0 * q3);
	r11 = 2 * (q0 * q0 + q2 * q2) - 1;
	r12 = 2 * (q2 * q3 - q0 * q1);

	% Third row of the rotation matrix
	r20 = 2 * (q1 * q3 - q0 * q2);
	r21 = 2 * (q2 * q3 + q0 * q1);
	r22 = 2 * (q0 * q0 + q3 * q3) - 1;

	out = [r00 r01 r02;
		r10 r11 r12;
		r20 r21 r22];
end
