function [u,teta] = rot2uth (R)

% Use: [u,teta] = rot2uth(R)
%
% Compute the axis u (3x1 vector)
% and the rotation angle th (in radians)

warning('rot2uth is not tested')

if isrot(R, 1e-6) == 0
	I = eye(3);
	
	if abs(norm(R-I)) < 1.0e-10 % rotazione di 0 gradi; il vettore u e` messo pari a [1 0 0]
			teta = 0;
			u(1) = 1;
			u(2) = 0;
			u(3) = 0;
	elseif abs(norm(R-R')) < 1.0e-6 % rotazione di 180 gradi
	   		teta = pi;
	   		M = 0.5*(R+I);
	   		us(1) = sqrt(M(1,1));
	   		us(2) = sqrt(M(2,2));
	   		us(3) = sqrt(M(3,3));
	   		if abs(us(1)) > 1.0e-6
	   			u(1) = us(1);
	   			u(2) = M(1,2)/us(1);
	   			u(3) = M(1,3)/us(1);
	   		elseif abs(us(2)) > 1.0e-6
	    		u(2) = us(2);
	   			u(1) = M(2,1)/us(2);
	   			u(3) = M(2,3)/us(2);
	   		elseif abs(us(3)) > 1.0e-6
	    		u(3) = us(3);
	   			u(1) = M(3,1)/us(3);
	   			u(2) = M(3,2)/us(3);
	  		end
	
	else
    	s(4) = .5 * sqrt(1 + trace(R));
		s(1) = (R(3,2) - R(2,3))/(4 * s(4));
    	s(2) = (R(1,3) - R(3,1))/(4 * s(4));
		s(3) = (R(2,1) - R(1,2))/(4 * s(4));
    	amez = acos(s(4));
		a = 2*amez;
    	teta = a;
	 	sa=sin(amez);
		u(1) = s(1)/sa;
		u(2) = s(2)/sa;
		u(3) = s(3)/sa;   
	end
    u = u';
	nu=norm(u);
	u=u/nu;

else
	disp ('Error in input matrix')
	x='ERROR';
end
