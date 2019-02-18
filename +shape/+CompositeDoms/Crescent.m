classdef Crescent< shape.C2boundary
    % Class for regular n-polygon shape 
    properties(SetAccess = protected)
        radius_r1 = 1;  
		radius_r2 = 1; 
        seperation_a = 0.5 ;
		phi = 0; % orientation
        degree = 4;
    end
    
    methods
      function obj = Crescent(scales,nbPoints,r1,r2,dspl)
            a = 0.5 ;
            if nargin < 3
               r1 = 0.7;
               r2 = 0.5; 
               dspl = 1;
            elseif nargin < 5
               dspl = 1; 
            end
            
            if abs(a) > r1+r2
                error('The disks must touch')
            elseif r1<= 0 || r2 <= 0
                errot('Must use positive radii')
            end
            
            
           
            phase1 = angle((r1^2-r2^2+a^2)/(2*a) + 1i*sqrt(r1^2-((r1^2-r2^2+a^2)/(2*a))^2)) ;
            phase2 = angle((r1^2-r2^2+a^2)/(2*a) - 1i*sqrt(r1^2-((r1^2-r2^2+a^2)/(2*a))^2)) ;
            
            phase3 = angle((r1^2-r2^2+a^2)/(2*a)-a + 1i*sqrt(r1^2-((r1^2-r2^2+a^2)/(2*a))^2)) ;
            phase4 = angle((r1^2-r2^2+a^2)/(2*a)-a - 1i*sqrt(r1^2-((r1^2-r2^2+a^2)/(2*a))^2)) ;
            
            fr1 =  (r1*(2*pi-2*phase1))/(r1*(2*pi-2*phase1) + r2*(2*phase3) ); 
            fr2 =  (r2*(2*phase3))/(r1*(2*pi-2*phase1) + r2*(2*phase3) ) ;
            
            selector1 = 0:1/(floor(nbPoints*fr1)):1-1/(floor(nbPoints*fr1)) ;
            selector2 = 0:1/(floor(nbPoints*fr2)):1 ;

            points = scales*[[r1*cos((2*pi+2*phase2)*selector1+phase1);r1*sin((2*pi+2*phase2)*selector1+phase1)],...
                [r2*cos((2*pi+2*phase4)*fliplr(selector2)+phase3)+a;r2*sin((2*pi+2*phase4)*fliplr(selector2)+phase3)]] ;
            points = circshift(points,[0,3*(floor(nbPoints/(14/11)))]); % this is so that the begining and end are centered on an edge and not a corner
            points = points - points(:,1) ;
            
            theta = 0:(1/length(points)):(1-(1/length(points))) ;
            [points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl); 
            
 
            obj = obj@shape.C2boundary(points, tvec, avec, normal, [0;0], 'Crescent');
            
            obj.radius_r1 = r1 ;  
            obj.radius_r2 = r2 ; 
            obj.seperation_a = a ;
            obj.phi = 0; 
                       
        end
        function obj = mtimes(obj, s)
			% Overload of the operator *
			obj = mtimes@shape.C2boundary(obj, s);
			obj.radius_r1 = obj.radius_r1*s ;  
            obj.radius_r2 = obj.radius_r2*s ; 
            obj.seperation_a = obj.seperation_a*s ;
		end
		
		function obj = lt(obj, phi)
			% Redefine the < operator as the rotation of the boundary
			obj = lt@shape.C2boundary(obj,phi);
			obj.phi = obj.phi+phi;
		end
            
    end
end