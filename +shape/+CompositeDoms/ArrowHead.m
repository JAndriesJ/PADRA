classdef ArrowHead< shape.C2boundary
    % Class for regular n-polygon shape 
    properties(SetAccess = protected)
        a= 1; 
		phi = 0; % orientation 
        degree = 5 ;
    end
    
    methods
      function obj = ArrowHead(a,nbPoints,dspl)  
            if nargin < 3
               dspl = 1 ; 
            end
            selector = 0:1/floor(nbPoints/3):1 ;

            edge1 = [a*cos(pi-(pi/2)*selector)+a; a*sin(pi-(pi/2)*selector)] ;
            edge2 = [a*(-2*selector+1);a*ones(1,floor(nbPoints/3)+1)] ;
            edge3 = [a*cos((pi/2)*fliplr(selector))-a; a*sin((pi/2)*fliplr(selector))] ; 
            
            points = [edge1, edge2, edge3] ;            
            points = circshift(points,[0,floor(nbPoints/5)]); % This is so that the begining and end are centered on an edge        
            points = points - points(:,floor(nbPoints/1.4)) ;
            
            theta = 0:(1/length(points)):(1-(1/length(points))) ;            
            [points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl);

            obj = obj@shape.C2boundary(points, tvec, avec, normal, [0;0], 'Arrowhead');
            
            obj.a = 0;
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