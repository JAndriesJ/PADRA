classdef Candy< shape.C2boundary
    % Class for regular n-polygon shape 
    properties(SetAccess = protected)
        a= 1 ; 
		phi = 0 ; % orientation 
        degree = 8 ;
    end
    
    methods
      function obj = Candy(a,nbPoints,dspl)  
            if nargin < 3
               dspl = 1 ; 
            end
            b = a/2 ;
            c = b/2 ;
            selector1 = 0:1/floor(nbPoints/8):1 ;
            nbSelc = length(selector1) ;
            
            edge1 = [b*cos((pi/2)*selector1); b*sin((pi/2)*selector1)] ;
            edge1 = edge1(:,2:(end-1)) ;
            edge2 = [zeros(1,nbSelc);(b-c*selector1)] ;
            
            edge3 = [c*cos((pi/2)*selector1+pi/2); c*sin((pi/2)*selector1+pi/2)] ;
            edge3 = edge3(:,2:(end-1)) ;
            edge4 = [(-c-c*selector1); zeros(1,nbSelc)] ;

            edge5 = [b*cos((pi/2)*selector1+pi); b*sin((pi/2)*selector1+pi)] ;
            edge5 = edge5(:,2:(end-1)) ;
            edge6 = [zeros(1,nbSelc);(-b+c*selector1)] ;

            edge7 = [c*cos((pi/2)*selector1+3*pi/2); c*sin((pi/2)*selector1+3*pi/2)] ;
            edge7 = edge7(:,2:(end-1)) ;
            edge8 = [(c+(b/2)*selector1); zeros(1,nbSelc)] ;

            points = [edge1, edge2, edge3, edge4, edge5, edge6, edge7, edge8] ;            
            points = circshift(points,[0,floor(nbPoints/pi)]); % This is so that the begining and end are centered on an edge        
            points = points - points(:,1) ;
            
            theta = 0:(1/length(points)):(1-(1/length(points))) ;            
            [points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl);

            obj = obj@shape.C2boundary(points, tvec, avec, normal, [0;0], 'Candy');
            
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