classdef Sinsquare < shape.C2boundary
 
    properties(SetAccess = protected)
		scale = 1;
        phi = 0;
        degree = 4 ;
    end
    
    methods
        
      function obj = Sinsquare(a,nbPoints,dspl)
            if nargin < 3
               dspl = 1 ; 
            end
            b = a/2.3 ;
            w = 1 ;        
          
            selector1  = floor(nbPoints/4 )  ;

            edge1 = [fliplr(linspace(-1,1,selector1));1+fliplr(0.1*sin(linspace(-w*pi,w*pi,selector1)))] ;
            edge1 = edge1(:,2:end-1) ;
            
            edge2 = [fliplr(0.1*sin(linspace(-w*pi,w*pi,selector1)))-1;fliplr(linspace(-1,1,selector1))]  ;
            edge2 =  edge2(:,2:end-1) ;
            
            edge3 =  fliplr(edge1 - repmat([0;2],1,length(edge1)))  ;
            edge4 =  fliplr(edge2 + repmat([2;0],1,length(edge2)))  ;
            
            points = b*[edge1,edge2,edge3,edge4] ;
            points = circshift(points,[0,selector1 /2]); 
            points = points - points(:,1) ;
            
            theta = 0:(1/length(points)):(1-(1/length(points))) ;            
            [points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl);            

            obj = obj@shape.C2boundary(points, tvec, avec, normal, [0;0], 'Sin square');
            
            obj.scale = b;
            obj.phi = 0;
      end
      
        function obj = mtimes(obj, s)
			% Overload of the operator *
			obj = mtimes@shape.C2boundary(obj, s);
			obj.size_a = obj.size_a * s;
		
		end
		
		function obj = lt(obj, phi)
			% Redefine the < operator as the rotation of the boundary
			obj = lt@shape.C2boundary(obj,phi);
			obj.phi = obj.phi+phi;
        end      
    end
end