classdef LogShield < shape.C2boundary
 
    properties(SetAccess = protected)
		scale = 1 ;
        phi = 0 ;
        degree= 5 ;
    end
    
    methods
        function obj =LogShield(a,nbPoints,dspl)
            if nargin < 3
                dspl = 1 ;
            end
            a = a/4 ;
            selector1  = floor(nbPoints/3 )  ;

            edgeP1 = [linspace(1,a*exp(2),selector1);log(linspace(1,a*exp(2),selector1))] ;
            edgeP2 = fliplr([edgeP1(2,end)*cos(linspace(0,pi,selector1) - pi/2) + edgeP1(1,end); edgeP1(2,end)*sin(linspace(0,pi,selector1) - pi/2) ]) ;
            edgeP3 = fliplr([edgeP1(1,:);-edgeP1(2,:)]) ;
            
            points = [edgeP1,edgeP2,edgeP3] ;
            points = circshift(points,[0,floor(selector1/3)]);
            points = points - points(:,1) ;
        
            theta = 0:(1/length(points)):(1-(1/length(points))) ;
            [points, tvec, avec, normal] = shape.C2boundary.rescale(points, theta, nbPoints, [], dspl);            

            obj = obj@shape.C2boundary(points, tvec, avec, normal, [0;0], 'sector');

            obj.scale = a;
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