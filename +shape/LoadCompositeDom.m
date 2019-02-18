% Visual of all shape in second paper
function Domain = LoadCompositeDom(DomNumb,a,angle,nbPoints,dspl)
    if nargin < 2
        a = 1 ;
        angle = 2*pi/5 ;
        nbPoints = 1000 ;
        dspl = 30 ;
    elseif nargin < 3
        angle = 2*pi/5 ;
        nbPoints = 1000 ;
        dspl = 30 ;
    elseif nargin < 4
        nbPoints = 1000 ;
        dspl = 30 ;
    elseif nargin < 5
        dspl = 30 ;
    end
    
    switch DomNumb
        case 1
            D = shape.CompositeDoms.HandFan(a,nbPoints,dspl) ;          
        case 2
            D = shape.CompositeDoms.Sector(a,angle,nbPoints,dspl) ; 
        case 3
            D = shape.CompositeDoms.Sector(a*0.7,(angle+pi),nbPoints,dspl) ;
        case 4
            D = shape.CompositeDoms.Conjoined(a,nbPoints,a*0.5,a*0.7,dspl) ;
        case 5
            D = shape.CompositeDoms.Crescent(a,nbPoints,a*0.5,a*0.7,dspl) ;
        case 6
            D = shape.CompositeDoms.Vesicapisces(a,nbPoints,a*0.5,a*0.7,dspl) ;        
        case 7
            D = shape.CompositeDoms.Shell(a,nbPoints,dspl) ;
        case 8
            D = shape.CompositeDoms.HalfPetal(1.5*a,nbPoints,dspl) ;                
        case 9
            D = shape.CompositeDoms.Capsual(a,nbPoints,dspl) ;                    
        case 10
            D = shape.CompositeDoms.Egg(0.9*a/2,nbPoints,dspl) ;       
        case 11
            D = shape.CompositeDoms.Sinsquare(a,nbPoints,dspl);     
        case 12
            D = shape.CompositeDoms.LogShield(a,nbPoints,dspl) ;      
        case 13 
            D = shape.CompositeDoms.Candy(a,nbPoints,dspl); 
        case 14 
            D = shape.CompositeDoms.TriloShell(a,nbPoints,dspl) ;                
        case 15 
            D = shape.CompositeDoms.ArrowHead(a,nbPoints,dspl*2) ;
        case 16            
            D = shape.CompositeDoms.Zircon(a*0.7,nbPoints,dspl*2) ; 
        case 17 
            D = shape.CompositeDoms.Kite(a*0.7,nbPoints,dspl*2) ;
        case 18 
            D = shape.CompositeDoms.Toenail(0.9*a/2,nbPoints,dspl*2) ;  
        case 19 
            D = shape.CompositeDoms.Toll(a,nbPoints,dspl) ; 
    end
Domain = D ;

end