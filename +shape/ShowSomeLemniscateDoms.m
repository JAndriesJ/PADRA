function LemniscateDom = ShowSomeLemniscateDoms(LemNum,nbPoints)
    Domain = cell(19,1) ;
    dpsl = 1;
    
    foci = cell(1,9);
    foci{1} = [ 1.028217386556707  -0.178663690857326;
                0.403957729709384   0.382422178075320];
    foci{2} = [ 0.810360585773191  -0.038036076757269;
                0.625308090908918   0.264849291346766];
    foci{3} = [ 1.266431972030214  -0.383478043790057;
                0.407818154899291   0.351390839037181];
    foci{4} = [ 0.888989313310634   0.338395131499682  -0.438784042334693;
                0.260481131512223   0.315406797883773   0.504719453341981];
    foci{5} = [ 0.681471981358612  -0.420917152579725   0.683768368529442;
                0.761265156042021   0.677878891153285  -0.087271383248153];
    foci{6} = [ 1.372859368087138   0.503519450922174  -0.548238253654143;
                0.486065323362062   0.103745570613663   0.523887412265618];
    foci{7} = [ 0.943948570372622   0.706842234760497  -0.278244559388281;
                0.467248343536628   0.324452110484066   0.498808054589172];
    foci{8} = [-0.502468829623819   0.974995620940245   0.973104643980599;
                0.255626486680377   0.557018870743030   0.198179037385723];
    foci{9} = [ 1.472640571692003  -0.119246149585907  -0.243481878483808;
                0.048572098416752   0.653795642874892   0.764312867406655];
    foci{10}= [ 0.222164712278007   0.229512169809647   0.837986073749624;
                0.719359574604910   0.144909870411321   0.555513094859472];
    if nargin < 1
        nbPoints = 1000 ;
        for iterator1 = 1:9
            Domain{iterator1} =...
                shape.Lemniscate.Lemniscate(foci{iterator1},nbPoints,dpsl);   
            hold on
            subtightplot(3,3,iterator1)
            plot(Domain{iterator1},'Linewidth',2,'color','red')
            axis([-1 1 -1 1]*2)
            axis square
            hold off
        end
    elseif nargin < 2
        nbPoints = 1000;
        LemniscateDom = ...
            shape.Lemniscate.Lemniscate(foci{LemNum},nbPoints,dpsl);
    elseif nargin < 3
        LemniscateDom = ...
            shape.Lemniscate.Lemniscate(foci{LemNum},nbPoints,dpsl);
    end
    
end


function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
%function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
%
% Functional purpose: A wrapper function for Matlab function subplot. Adds the ability to define the gap between
% neighbouring subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the gap between
% subplots can reach 40% of figure area, which is pretty lavish.  
%
% Input arguments (defaults exist):
%   gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
%            is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%            relatively large axis. 
%   marg_h  margins in height in normalized units (0...1)
%            or [lower uppper] for different lower and upper margins 
%   marg_w  margins in width in normalized units (0...1)
%            or [left right] for different left and right margins 
%
% Output arguments: same as subplot- none, or axes handle according to function call.
%
% Issues & Comments: Note that if additional elements are used in order to be passed to subplot, gap parameter must
%       be defined. For default gap value use empty element- [].      
%
% Usage example: h=subtightplot((2,3,1:2,[0.5,0.2])

if (nargin<4) || isempty(gap),    gap=0.01;  end
if (nargin<5) || isempty(marg_h),  marg_h=0.05;  end
if (nargin<5) || isempty(marg_w),  marg_w=marg_h;  end
if isscalar(gap),   gap(2)=gap;  end
if isscalar(marg_h),  marg_h(2)=marg_h;  end
if isscalar(marg_w),  marg_w(2)=marg_w;  end
gap_vert   = gap(1);
gap_horz   = gap(2);
marg_lower = marg_h(1);
marg_upper = marg_h(2);
marg_left  = marg_w(1);
marg_right = marg_w(2);

%note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
[subplot_col,subplot_row]=ind2sub([n,m],p);  

% note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot 
subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot   

% single subplot dimensions:
%height=(1-(m+1)*gap_vert)/m;
%axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
height=(1-(marg_lower+marg_upper)-(m-1)*gap_vert)/m;
%width =(1-(n+1)*gap_horz)/n;
%axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
width =(1-(marg_left+marg_right)-(n-1)*gap_horz)/n;

% merged subplot dimensions:
merged_height=subplot_rows*( height+gap_vert )- gap_vert;
merged_width= subplot_cols*( width +gap_horz )- gap_horz;

% merged subplot position:
merged_bottom=(m-max(subplot_row))*(height+gap_vert) +marg_lower;
merged_left=(min(subplot_col)-1)*(width+gap_horz) +marg_left;
pos_vec=[merged_left merged_bottom merged_width merged_height];

% h_subplot=subplot(m,n,p,varargin{:},'Position',pos_vec);
% Above line doesn't work as subplot tends to ignore 'position' when same mnp is utilized
h=subplot('Position',pos_vec,varargin{:});

if (nargout < 1),  clear h;  end

end