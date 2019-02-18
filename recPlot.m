function [] = recPlot(poly, BifurcationPoints,SegmentationPoints,DomainCandidateCell,Domain)
    if nargin < 2
        BifurcationPoints  = [] ;
        SegmentationPoints = [] ;
        DomainCandidateCell = {} ;
        Domain = [] ;
    elseif nargin < 3
        SegmentationPoints = [] ;
        DomainCandidateCell = {} ;
        Domain = [] ;
    elseif  nargin < 4
        DomainCandidateCell = {} ;
        Domain = [] ;
    elseif  nargin < 5
        Domain = [] ;
    end
    if ~isempty(DomainCandidateCell)
        nbSubPlots = length(DomainCandidateCell) ;
    else
        nbSubPlots = 1 ;
    end
    
    
    % first plot of segmentation points and such
    subtightplot(ceil(sqrt(nbSubPlots)),ceil(sqrt(nbSubPlots)),1)
    hold on 
    polyplot(poly,2)
    if ~isempty(BifurcationPoints)
        scatter(BifurcationPoints(1,:), BifurcationPoints(2,:), 50, 'k', 'filled')
    end
    if ~isempty(SegmentationPoints)
        scatter(SegmentationPoints(1,:), SegmentationPoints(2,:),50,'g','<','filled')
        text(SegmentationPoints(1,:), SegmentationPoints(2,:),cellstr(num2str(SegmentationPoints(3,:)')))
    end
    scatter(0,0,70,'k','p','filled')
    axis([-1 1 -1 1]*2)
    axis square
    hold off
    if ~isempty(DomainCandidateCell)
        for iter = 2:(nbSubPlots+1)
            subtightplot(ceil(sqrt(nbSubPlots)),ceil(sqrt(nbSubPlots)),iter)
            hold on
            %polyplot(poly,2)
            %plot(Domain,'Linewidth',2,'color','magenta')
            scatter(Domain.points(1,:),Domain.points(2,:),7,'magenta','filled') 
            scatter(DomainCandidateCell{iter-1}.points(1,:),DomainCandidateCell{iter-1}.points(2,:),7,linspace(1,25,1000),'filled') 
            scatter(0,0,70,'k','p','filled')
            axis([-1 1 -1 1]*2)
            axis square
            hold off
        end
    end      
end

function [] = polyplot(poly,rannge)
    if nargin < 2
        rannge = 2 ;
    end
    f = matlabFunction(poly) ;
    fc = fcontour( f,'r') ;
    fc.LineWidth = 2 ;
    fc.LevelList = 0 ;
    fc.MeshDensity = 500 ;
    fc.XRange = [-1,1]*rannge ;
    fc.YRange = [-1,1]*rannge ;
end
