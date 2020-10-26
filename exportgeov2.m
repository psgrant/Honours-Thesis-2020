%% To export the .geo file

% define boundary bpoins
pointlocs = [0 0;
    500 0;
    500 50;
    500 100;
    300 100;
    0 100;
    0 75;
    0 40];

% Define embedded points
embpointlocs = [50,  40;
                350, 40;
                350, 50;
                300, 50;
                50,  50;
                50,  55;
                50,  75];

% Define embedded lines
emblines = [8 9;
            9 10;
            10 11;
            11 12;
            12 13;
            13 9;
            11 3;
            12 5;
            14 15];

        
% Number of nodes along each embedd line
emblinenodes = [12 100 3 15 70 3 50 16 8];
% NUmber of nodes along boudary
bnn = [120 15 15 60 90 10 10 12];

% Number of embedded points
nop = size(pointlocs,1);
noep = size(embpointlocs,1);
% Setting up GMSH formatting
pointlocs = [(1:nop)',pointlocs,ones(nop,1).*0.1];
embpointlocs = [((1:noep)+nop)',embpointlocs,ones(noep,1).*0.1];

% Define booundary lines 
linenum = [(1:nop)',[(2:nop),1]'];
curvenum = 1:nop;

% Print geo string%

geoStrPoints = cell(nop,1);
%
for pointi = 1:nop
    geoStrPoints{pointi} = sprintf('Point(%d) = {%f,%f,0,%f};\r\n',...
        pointlocs(pointi,1),pointlocs(pointi,2),pointlocs(pointi,3),pointlocs(pointi,4));
end
geostrP = strjoin(geoStrPoints);

geoStrLines = cell(nop-1,1);
% export Lines
for  linei = 1:nop
    geoStrLines{linei} = sprintf('Line(%d) = {%d,%d};\r\nTransfinite Line {%d} = %d Using Progression 1;\r\n',...
        linei,linenum(linei,1),linenum(linei,2),linei,bnn(linei));
end
geostrL = strjoin(geoStrLines);

curvecell = cell(3,1);
curvecell{1} = 'Curve Loop(1) = {';
curvecell{2} = sprintf('%d,',curvenum);
curvecell{2}(end) = [];
curvecell{3} = sprintf('};\r\n Plane Surface(%d) = {%d};\r\n',1,1);
geostrC = strjoin(curvecell);


% Embeded curves
% 
geoStrPoints = cell(noep,1);
% %
for pointi = 1:noep
    geoStrPoints{pointi} = sprintf('Point(%d) = {%f,%f,0,%f};\r\n',...
        embpointlocs(pointi,1),embpointlocs(pointi,2),embpointlocs(pointi,3),embpointlocs(pointi,4));
end
embgeostrP = strjoin(geoStrPoints);


geoStrLinesEMB = cell(noep-1,1);
lineoffset = size(linenum,1);
emblineids = zeros(size(emblinenodes));
% Embeded Lines
for  linei = 1:length(emblinenodes)
    geoStrLinesEMB{linei} = sprintf('Line(%d) = {%d,%d};\r\nTransfinite Line {%d} = %d Using Progression 1;\r\n',...
        linei+lineoffset,emblines(linei,1),emblines(linei,2),linei+lineoffset,emblinenodes(linei));
    emblineids(linei) = linei+lineoffset;
end
geoStrLinesEMB = strjoin(geoStrLinesEMB);

curvecell = cell(2,1);
curvecell{1} = 'Curve{';
curvecell{2} = sprintf('%d,',emblineids);
curvecell{2}(end) = [];
curvecell{3} = sprintf('} In Surface{1};',1,1);
geoEmbLinesC = strjoin(curvecell);



% Export to .geo
totGeoStr = cat(2,geostrP,embgeostrP,geostrL,geostrC,geoStrLinesEMB,geoEmbLinesC);
fileID = fopen('groundwatertest.geo','w');
fprintf(fileID, totGeoStr);
fclose(fileID);
