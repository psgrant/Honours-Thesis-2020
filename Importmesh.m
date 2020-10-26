%% Import .msh data
clear all
clc
% close all
file = ('Groundwater_coarse.msh');
mshfile = fileread(file);
 
% find and classify entities
% # points - # curves - # surfaces - # volumes
% Point tag - X loc - Y loc - Z loc -# phys tag - phys tag...note
%
% Curve tag - minX - minY - minZ - maxX - maxY - maxZ -
%            #phys tags - phystag - #bounding curves - pointTag
%


% Use regular expressions to isolate each section
expression = '\$Entities\s\s(.+)\$EndEntities';
% [startind,endind]
p = regexp(mshfile,expression,'tokens');
entstr = string(p(1));
p = regexp(mshfile,'\$Nodes(.+)\$EndNodes','tokens');
nodestr = string(p(1));
p = regexp(mshfile,'\$Elements(.+)\$EndElements','tokens');
elestr = string(p(1));
%% Number of entities

% Scan the string and create cell array with each row
entcell = textscan(entstr,'%s','Delimiter','\r');
numents = str2num(entcell{1}{1})';
Entity_Type = {'Points';'Curves';'Surfaces';'Volumes'};
% Put data in table:
E = table(Entity_Type,numents);


%% Points - no tags

%
Poffset = 1;

points = zeros(numents(1),5);
for i = 1:numents
    points(i,:) = str2num(entcell{1}{i+Poffset});
end

P = array2table(points,'variablenames',...
    {'Point_Number', 'Xpos', 'Ypos', 'Zpos', 'Num_Tags'});
Coffset = Poffset + numents(1);

%% Curves - not tags
curves = zeros(numents(2),11);
for i = 1:numents(2)
    curves(i,:) = str2num(entcell{1}{i+Coffset});
end
C = array2table(curves, 'variablenames',{'Curve_Num','Xmin','Ymin','Zmin'...
    ,'Xmax','Ymax','Zmax','Num_Tags','Tags','Num_Bounding_Points','Point_Tag'});


%% isolate nodedata from node string.

[Nodedata, cornernodes, linenodes, surfnodes,numnodes] = isolateNodeData(nodestr,numents);

%% Reallocating corner and line nodes into the surfnodes category for the 
% internal structure of the aquifer

% Corner nodes
BoundaryPoints = 9;

temp = cornernodes(BoundaryPoints:end,:);
cornernodes(BoundaryPoints:end,:) = [];
numents(1) = size(cornernodes,1);

% Line nodes
for i = BoundaryPoints:size(linenodes,1)
    
   temp = [temp; linenodes{i}]; 
    
    
end
linenodes(BoundaryPoints:end) = [];

% Interior Nodes

surfnodes = [temp, zeros(size(temp,1),1);surfnodes];




%% Elements

[elements, BoundaryElements] = getElementInfo(elestr,numents);
numelements = size(elements,1);
%% RCM ORDERING

[Nodedata,RCM,connectionMatrix] = applyRCM(numnodes,numelements,Nodedata,elements,0);



%% NODE INDEXES

% Adjust ther indicies of the corner nodes using RCM ordering

for i = 1:numents(1)
    cornernodes(i,1) = Nodedata(RCM==cornernodes(i,1),1);
end

% RCM ordeing of boundarynodes
for i = 1:size(linenodes,1)
    for j = 1:size(linenodes{i},1)
        %         keyboard
        linenodes{i}(j,1) = Nodedata(RCM==linenodes{i}(j,1),1);
    end
end

for i = 1:size(surfnodes,1)
    surfnodes(i,1) = Nodedata(RCM==surfnodes(i,1),1);
    
end



%% - ADJUST ELEMENTS LIST

etemp = elements(:,2:4);
elist = etemp(:);
etemp2 = etemp;

for i = 1:numnodes
    Index = sum((RCM == i).*(1:numnodes));
    etemp(etemp2 == i) = Index;
end

elements(:,2:4) = etemp;


% Apply RCM ordering to the 1D elements on the boundaries

for i = 1:numents(2)
    for j = 1:size(BoundaryElements{i},1)
        BoundaryElements{i}(j,2) = nonzeros((RCM == BoundaryElements{i}(j,2)).*((1:numnodes)));
        BoundaryElements{i}(j,3) = nonzeros((RCM == BoundaryElements{i}(j,3)).*((1:numnodes)));
        
        % To get the correct element numbers for the 1D boundary elements :D
        % This is done by searching through the elements array to see which
        % element has the the 2 nodes in it.
        
        ElementCheck = sum((elements == BoundaryElements{i}(j,2)) + (elements == BoundaryElements{i}(j,3)),2);
        [~,I] = max(ElementCheck);
        BoundaryElements{i}(j,1) = elements(I,1);
        
        
    end
end








%% Node Class Allocations

% Corner nodes get index 1, boundary nodes get index 2 and internal nodes
% get index 3.

% Corner nodes
Nodedata(cornernodes(:,1)) = 1;
% Boundary Nodes

for i = 1:size(linenodes,1)
    if i == 3
        Nodedata(linenodes{i}(:,1)) = 2;
    else
        Nodedata(linenodes{i}(:,1)) = 2;
    end
end
% Internal Nodes
Nodedata(surfnodes(:,1)) = 3;



%% This code find what elements are connected to each node

% Corner nodes
corlen = length(cornernodes);
corele = cell(corlen,1);

for i = 1:corlen
    j = cornernodes(i,1);
    index = nonzeros(sum(elements(:,2:4) == j,2).*(1:length(elements))');
    corele{i} = [j.*ones(size(elements(index,:),1),1), elements(index,:)];
    
end

% Boundary line time :)

linelen = size(linenodes,1);
linele = cell(linelen,1);
temple = linele;


for i = 1:linelen % For each boundary
    
    for j = 1:size(linenodes{i},1) % For each node on boundary
        k = linenodes{i}(j,1);
        index = nonzeros(sum(elements(:,2:4) == k,2).*(1:length(elements))');
        temple = [k.*ones(size(elements(index,:),1),1), elements(index,:)];
        linele{i} = [linele{i}; temple];
    end
end

% Internal node boyos


surflen = length(surfnodes);
tempse = [];
surfele = zeros(100,5);
curnode = surfnodes(1,1);
pointer1 = 1;
pointer2 = 0;
for i = 1:length(surfnodes(:,1))
    curnode = surfnodes(i,1);
    
    index = nonzeros(sum(elements(:,2:4) == curnode,2).*(1:length(elements))');
    tempse = [(surfnodes(curnode)).*ones(size(elements(index,:),1),1), elements(index,:)];
    pointer2 = pointer1 + length(index);
    
    if pointer2 > length(surfele)
        surfele(length(surfele)+100,1) = 0;
    end
    
    surfele(pointer1:pointer2-1,:) = [curnode.*ones(size(elements(index,:),1),1), elements(index,:)];
    pointer1 = pointer2;
    curnode = curnode + 1;
    
end

surfele = surfele(1:(pointer2-1),:);


% [NODE NUMBER, ELEMENT NUMBER, NODE 1, NODE 2, NODE 3]
eledata.corners = corele;
eledata.lines = linele;
eledata.internal = surfele;

%% Find out how many elements are in a CV

instances = sparse(zeros(length(elements),numnodes));
for i = 1:numnodes
    
    noi = i;
    instance = elements(:,2:end) == noi;
    instances(:,i) = sum(instance,2).*((1:length(elements))');
    
    
end

elecentx = zeros(length(elements),1);
elecenty = elecentx;
nodexmid = zeros(length(elements),3);
nodeymid = nodexmid;

for j = 1:length(elements)
    
    % Centroid of all the elements
    nodeeles = elements(j,2:4);
    nodex = Nodedata(nodeeles,2);
    nodey = Nodedata(nodeeles,3);
    elecentx(j) = mean(nodex);
    elecenty(j) = mean(nodey);
    
    % middle point along element borders
    nodexmid(j,1) = 0.5 * (nodex(1) + nodex(2));
    nodexmid(j,2) = 0.5 * (nodex(3) + nodex(2));
    nodexmid(j,3) = 0.5 * (nodex(1) + nodex(3));
    
    nodeymid(j,1) = 0.5 * (nodey(1) + nodey(2));
    nodeymid(j,2) = 0.5 * (nodey(3) + nodey(2));
    nodeymid(j,3) = 0.5 * (nodey(1) + nodey(3));
end

nodexmid = [nodexmid(:,1) nodexmid(:,2) nodexmid(:,3)];
nodeymid = [nodeymid(:,1) nodeymid(:,2) nodeymid(:,3)];

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% NORMAL VECTORS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scaleFactor = 1;


elementNormalsx = zeros(numelements,2,3);
elementNormalsy = zeros(numelements,2,3);
elementNormalsx2 = zeros(numelements,2,3);
elementNormalsy2 = zeros(numelements,2,3);
SCVFaceLengths = zeros(numelements,3);
SCVFaceLengths2 = SCVFaceLengths;
SCVF1 = [0 0 0];
SCVF2 = [0 0 0];
SCVNX1 = [0 0 0];
SCVNX2 = [0 0 0];
SCVNY1 = [0 0 0];
SCVNY2 = [0 0 0];
SCVNx = zeros(numelements,3);
SCVNy = zeros(numelements,3);
SCVNfinal2x = zeros(numelements,3);
SCVNfinal2y = zeros(numelements,3);
% for each element
% tic
for e = 1:size(elements,1)
    % Extract the nodes on element e
    elementNodes = elements(e,2:end);
    
    SCVF1 = [0 0 0];
    SCVNX1 = [0 0 0];
    SCVN1 = zeros(3,2);
    SCVN12 = zeros(3,2);
    SCVNX2 = [0 0 0];
    SCVNY1 = [0 0 0];
    SCVNY2 = [0 0 0];
    % For each Sub-Control Volume in the element
    for SCV = 1:3
        
        % Declare borer and midpoint x, y data
        elementBorder = [nodexmid(e,SCV), nodeymid(e,SCV)];
        elementMidpoint = [elecentx(e), elecenty(e)];
        vectorLength = sqrt(sum((elementBorder - elementMidpoint).^2));
        SCVF1(SCV) = vectorLength;
        
        % The vector of ther SCV face
        FaceVec = elementBorder - elementMidpoint;
        
        if FaceVec(2) == 0
            FaceNorm1 = [0 -1];
            FaceNorm2 = [0 1];
            
        else
            
            % This actually calculates the normal vectors
            FaceNormy = FaceVec(1) / FaceVec(2);
            FaceNorm1 = [1 -FaceNormy];
            FaceNorm1 = FaceNorm1 ./ (norm(FaceNorm1));
            FaceNorm2 = [-1 FaceNormy];
            FaceNorm2 = FaceNorm2 ./ (norm(FaceNorm2));
            
        end
        FaceNorm = [FaceNorm1; FaceNorm2];
        AdditonConst = elementMidpoint - Nodedata(elementNodes(SCV),2:3);
        
        [~,I] = max([norm(FaceNorm1 + AdditonConst),...
            norm(FaceNorm2 + AdditonConst)]);
        SCVN1(SCV,:) = FaceNorm(I,:);
        SCVN2(SCV,:) = FaceNorm(3-I,:);
    end
    
    % SCVFINAL is the array of all outwar pointing normal vectors
    SCVNx(e,:) = SCVN1(:,1)';
    SCVNy(e,:) = SCVN1(:,2)';
    SCVN2 = SCVN2([3 1 2],:);
    SCVNfinal2x(e,:) = SCVN2(:,1)';
    SCVNfinal2y(e,:) = SCVN2(:,2)';
    SCVF2 = SCVF1(:,[3 1 2]);
    SCVNX2 = [SCVNX2(2:3), SCVNX2(1)];
    SCVNY2 = [SCVNY2(2:3), SCVNY2(1)];
    SCVFaceLengths(e,:) = SCVF1;
    SCVFaceLengths2(e,:) = SCVF2;
    
    
end

% Scale Factor
sf = 0.5;

% Scaled normal vectors for plotting
SCVNormalsx1 = SCVNx * sf .* SCVFaceLengths;
SCVNormalsy1 = SCVNy * sf .* SCVFaceLengths;
SCVNormalsx2 = SCVNfinal2x * sf .* SCVFaceLengths2;
SCVNormalsy2 = SCVNfinal2y * sf .* SCVFaceLengths2;

% Translate vectors around the middle of the element for BLUE vectors
SCVNormalsx2 = SCVNormalsx2(:,[2 3 1]);
SCVNormalsy2 = SCVNormalsy2(:,[2 3 1]);


% SCV normal vectors and facelengths multiplied together
SCVNLenx = SCVNx .* SCVFaceLengths;
SCVNLeny = SCVNy .* SCVFaceLengths;



% Calculate the face lengths of the 1D boundary elements
% This is given as half the distance between the nodes
BoundaryLengths = cell(size(BoundaryElements));

for i = 1:numents(2)
    Nodes = BoundaryElements{i}(:,2:3);
    xdist = Nodedata(Nodes(:,1),2) - Nodedata(Nodes(:,2),2);
    ydist = Nodedata(Nodes(:,1),3) - Nodedata(Nodes(:,2),3);
    BoundaryLengths{i} = 0.5 * sqrt(xdist.^2+ydist.^2);
end

%%
CVAreas = zeros(numnodes,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  CV AREAS  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is for the corner nodes my man
for j = 1:length(cornernodes)
    
    % Node of interest (noi) location - CORNER NODE
    noilocs = cornernodes(j,2:3);
    
    for i = 1:size(corele{j},1) % for each element that is connected via the noi
        
        % Nodes of a selected element
        noi = corele{j}(i,:);
        SCVAreasnode = sum((noi(1) == noi(3:5)).*(1:3));
        % Finding the distance of all the midpointsd along the element
        % border.
        diffx = abs(noilocs(1) - nodexmid(noi(2),:));
        diffy = abs(noilocs(2) - nodeymid(noi(2),:));
        dist = (diffx.^2 + diffy.^2).^0.5;
        
        
        [~,I] = sort(dist); % Sort the distances in acending order
        
        % Set up the infill, selected the 2 closest element midpoints (from noi)
        % as they will be the SCV points needed.
        patchx = [noilocs(1), nodexmid(noi(2),I(1)),...
            elecentx(noi(2)),  nodexmid(noi(2),I(2))];
        
        patchy = [noilocs(2), nodeymid(noi(2),I(1)),...
            elecenty(noi(2)), nodeymid(noi(2),I(2))];
        
        
        px =  elecentx(noi(2)) - noilocs(1);
        
        py =  elecenty(noi(2)) - noilocs(2);
        
        
        qx = patchx(4) - patchx(2);
        
        qy =  patchy(4) - patchy(2);
        AC = 0.5 * abs((qx*py - qy*px));
        CVAreas((corele{j}(1))) = CVAreas((corele{j}(1))) + AC;
%         SCVAreas(noi(2),SCVAreasnode) = AC;
    end
end


for j = 1:size(linenodes,1)
    
    % Boundary Of Interest
    BOI = linenodes{j};
    % Gets the NOI Location
    
    % For each node on the boundary
    %     for curnode = 1:(size(BOI,1))
    % Node of interest (noi) location - BOUNDARY NODES
    
    
    for i = 1:size(linele{j},1) % for each element that is connected via the noi
        
        % Element Of Interest in BOI
        EOI = linele{j}(i,:);
        noilocs = Nodedata(EOI(1),2:3);
        % retrieve Node of interest
        %  NOI = sum((BOI(:,1) == EOI(1)).*(1:size(BOI,1))');
        
        % Nodes in the slected EOI
        elenodes = EOI(3:5);
        SCVAreasnode = sum((EOI(1) == elenodes).*(1:3));
        % Removes the NOI out of the node list for the EOI
        elenodes = elenodes(elenodes ~=  EOI(1));
        
        % Calculates the midpoints between the NOI and the other two
        % nodes in the EOI
        kitepoints(1,:) = 0.5 * (noilocs + Nodedata(elenodes(1),2:3));
        kitepoints(2,:) = 0.5 * (noilocs + Nodedata(elenodes(2),2:3));
        
        dist = (kitepoints(:,1).^2+kitepoints(:,2).^2).^0.5;
        
        [~,I] = sort(dist); % Sort the distances in acending order
        
        
        % Set up the infill, selected the 2 closest element midpoints (from EOI)
        % as they will be the SCV points needed.
        patchx = [noilocs(1), kitepoints(I(1),1), elecentx(EOI(2)),  kitepoints(I(2),1)];
        patchy = [noilocs(2), kitepoints(I(1),2), elecenty(EOI(2)), kitepoints(I(2),2)];
        
        
        % Apply infill
        px =  elecentx(EOI(2)) - noilocs(1);
        
        py =  elecenty(EOI(2)) - noilocs(2);
        
        
        qx = patchx(4) - patchx(2);
        
        qy =  patchy(4) - patchy(2);
        AC = 0.5 * abs((qx*py - qy*px));
        CVAreas(EOI(1)) = CVAreas(EOI(1)) + AC;
%         SCVAreas(EOI(2),SCVAreasnode) = AC;
    end
    %     end|
    
    
end

% Internal Nodes


for j = 1:size(surfnodes,1)
    
    
    % INTERNAL NODES
    NOI = surfnodes(j,:);
    noilocs = NOI(2:3);
    % Indexes of elements on node j
    elindex = nonzeros(sum(surfele(:,1) == NOI(1),2).*(1:size(surfele,1))');
    
    for i = 1:size(elindex,1) % for each element that is connected via the noi
        
        % Element Of Interest in BOI
        EOI = surfele(elindex(i),:);
        
        % Nodes in the slected EOI
        elenodes = EOI(3:5);
        SCVAreasnode = sum((EOI(1) == elenodes).*(1:3));
        % Removes the NOI out of the node list for the EOI
        elenodes = elenodes(elenodes ~=  EOI(1));
        
        % Gets the NOI Location
        
        
        % Calculates the midpoints between the NOI and the other two
        % nodes in the EOI
        kitepoints(1,:) = 0.5 * (noilocs + Nodedata(elenodes(1),2:3));
        kitepoints(2,:) = 0.5 * (noilocs + Nodedata(elenodes(2),2:3));
        
        % Distances of the midpoints to the node of interetst
        dist = (kitepoints(:,1).^2+kitepoints(:,2).^2).^0.5;
        
        [~,I] = sort(dist); % Sort the distances in acending order
        
        
        % Set up the infill, selected the 2 closest element midpoints (from EOI)
        % as they will be the SCV points needed.
        patchx = [noilocs(1), kitepoints(I(1),1), elecentx(EOI(2)),  kitepoints(I(2),1)];
        patchy = [noilocs(2), kitepoints(I(1),2), elecenty(EOI(2)), kitepoints(I(2),2)];
        
        % Apply infill
        
        
        px =  elecentx(EOI(2)) - noilocs(1);
        
        py =  elecenty(EOI(2)) - noilocs(2);
        
        
        qx = patchx(4) - patchx(2);
        
        qy =  patchy(4) - patchy(2);
        SCVAreas(EOI(1),SCVAreasnode) = AC;
        AC = 0.5 * abs((qx*py - qy*px));
        CVAreas(EOI(1)) = CVAreas(EOI(1)) + AC;
%         SCVAreas(EOI(2),SCVAreasnode) = AC;
    end
end

%% SCV areas

SCVAreas = zeros(size(elements,1),3);
nodeOrder = [3 1 2];
for e = 1:numelements
    
    for n = 1:3
        
       qx = (elecentx(e) - Nodedata(elements(e,n+1),2));
       qy = (elecenty(e) - Nodedata(elements(e,n+1),3));
       
       px = (nodexmid(e,n) - nodexmid(e,nodeOrder(n)));
       py = (nodeymid(e,n) - nodeymid(e,nodeOrder(n)));
       
       
       SCVAreas(e,n) = 0.5 * abs((px*qy - py*qx));
        
    end
    
    
    
end

sum(sum(SCVAreas))


%% Get node numbers and the elements connecting them

% Node number, node index, element number.

% Node index 1 - corner node: 2 - Boundary node: 3 - Internal Node
NodeElements = zeros(2*size(eledata.internal,1),3);
counter = 1;
ExitCond = 0;



% This section is for the corner nodes
for c = 1:size(eledata.corners,1)
    for i = 1:size(eledata.corners{c},1)
        NodeElements(counter,:) = [eledata.corners{c}(i,[1 2]), 1];
        counter = counter + 1;
    end
end

% This section is for the boundary nodes
for c = 1:size(eledata.lines,1)
    for i = 1:size(eledata.lines{c})
        NodeElements(counter,:) =  [eledata.lines{c}(i,[1 2]), 2];
        counter = counter + 1;
    end
end

% This section is for the internal nodes
CounterConst = size(eledata.internal,1);
NodeElements(counter:counter + CounterConst - 1,:) = [eledata.internal(:,[1 2]),3*ones(size(eledata.internal,1),1)];
InternalNodeConnections = sortrows([eledata.internal(:,[1 2]),3*ones(size(eledata.internal,1),1)]);
NodeElements(counter + CounterConst:end,:) = [];
%%
BoundaryElements(9:end) = [];
BoundaryLengths(9:end) = [];
clc
plotdata.nodepoints = cornernodes;
plotdata.nodelines = linenodes;
plotdata.surfnodes = surfnodes;



% clf

ploto = 0;

if ploto == 1
    plotmesh(C.Xmin,C.Xmax,C.Ymin,C.Ymax,C.Point_Tag,C.Num_Bounding_Points,...
        P,Nodedata,elements,nodexmid,nodeymid,elecentx,elecenty,plotdata,eledata...
        ,0,0,0,SCVNormalsx1,SCVNormalsx2,SCVNormalsy1,SCVNormalsy2);
%     pbaspect([1 0.2 1])
end



save('mesh.mat','Nodedata','elements','BoundaryElements', 'BoundaryLengths', 'cornernodes', 'CVAreas', 'linenodes','SCVFaceLengths', 'SCVNx','SCVNy','RCM','connectionMatrix','elecentx','elecenty','SCVAreas')





