function plotmesh(xmin,xmax,ymin,ymax,pointtag1,pointtag2,P,Nodedata,...
    Elements,nodexmid,nodeymid,elecentx,elecenty,plotdata,eledata...
    ,nodeLabels,ElementLabels,legendmyguy,SCVNormalsx1,SCVNormalsx2,SCVNormalsy1,SCVNormalsy2)


hold on
nodex = Nodedata(:,2);
nodey = Nodedata(:,3);
pointno = P.Point_Number;
pointtag1 = abs(pointtag1);
pointtag2 = abs(pointtag2);
pointloc = [P.Xpos P.Ypos];


areac = 0;
areab = 0;
areai = 0;
p = 0;
q = 0;
px = 0;
py = 0;
hehe = 1;
nodexmid = [nodexmid(:,1) elecentx nodexmid(:,2) elecentx nodexmid(:,3)];
nodeymid = [nodeymid(:,1) elecenty nodeymid(:,2) elecenty nodeymid(:,3)];
elemidsx = nodexmid(:,[1,3,5]);
elemidsy = nodeymid(:,[1,3,5]);



colvar = [0         0    1.0000
    1.0000         0         0
    0    1.0000         0
    0         0    0.1724
    1.0000    0.1034    0.7241
    1.0000    0.8276         0
    0    0.3448         0
    0.5172    0.5172    1.0000
    0.6207    0.3103    0.2759
    0    1.0000    0.7586
    0    0.5172    0.5862
    0         0    0.4828
    0.5862    0.8276    0.3103
    0.9655    0.6207    0.8621
    0.8276    0.0690    1.0000
    0.4828    0.1034    0.4138
    0.9655    0.0690    0.3793
    1.0000    0.7586    0.5172
    0.1379    0.1379    0.0345
    0.5517    0.6552    0.4828
    0.9655    0.5172    0.0345
    0.5172    0.4483         0
    0.4483    0.9655    1.0000
    0.6207    0.7586    1.0000
    0.4483    0.3793    0.4828
    0.6207         0         0
    0    0.3103    1.0000
    0    0.2759    0.5862
    0.8276    1.0000         0
    0.7241    0.3103    0.8276
    0.2414         0    0.1034
    0.9310    1.0000    0.6897
    1.0000    0.4828    0.3793
    0.2759    1.0000    0.4828
    0.0690    0.6552    0.3793
    0.8276    0.6552    0.6552
    0.8276    0.3103    0.5172
    0.4138         0    0.7586
    0.1724    0.3793    0.2759
    0    0.5862    0.9655
    0.0345    0.2414    0.3103
    0.6552    0.3448    0.0345
    0.4483    0.3793    0.2414
    0.0345    0.5862         0
    0.6207    0.4138    0.7241
    1.0000    1.0000    0.4483
    0.6552    0.9655    0.7931
    0.5862    0.6897    0.7241
    0.6897    0.6897    0.0345
    0.1724         0    0.3103];

ep = eledata.corners;
el = eledata.lines;
se = eledata.internal;

np = plotdata.nodepoints;
nl = plotdata.nodelines;
ns = plotdata.surfnodes;
pcnt = 1;








% Elemnts/volumes

% % Corner Nodes
% 
for j = 1:length(np)
    colsel = randi([1,length(colvar)]);
    h(pcnt) = plot(np(j,2),np(j,3),'o','color',colvar(colsel,:),'markersize',8,'linewidth',1.7);
    
    if nodeLabels == 1
        text(np(j,2),np(j,3),num2str(np(j,1)),'Color','red','FontSize',14,'FontWeight','Bold')
    end
    
    pcnt = pcnt+1;
    
    % Node of interest (noi) location - CORNER NODE
    noilocs = np(j,2:3);
    
    for i = 1:size(ep{j},1) % for each element that is connected via the noi
        
        % Nodes of a selected element
        noi = ep{j}(i,:);
        
        % Finding the distance of all the midpointsd along the element
        % border.
        diffx = abs(noilocs(1) - elemidsx(noi(2),:));
        diffy = abs(noilocs(2) - elemidsy(noi(2),:));
        dist = (diffx.^2 + diffy.^2).^0.5;
        
        
        [~,I] = sort(dist); % Sort the distances in acending order
        
        % Set up the infill, selected the 2 closest element midpoints (from noi)
        % as they will be the SCV points needed.
        patchx = [noilocs(1), elemidsx(noi(2),I(1)), elecentx(noi(2)),  elemidsx(noi(2),I(2))];
        patchy = [noilocs(2), elemidsy(noi(2),I(1)), elecenty(noi(2)), elemidsy(noi(2),I(2))];
        % Apply infill
        fill(patchx,patchy,colvar(colsel,:),'facealpha',0.25)
        
        
        px =  elecentx(noi(2)) - noilocs(1);
        
        py =  elecenty(noi(2)) - noilocs(2);
        
        
        qx = patchx(4) - patchx(2);
        
        qy =  patchy(4) - patchy(2);
        AC = 0.5 * abs((qx*py - qy*px));
        areac = areac + AC;
    end
    colvar(colsel,:) = [];
end

% Nodes on boundries

for j = 1:size(nl,1)
    colsel = randi([1,length(colvar)]);
    
    if size(nl{j} ~= 0)
        h(pcnt) = plot(nl{j}(:,2),nl{j}(:,3),'x','color',colvar(colsel,:),'markersize',8,'linewidth',1.7);
        pcnt = pcnt+1;
    end
    
    for i = 1:size(nl{j},1)
        if nodeLabels == 1
            text(nl{j}(i,2),nl{j}(i,3),num2str(nl{j}(i)),'Color','red','FontSize',14,'FontWeight','Bold')
        end
    end
    
    % Boundary Of Interest
    BOI = nl{j};
    % Gets the NOI Location
    
    % For each node on the boundary
    %     for curnode = 1:(size(BOI,1))
    % Node of interest (noi) location - BOUNDARY NODES
    
    
    for i = 1:size(el{j},1) % for each element that is connected via the noi
        
        % Element Of Interest in BOI
        EOI = el{j}(i,:);
        noilocs = Nodedata(EOI(1),2:3);
        % retrieve Node of interest
        %  NOI = sum((BOI(:,1) == EOI(1)).*(1:size(BOI,1))');
        
        % Nodes in the slected EOI
        elenodes = EOI(3:5);
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
        
        
        fill(patchx,patchy,colvar(colsel,:),'facealpha',0.25)
        % Apply infill
        px =  elecentx(EOI(2)) - noilocs(1);
        
        py =  elecenty(EOI(2)) - noilocs(2);
        
        
        qx = patchx(4) - patchx(2);
        
        qy =  patchy(4) - patchy(2);
        
        ACb = 0.5 * abs((qx*py - qy*px));
        areab = areab + ACb;
    end
    %     end
    colvar(colsel,:) = [];
    
    
end

% Internal Nodes

colsel = randi([1,length(colvar)]);
h(pcnt) = plot(ns(:,2),ns(:,3),'s','color',colvar(colsel,:),'markersize',8,'linewidth',1.7);
uistack(h(pcnt),'top')
pcnt = pcnt+1;

legvar = length(h);



for j = 1:size(ns,1)
    
    if nodeLabels == 1
        text(ns(j,2),ns(j,3),num2str(ns(j,1)),'Color','red','FontSize',14,'FontWeight','Bold')
    end
    
    % INTERNAL NODES
    NOI = ns(j,:);
    noilocs = NOI(2:3);
    % Indexes of Elements on node j
    elindex = nonzeros(sum(se(:,1) == NOI(1),2).*(1:size(se,1))');
    
    for i = 1:size(elindex,1) % for each element that is connected via the noi
        
        % Element Of Interest in BOI
        EOI = se(elindex(i),:);
        
        % Nodes in the slected EOI
        elenodes = EOI(3:5);
        
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
        fill(patchx,patchy,colvar(colsel,:),'facealpha',0.25)
        
        px =  elecentx(EOI(2)) - noilocs(1);
        
        py =  elecenty(EOI(2)) - noilocs(2);
        
        
        qx = patchx(4) - patchx(2);
        
        qy =  patchy(4) - patchy(2);
        
        ACs = 0.5 * abs((qx*py - qy*px));
        areai = areai + ACs;
        
    end
end
colvar(colsel,:) = [];

% Elements

xele = [0 0 0 0];
yele = [0 0 0 0];
for i = 1:length(Elements)
    elenodes = [Elements(i,2:4) Elements(i,2)];
    
    
    for j = 1:4
        
        z = 1:size(Nodedata,1)  == elenodes(j);
        xele(j) =  Nodedata(z,2);
        
        z = 1:size(Nodedata,1) == elenodes(j);
        yele(j) =  Nodedata(z,3);
        
        
    end
    
    plot(xele,yele,'color',[0.5 0 0.5],'linewidth',1.25)
    if ElementLabels == 1
        text(mean(xele(1:3)),mean(yele(1:3)),num2str(Elements(i)))
    end
    
end




% Border
for i = 1:length(xmin)
    xarr = [pointloc(pointtag1(i),1) pointloc(pointtag2(i),1)];
    yarr = [pointloc(pointtag1(i),2) pointloc(pointtag2(i),2)];
    
    plot(xarr,yarr,'k','linewidth',2)
end
% 
% 
% SCV
for i = 1:size(nodeymid,1)
    plot(nodexmid(i,:),nodeymid(i,:),'k','linewidth',1.25)
    
    
end
xlim([min(xmin) max(xmax)])
ylim([min(ymin) max(ymax)])
pbaspect([(max(xmax)-min(xmin)) (max(ymax)-min(ymin)) 1])

legendstrings = cell(1, length(np)+length(nl)+1);
for i = 1:(length(np)+length(nl)+1)
    
    if i <= length(np)
        legendstrings{i} = sprintf(' Corner Node %d',i);
    elseif i <= (length(np)+length(nl))
        legendstrings{i} = sprintf(' Boundary Nodes %d',i-length(np));
    else
        legendstrings{i} = sprintf(' Internal Nodes');
    end
end

% Scaled Normal Vector Plot
% SCVNormalsx1;
% For each element
for e = 1:size(Elements)
    
    % For each SCV
    for SCV = 1:3
        
        xnorm1 = [0.5*(elecentx(e)+elemidsx(e,SCV)),...
            0.5*(elecentx(e)+elemidsx(e,SCV)) + SCVNormalsx1(e,SCV)];
        
        ynorm1 = [0.5*(elecenty(e)+elemidsy(e,SCV)),...
            0.5*(elecenty(e)+elemidsy(e,SCV)) + SCVNormalsy1(e,SCV)];
        plot(xnorm1,ynorm1,'r','linewidth',1.5)
        
        xnorm1 = [0.5*(elecentx(e)+elemidsx(e,SCV)),...
            0.5*(elecentx(e)+elemidsx(e,SCV)) + SCVNormalsx2(e,SCV)];
        
        ynorm1 = [0.5*(elecenty(e)+elemidsy(e,SCV)),...
            0.5*(elecenty(e)+elemidsy(e,SCV)) + SCVNormalsy2(e,SCV)];
        plot(xnorm1,ynorm1,'b','linewidth',1.5);
    end
    
    
    
end







areacv = areac + areab + areai;
titlestr = {strcat('Total Area=  ',num2str(areacv),', Corner area=  ', num2str(areac)),...
    strcat('Boundary Area=  ', num2str(areab),', Internal Area=  ', num2str(areai))};
title(titlestr)
if legendmyguy == 1
    hlegend = legend((h(1:pcnt-1)),legendstrings,'location','eastoutside','FontSize',12);
    pbaspect([1 1 1])
end
shg



