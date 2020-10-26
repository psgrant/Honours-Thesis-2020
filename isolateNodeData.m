function [Nodedata, cornernodes, linenodes, surfnodes,numnodes] = isolateNodeData(nodestr,numents)



NBcell = textscan(nodestr,'%s','Delimiter','\r');
numnodes = str2num(NBcell{1}{2});
numblocks = numnodes(1);
numnodes = numnodes(4);
seed = 0;
NBcell{1}(1:2) = [];
prev = 1;
pointer = 1;
Nodes = [];
counter = 100;

cornernodes = zeros(numents(1),3);
linenodes = cell(numents(2),1);


infodata = str2num(NBcell{1}{1});

for pn = 1:sum(numents(1:3))
    
    infodata = str2num(NBcell{1}{1})
    % Corner nodes!
    
    
    
    switch infodata(1)
        
        case 0 % Points
            
            cpt = str2num(NBcell{1}{2});
            nodeloc = str2num(NBcell{1}{3});
            cornernodes(pn,:) = [cpt , nodeloc(1:2)];
            NBcell{1}(1:3) = [];
            
        case 1 % Lines
            tl = infodata(4);
            cpt = zeros(tl,4);
            
            for ln = 1:tl
                cpt(ln,:) = [str2num(NBcell{1}{ln+1}), str2num(NBcell{1}{ln+1+tl})];
            end
            linenodes{infodata(2)} = cpt(:,1:3);
            NBcell{1}(1:2*tl+1) = [];
            
        case 2 % Surface
            ts = infodata(4);
            surfnodes = zeros(ts,4);
            for sn = 1:ts
                surfnodes(sn,:) = [str2num(NBcell{1}{sn+1}), str2num(NBcell{1}{ts+1+sn})];
            end
        otherwise
            error('This should not happen. There are volumes in the defined mesh')
            
    end
    
end


Nodedata = [cornernodes];

for i = 1:length(linenodes)
    Nodedata = [Nodedata; linenodes{i}];
end

Nodedata = [Nodedata; surfnodes(:,1:3)];