function [elements, BoundaryElements] = getElementInfo(elestr,numents)


Elecell = textscan(elestr,'%s','delimiter','\r');
totelements = str2num(Elecell{1}{2});
totelements = totelements(4);
Elecell{1}(1:2) = [];
pointer = 0;
elements = zeros(1,4);
BoundaryElements = cell(numents(2),1);

for surfs = 1:numents(3)
    arenodes2d = 0;
    
    while arenodes2d == 0
        
        teststr = str2num(Elecell{1}{1});
        strlen = length(teststr);
        strindex = teststr(1);
        
        
        % If there current row is a BLOCK HEADER
        if strlen == 4
            
            % See if the current element is a
            switch strindex
                
                
                case 1
                    % Pre allocate space in the cell array
                    BoundaryElements{teststr(2)} = zeros(teststr(4),3);
                    
                    
                    % Loop through all 1D boundary elements converting the string into
                    % an array
                    for i = 1:teststr(4)
                        BoundaryElements{teststr(2)}(i,:) = str2num(Elecell{1}{i+1});
                    end
                    
                    % Remove elements from array
                    Elecell{1}(1:(i+1),:) = [];
                    
                    % Tests to see if the element is 2d (3 node element)
                case  2
                    numelements = teststr(4);
                    Elecell{1}(1) = [];
                    arenodes2d = 1;
                otherwise
                   Elecell{1}(1) = []; 
                    
            end
        else
            Elecell{1}(1) = [];
        end
    end
    
    for i = 1:numelements
        elements(i+pointer,:) = str2num(Elecell{1}{i});
    end
    
    pointer =  pointer + numelements;
    
end
elements(:,1) = elements(:,1) - elements(1,1) + 1;