function surf2stl(surfhandle, filename, mode)
%SURF2STL   Write STL file from surface data.
%   SURF2STL(surfhandle, 'filename') writes a stereolithography (STL) file
%   for a surface with the given handle surfhandle.
%
%   SURF2STL(...,'mode') may be used to specify the output format.
%
%     'binary' - writes in STL binary format (default)
%     'ascii'  - writes in STL ASCII format
%
%   Example:
%        [X, Y, Z] = peaks;
%        s = surf(X, Y, Z);
%        surf2stl(s, 'test.stl');
%   See also SURF.

    if (nargin < 3)
        mode = 0;
    end
    X = surfhandle.XData; Y = surfhandle.YData; Z = surfhandle.ZData;
    [N, M] = size(X);
    fileID = fopen([filename,'.stl'],'w');
    if(~mode)
        fileID = fopen([filename,'.stl'],'wb+');
    end
    if(mode)
        fprintf(fileID,'%1s\n','solid MYSOLID');
    else
        str = sprintf('%-80s','');
        fwrite(fileID,str,'uchar');
        fwrite(fileID,0,'int32'); 
    end
    nfacets = 0;
    for n = 1:N-1
        np1 = n+1;
        for m = 1:M-1
            mp1 = m+1;
            p1 = [X(n, m), Y(n, m), Z(n, m)];
            p2 = [X(np1, m), Y(np1, m), Z(np1, m)];
            p3 = [X(np1, mp1), Y(np1, mp1), Z(np1, mp1)];
            if(~any(isnan(p1) | isnan(p2)| isnan(p3)))
                nfacets = nfacets + stlwrite(p1, p2, p3, fileID, mode);
            end
            
            p1 = [X(np1, mp1), Y(np1, mp1), Z(np1, mp1)];
            p2 = [X(n, mp1), Y(n, mp1), Z(n, mp1)];
            p3 = [X(n, m), Y(n, m), Z(n, m)];
            if(~any(isnan(p1) | isnan(p2)| isnan(p3)))
                nfacets = nfacets + stlwrite(p1, p2, p3, fileID, mode);
            end
        end
    end
    if(mode)
        fprintf(fileID,'%1s\n','endsolid MYSOLID');
    else
        fseek(fileID, 0,'bof');
        fseek(fileID,80,'bof');
        fwrite(fileID, nfacets,'int32');
    end
    fclose(fileID);
end
function n = stlwrite(p1, p2, p3, fileID, mode)
    v12 = p2 - p1; v23 = p3 - p2; v31 = p1 - p3; 
    s = (norm(v12) + norm(v23) + norm(v31))/3;
    normal = cross(v12, v23); nnorm = norm(normal); 
    normal = normal/nnorm; relheight = nnorm/s^2;
    n = 0;
    if (relheight > 1e-6)
        n = 1;
        if(mode)
            fprintf(fileID,'%14s','facet normal');
            fprintf(fileID,'%12.4e %12.4e %12.4e\n',normal);
            fprintf(fileID,'%14s\n','outer loop');
            fprintf(fileID,'%14s','vertex  ');
            fprintf(fileID,'%12.4e %12.4e %12.4e\n',p1);
            fprintf(fileID,'%14s','vertex  ');
            fprintf(fileID,'%12.4e %12.4e %12.4e\n',p2);
            fprintf(fileID,'%14s','vertex  ');
            fprintf(fileID,'%12.4e %12.4e %12.4e\n',p3);
            fprintf(fileID,'%11s\n','endloop');
            fprintf(fileID,'%10s\n','endfacet');
        else
            fwrite(fileID,normal,'float32');
            fwrite(fileID,p1,'float32');
            fwrite(fileID,p2,'float32');
            fwrite(fileID,p3,'float32');
            fwrite(fileID,0,'int16');
        end
    end
end