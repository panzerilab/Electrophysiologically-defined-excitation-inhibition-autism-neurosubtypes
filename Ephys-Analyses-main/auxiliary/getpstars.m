function [ pstars ] = getpstars( pvalue )
if pvalue <0.001
    pstars=' ***';
elseif pvalue<0.01
    pstars=' ** ';
elseif pvalue<0.05
    pstars='  * ';
else
    pstars='  n.s ';
end
end

