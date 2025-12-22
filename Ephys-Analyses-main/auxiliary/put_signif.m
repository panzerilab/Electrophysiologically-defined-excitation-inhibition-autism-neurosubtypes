function put_signif(p_vals)

[r,c] = size(p_vals);
for i=1:r
    for j=1:c
        if i~=j
           str = getpstars(p_vals(i,j)); 
           x = j-0.1;
           y = i-0.1;
           text(x,y,str);
        end
    end
end

end

