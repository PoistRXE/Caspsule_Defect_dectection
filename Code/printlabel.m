function printlabel(boundingBox,flag)
    if flag~=0
        rectangle('Position',boundingBox,'edgecolor','r','LineWidth', 1.000);
    else
        rectangle('Position',boundingBox,'edgecolor','g','LineWidth', 1.000);
    end
end