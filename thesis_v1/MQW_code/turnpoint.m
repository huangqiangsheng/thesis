function [indx] = turnpoint(y)
dy = y(2:end) - y(1:end-1);
dy = dy./(abs(dy));
ddy = dy(2:end)+dy(1:end-1);
indx = find((ddy == 0 | ddy == inf));
end