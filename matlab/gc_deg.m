function d = gc_deg(la1, lo1, la2, lo2)
% distance between two points in degrees (haversine formula)
la1 = deg2rad(la1); la2 = deg2rad(la2); dlo = deg2rad(lo2 - lo1);
a = sin((la2-la1)/2).^2 + cos(la1).*cos(la2).*sin(dlo/2).^2;
d = rad2deg(2*asin(sqrt(a)));
end
