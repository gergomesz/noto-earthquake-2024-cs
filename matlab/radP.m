function R = radP(s, d, r, az, i)
% P wave radiation pattern (Aki & Richards)
% s,d,r = strike, dip, rake, az = azimuth, i = takeoff angle
R = cosd(r).*sind(d).*sind(i).^2.*sind(2*(az-s)) - cosd(r).*cosd(d).*sind(2*i).*cosd(az-s) ...
  + sind(r).*sind(2*d).*(cosd(i).^2 - sind(i).^2.*sind(az-s).^2) + sind(r).*cosd(2*d).*sind(2*i).*sind(az-s);
end
