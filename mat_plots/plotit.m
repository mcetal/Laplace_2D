close all
clear

xgrid;
xg = a;
ygrid;
yg = a;
igrid;
ig = a;
igrid_c;
igc = a;
ugrid;
ug = a;

figure(1)
mesh(xg,yg,ig)
hold on

figure(4)
mesh(xg,yg,igc)

figure(2)
geo_contours
target_points


figure(5)
du = ulim(2) - ulim(1);
nc = 40;
vc = ulim(1):du/nc:ulim(2);
contour(xg, yg, ug, vc)
hold on
geo_contours
hold on
geo_bad_contours
hold on

xgrd_plot;
xg = sol';
ygrd_plot;
yg = sol';

bad_param;

for ik = 1:K
a = zeros(NT, NR);
ugrd_bad_plot;
ug = sol((ik-1)*NR*NT + 1:ik*NR*NT);
a(:) = ug(:);


xgrid_bad_plot;
xg = sol((ik-1)*NR*NT + 1:ik*NR*NT);
b = zeros(NT, NR);
b(:) = xg(:);


ygrid_bad_plot;
yg = sol((ik-1)*NR*NT + 1:ik*NR*NT);
c = zeros(NT, NR);
c(:) = yg(:);


du_bad = ulim(2) - ulim(1);
nc = 40;
vc = ulim(1):du_bad/nc:ulim(2);
contour(b, c, a)
hold on
%geo_contours
%hold on
%geo_bad_contours

end 
