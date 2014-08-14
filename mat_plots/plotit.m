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




xgrid_bad;
xg_bad = sol;
ygrid_bad;
yg_bad = sol;
ugrid_bad;
ug_bad = sol;


figure(5)
du_bad = ulim(2) - ulim(1);
nc = 40;
vc = ulim(1):du_bad/nc:ulim(2);
contour(xg_bad, yg_bad, ug_bad, vc)
hold on
geo_contours
