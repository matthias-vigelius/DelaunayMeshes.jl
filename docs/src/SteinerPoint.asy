import graph;

size(20cm,20cm);
pair cc = (0.,0.);
real a0 = pi/24.;
real ar = pi/4.;
real r0 = 2.;
real beta = sqrt(2.);

path outcircle=shift(cc)*scale(r0)*unitcircle;
draw(outcircle, dotted);


pair p = r0 * (cos(pi - a0), sin(pi - a0));
pair q = r0 * (cos(pi + a0), sin(pi + a0));
pair r = r0 * (cos(ar), sin(ar));
draw(p--q);
draw(q--r);
draw(r--p);

label("$p$", p, NW);
label("$q$", q, SW);
label("$r$", r, E);

pair bs = 0.5 * (p + q);
draw(bs--cc, dashed+blue);

real lpq = length(p - q);
real ccr = beta * lpq;
real dobsoc = sqrt(ccr*ccr - (lpq/2.)*(lpq/2.) );
pair occ = bs + dobsoc * unit(cc - bs);
path occircle = shift(occ)*scale(ccr)*unitcircle;
draw (occircle, dashed+red);

pair ocr = occ + ccr * unit(cc - bs);
draw(p--ocr, red);
draw(q--ocr, red);
label("$c$", ocr, SE);
