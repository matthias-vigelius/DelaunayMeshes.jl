import graph;

size(20cm,20cm);
pair v1 = (1,1);
pair v2 = (2,1);
pair v3 = (1.5,2);
draw(v1--v2, L="ea=1", Arrow);
draw(v2--v3, L="eb=5", Arrow);
draw(v3--v1, L="ec=9", Arrow);

pair v4 = (1.375, 1.25);
pair v5 = (1.5,1.5);

draw(v1--v4, L="ed=13", Arrow);
draw(v2--v4, L="ee=17", Arrow);
draw(v3--v4, L="ef=21", Arrow);

draw(v4--v5, L="eg=25", Arrow);
draw(v2--v5, L="eh=29", Arrow);
draw(v3--v5, L="ei=33", Arrow);

pair f1 = (v1+v2+v4)/3.;
pair f3 = (v2+v4+v5)/3.;
pair f4 = (v1+v3+v4)/3.;
pair f5 = (v2+v3+v5)/3.;
pair f6 = (v3+v4+v5)/3.;
label("1", f1);
label("3", f3);
label("4", f4);
label("5", f5);
label("6", f6);

draw(v1--v5, L="ei'", dashed, Arrow); 
