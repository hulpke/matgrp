gap> START_TEST("matgrp basic.tst");
gap> SetAssertionLevel(0);SetInfoLevel(InfoRecog,0);
gap> gens:=[ [ [ -28, 3, 9 ], [ 0, -1, -6 ], [ 3, 0, 1 ] ],
> [ [ -1, 0, 0 ], [ -9, 1, 3 ], [ -3, 0, -1 ] ],
> [ [ 0, 0, 1 ], [ 1, 0, 9 ], [ 0, 1, 0 ] ] ];;
gap> g:=Group(List(gens,x->x*One(Integers mod 89)));
<matrix group with 3 generators>
gap> ffs:=FittingFreeLiftSetup(g);;
gap> Size(Image(ffs.factorhom));
704880
gap> Product(RelativeOrders(ffs.pcgs));
697048
gap> h:=Group(List(gens,x->x*One(Integers mod 2403)));
<matrix group with 3 generators>
gap> ffs:=FittingFreeLiftSetup(g);;
gap> Size(Image(ffs.factorhom));
704880
gap> ffs:=FittingFreeLiftSetup(h);;
gap> Size(Image(ffs.factorhom));
704880
gap> Product(RelativeOrders(ffs.pcgs));
164639949408
gap> cl:=ConjugacyClasses(g);;
gap> Collected(List(cl,Size));
[ [ 1, 1 ], [ 7920, 1 ], [ 7921, 87 ], [ 704880, 1 ], [ 712890, 87 ], 
  [ 62029440, 1 ], [ 62037272, 3916 ], [ 62734320, 174 ], [ 63447210, 3741 ] ]
gap> STOP_TEST("matgrp basic.tst",0);
