gap> START_TEST("matgrp fixes.tst");
gap> SetAssertionLevel(0);SetInfoLevel(InfoRecog,0);
gap> gens:=[ [ [ -28, 3, 9 ], [ 0, -1, -6 ], [ 3, 0, 1 ] ],
> [ [ -1, 0, 0 ], [ -9, 1, 3 ], [ -3, 0, -1 ] ],
> [ [ 0, 0, 1 ], [ 1, 0, 9 ], [ 0, 1, 0 ] ] ];;

# residue class ring, generators given as lists of lists
gap> h:=Group(List(gens,x->x*One(Integers mod 9)));;
gap> ffs:=FittingFreeLiftSetup(h);;
gap> IsIdenticalObj(Source(ffs.factorhom),h);
true
gap> ForAll(ffs.pcgs,x->x in h);
true
gap> Length(ConjugacyClasses(h));
13

# trivial group
gap> t:=Group(One(GL(3,5)));;
gap> Size(Image(FittingFreeLiftSetup(t).factorhom));
1
gap> STOP_TEST("matgrp fixes.tst",0);
