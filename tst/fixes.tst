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

# base point candidates from quotient modules (MSSFBPC, factor steps)
gap> F:=GF(3);; k:=6;; n:=2*k;;
gap> Reset(GlobalMersenneTwister,1);; Reset(GlobalRandomSource,1);;
gap> mk:=function(A) local M,i,j; M:=MutableCopyMat(IdentityMat(n,F));
>   for i in [1..k] do M{[2*i-1,2*i]}{[2*i-1,2*i]}:=A;
>     for j in [1..i-1] do M{[2*i-1,2*i]}{[2*j-1,2*j]}:=RandomMat(2,2,F); od;
>   od; return ImmutableMatrix(F,M); end;;
gap> G:=Group(List(GeneratorsOfGroup(SL(2,3)),mk));;
gap> opt:=rec();; GENSS_CopyDefaultOptions(GENSS,opt);;
gap> c:=MSSFBPC(G,opt,1,false:Subrecurse:=false);;
gap> ForAll([1..Length(c.points)],i->ForAll(GeneratorsOfGroup(G),x->
>   ForAll(GeneratorsOfGroup(G),y->c.ops[i](c.ops[i](c.points[i],x),y)
>     =c.ops[i](c.points[i],x*y))));
true
gap> Size(StabilizerChain(G,rec(Cand:=c,StrictlyUseCandidates:=true,
>   Reduced:=false)))=Size(StabilizerChain(G));
true
gap> STOP_TEST("matgrp fixes.tst",0);
