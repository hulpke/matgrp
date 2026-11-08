#############################################################################
##
#W  recograt.gi                 matgrp package               Alexander Hulpke
##
##
#Y  Copyright (C)  2014-18, Alexander Hulpke
##
##  basic setup for matrix fitting free.
##

SetInfoLevel(InfoRecog,0); # recog will print status messages otherwise

# library change
if not IsBound(RUN_IN_GGMBI) then RUN_IN_GGMBI:=fail;fi;


# newer recog has some changes to these internal APIs
if not IsBound(ImageRecogNode) then
  ImageRecogNode := RIFac;
  KernelRecogNode := RIKer;
fi;

OnSubmoduleCosets:=function(cset,g)
  return [SiftedVector(cset[2],cset[1]*g),cset[2]];
end;

MakeSubmoduleCosetAction:=function(basis)
  return function(vec,g)
    return SiftedVector(basis,vec*g);
  end;
end;

MakeSubmoduleColineAction:=function(basis)
  return function(vec,g)
  local
    c;  # position of the first nonzero entry of the sifted vector, used to normalize (scale) the coline representative
    vec:=SiftedVector(basis,vec*g);
    c:=PositionNonZero(vec);
    if c<=Length(vec) then
      vec := Inverse( vec[c] ) * vec;
    fi;
    return vec;
  end;
end;

FUNCSPACEHASH:=[];
# mod to genss -- use submodules and quotients for base points

MSSFBPC:=function( grp, opt, ii, parentS ) # owf fct to call easily
local
  # -- local helper function --
  trysel,     # helper: test a coordinate range as a subquotient module and, if usable, add base-point candidates from it
  # -- group / field setup --
  limit,      # orbit-length limit; spaces whose orbits would exceed this are considered too big and skipped
  F,          # the field over which the matrix group `grp` is defined
  cand,       # record (fields ops/points/used) accumulating the base-point candidates that get returned
  # -- module and adapted basis --
  mo,         # WARNING: multi-role -- (1) the GModule of `grp`, then (2) the list of generator matrices rewritten in
              #          the adapted basis, then (3) an action function taken from opt.PCand.ops in the parent fallback
  cs,         # WARNING: multi-role -- (1) the group generators, then (2) the bases of the composition series,
              #          then (3) an orbit image point during the parent-candidate orbit enumeration
  dim,        # dimension of the module `mo`
  dims,       # list of the dimensions of the composition-series terms
  bas,        # new basis of the underlying space, adapted to the composition series
  basc,       # accumulator of subspace+factorspace vectors while `bas` is being built up
  basinv,     # inverse of the basis-change matrix `bas`
  lastdim,    # dimension of the submodule processed so far
  nb,         # WARNING: two roles -- (1) a BaseSteinitzVectors basis-extension record while building `bas`,
              #          and (2) a freshly constructed base-point (vector or subspace) inside trysel
  # -- subquotient action (used inside trysel) --
  fcand,      # base-point candidates found (recursively) for the current subquotient factor
  fct,        # action function (submodule coset or coline action) built for a submodule
  submodule,  # semi-echelonized basis of the submodule used for the coset/coline action
  sel,        # WARNING: two roles -- (1) the coordinate range of the current subquotient, and
              #          (2) the list of indices of parent candidates that turned out to have short orbits
  # -- parent-candidate orbit fallback --
  orb,        # orbit currently being enumerated
  orbs,       # the same orbit stored as a Set for fast membership tests
  # -- loop counters --
  i,          # loop index over parent candidate points
  j,          # loop index over composition-series terms / dimension steps
  k;          # loop index over generators (and over factor candidates inside trysel)

  trysel:=function(recsub,recfac)
  local
    lgens;  # the module generators restricted to the selected coordinate range `sel` (the subquotient action)
    # nor the trivial action
    if ForAny(mo,x->Order(x{sel}{sel})>1) then
    Info(InfoFFMat,2,"range ",sel," have ",Length(cand.points));
      lgens:=List(mo,x->x{sel}{sel});
      fcand:=FindBasePointCandidates(Group(lgens),opt,ii,
               false:Subrecurse:=recsub,Facrecurse:=recfac);
      Info( InfoGenSS, 3, "Subfactor module of range ",sel,", ",Length(fcand.ops),
            " candidates");
      for k in [1..Length(fcand.ops)] do
        if ForAll(lgens,x->fcand.ops[k](fcand.points[k],x)=fcand.points[k]) then
          Info(InfoFFMat,2,"Ignoring fixed element for base");

        elif fcand.ops[k]=OnRight or fcand.ops[k]=OnPoints or fcand.ops[k]=OnLines then
          nb:=fcand.points[k]*bas{sel};
          if lastdim=0 then
            # proper subspace -- just vectors
            Add(cand.points,nb);
            Add(cand.ops,fcand.ops[k]);
          elif true then
            # # action on cosets
            submodule:=SemiEchelonBasis(VectorSpace(F,bas{[1..lastdim]},Zero(bas[1])));
            nb:=SiftedVector(submodule,nb);
            if fcand.ops[k]=OnLines then
              fct:=MakeSubmoduleColineAction(submodule);
              Add(FUNCSPACEHASH,[fct,submodule]);
            else
              fct:=MakeSubmoduleCosetAction(submodule);
              Add(FUNCSPACEHASH,[fct,submodule]);
            fi;
            Add(cand.points,nb);
            Add(cand.ops,fct);
          elif Length(sel)=1 then
            # TODO: 1-dim factor -- need to do cosets
            Info(InfoWarning,1,"Case not yet implemented");
          else
            # subfactor -- take subspace preimage
            nb:=OnSubspacesByCanonicalBasis(Concatenation(bas{[1..lastdim]},[nb]),
                  One(grp));
            Add(cand.points,nb);
            Add(cand.ops,OnSubspacesByCanonicalBasis);
          fi;
        elif ForAny(FUNCSPACEHASH,x->x[1]=fcand.ops[k]) then
          fct:=First(FUNCSPACEHASH,x->x[1]=fcand.ops[k]);
          submodule:=SemiEchelonBasis(VectorSpace(F,
               Concatenation(bas{[1..lastdim]},
                 BasisVectors(fct[2])*bas{sel})));
          Add(cand.points,fcand.points[k]*bas{sel});
          fct:=MakeSubmoduleCosetAction(submodule);
          Add(FUNCSPACEHASH,[fct,submodule]);
          Add(cand.ops,fct);
          Info(InfoFFMat,2,"ACTPOP");
        elif fcand.ops[k]=OnSubspacesByCanonicalBasis then
          nb:=fcand.points[k]*bas{sel};
          nb:=OnSubspacesByCanonicalBasis(Concatenation(bas{[1..lastdim]},nb),
                One(grp));
          Add(cand.points,nb);
          Add(cand.ops,OnSubspacesByCanonicalBasis);
        else
          Info(InfoWarning,1,"Action not recognized");
        fi;
      od;
      return true;
    else
      return false;
    fi;
  end;

  if IsBound(opt.VeryShortOrbLimit) then
    limit:=2*opt.VeryShortOrbLimit;
  else
    limit:=10^4;
  fi;

  F := DefaultFieldOfMatrixGroup(grp);

  # don't bother in small cases
  if Size(F)^Length(One(grp))<=opt.ShortOrbitsOrbLimit then
    TryNextMethod();
  fi;

  cs:=GeneratorsOfGroup(grp);
  if ForAny(cs,IsObjWithMemory) then
    cs:=Concatenation(Filtered(cs,x->not IsObjWithMemory(x)),
           List(Filtered(cs,x->IsObjWithMemory(x)),
                x->x!.el));
  fi;

  cand:=rec(ops:=[],points:=[],used:=0);
  mo:=GModuleByMats(cs,F);
  dim:=mo.dimension;
  if MTX.IsIrreducible(mo) or Size(F)^dim<=limit then
    TryNextMethod();
  fi;

  # build new basis corresponding to comp ser.
  cs:=MTX.BasesCompositionSeries(mo);
  Info(InfoFFMat,2,"dims=",List(cs,Length));
  dims:=List(cs,Length);
  bas:=[];
  basc:=[];
  for j in [2..Length(cs)] do
    nb:=BaseSteinitzVectors(cs[j],basc);
    Append(bas,nb.factorspace);
    basc:=Concatenation(nb.subspace,nb.factorspace);
    Sort(basc,function(a,b) return PositionNonZero(a)<PositionNonZero(b);end);
  od;
  basinv:=bas^-1;

  mo:=List(mo.generators,x->bas*x*basinv);

  # now step up in sizes indicated by short orbit lengths
  lastdim:=0;
  j:=2;
  if ValueOption("Subrecurse")<>false then
    while j<=Length(dims) and Length(cand.points)<=5 do
      # don't bother is the space is too small
      if (j=Length(dims) and lastdim>0) or
        (j<Length(dims) and Size(F)^(dims[j+1]-lastdim)>limit) then
        sel:=[lastdim+1..dims[j]];
        if trysel(false,true) then
          # we tried a space
          lastdim:=dims[j];
        fi;

      fi;

      j:=j+1;
    od;
  fi;

  if lastdim=0 and ValueOption("Facrecurse")<>false then

    # all but last submodule was trivial -- step down on factors
    j:=Length(dims)-1;
    while j>0 and Length(cand.points)<=5 do
      lastdim:=dims[j]-1;
      if (j>1 and Size(F)^(dim-dims[j-1])>limit) then
        sel:=[lastdim+1..dim];
        if trysel(false,false) then
          j:=0;
        fi;
      fi;
      j:=j-1;
    od;

  fi;

  if Length(cand.points)>0 then return cand; fi;

  # refer to parent
  if IsBound(opt.PCand) then
    # try which ones are short
    sel:=[];
    for i in [1..Length(opt.PCand.points)] do
      orb:=[opt.PCand.points[i]];
      mo:=opt.PCand.ops[i];
      orbs:=Set(orb); # as short, set is fine.
      j:=1;
      while j<=Length(orb) and Length(orb)<=limit do
        for k in GeneratorsOfGroup(grp) do
          cs:=mo(orb[j],k);
          if not cs in orbs then
            Add(orb,cs);
            AddSet(orbs,cs);
          fi;
        od;
        j:=j+1;
      od;
      if Length(orb)<=limit then
        Add(sel,i);
        Add(cand.points,orb[1]);
        Add(cand.ops,mo);
      fi;
    od;
    Info(InfoFFMat,2,"Selected ",Length(sel)," of ",Length(opt.PCand.points)," group basis vectors");
    if Length(cand.points)>0 then return cand; fi;
  fi;

  TryNextMethod();
end;

InstallMethod(FindBasePointCandidates,
  "for reducible matrix group over a FF, use submodules and quotients",
  [ IsGroup and IsMatrixGroup and IsFinite, IsRecord, IsInt, IsObject ], 21,
  MSSFBPC);

GetInformationFromRecog:=function(recog)
local
  # -- local helper function --
  treerecurse,  # helper: recursively traverse the recognition tree, filling the per-factor arrays below
  # -- running indices / counters --
  n,            # running index of the current composition factor (also ends up as the total number of factors)
  genum,        # running count of nice generators consumed so far (generator-number offset)
  fn,           # running index into `furtherhom`
  extranum,     # running index/counter for `extragens`
  # -- per-factor result arrays (indexed by factor number n) --
  leafs,        # for each factor: the recognition (leaf) node, or false
  homs,         # for each factor: the homomorphism path `h` leading to it
  factors,      # for each factor: the factor group, or false
  sz,           # for each factor: the order of the factor
  leafgens,     # for each factor: the number of nice generators at the leaf, or false
  niceranges,   # for each factor: the range of nice-generator numbers, or false
  minranges,    # for each factor: the reduced generator-number range (minimal generating subset)
  mingens,      # for each factor: the indices of the chosen minimal generators
  permap,       # for each factor: a map to a permutation representation (or identity mapping)
  furtherhom,   # additional homomorphisms used when a non-simple factor has to be split
  extragens,    # extra generators added to compatibilize generating sets across a split
  # -- whole-tree data --
  allgens,      # the full list of nice generators of the whole recognition tree
  # -- variables shared into treerecurse via the closure --
  g,            # the current factor group being worked on (matrix group, then its permutation image)
  gd,           # the perfect residuum / smaller-degree image of `g`, used to improve the degree
  orbit,        # union of stabilizer-chain orbits, used as the domain of a nice action homomorphism
  cnt,          # counter for the random-vector search when building a projective action
  extra,        # list of extra-generator indices contributed to the current factor
  nicegp,       # the group generated by the nice-generator images `ng` (used for factorization into words)
  map;          # epimorphism from a free group onto `nicegp`, for decomposing elements into words

  # the main worker function -- traverse the tree.
  treerecurse:=function(r,h)
  # NB: g, nicegp, orbit, gd, cnt and map are NOT declared here -- they are locals of the
  #     enclosing GetInformationFromRecog and are shared into this function via the closure.
  local
    # -- current-node data --
    nicehom,  # the nice homomorphism from the current node's group to a permutation representation
    stbc,     # base/stabilizer-chain information extracted from the node
    act,      # the action (op) taken from the stabilizer chain
    v,        # a vector/line used to build the orbit for a projective action homomorphism
    dom,      # domain (orbit of lines) for the projective action homomorphism
    # -- factor splitting --
    cs,       # the composition series of the current group `g`
    ng,       # the list of nice-generator images under `nicehom`
    hom,      # NaturalHomomorphismByNormalSubgroup between successive composition factors
    u,        # subgroup built up by closing over generators, to find which generators are needed
    ngi,      # image of a single nice generator under `nicehom` (for the minimal generating set)
    x,        # a generator/element currently being tested or factorized
    xf,       # factorization of `x` as a word in the nice generators
    # -- range / selection --
    f,        # WARNING: two roles -- (1) the number of nice generators at a leaf, and
              #          (2) the image recognition node when recursing into the tree
    k,        # WARNING: two roles -- (1) the range [genum+1..genum+f] of nice-generator numbers, and
              #          (2) the kernel recognition node when recursing into the tree
    sel,      # indices of the selected (needed) generators for the current factor
    # -- loop counters --
    i,        # loop index (over composition series and over nice generators)
    j;        # loop index (over generators)

    if IsLeaf(r) then
      f:=Length(NiceGens(r));
      Info(InfoFFMat,2,"Leaf",Size(r)," ",genum," ",f);
      k:=[genum+1..genum+f];
      if Size(r)=1 then
        Info(InfoFFMat,2,"ignoring trivial factor");
      elif Length(Set(Factors(Size(r))))<3 then
        Info(InfoFFMat,2,"ignoring solvable factor of order ",Size(r));

        # store info
        n:=n+Length(Factors(Size(r)));
        SetIsSolvableGroup(Grp(r),true);
        permap[n]:=IdentityMapping(Grp(r));
        leafs[n]:=r;
        homs[n]:=h;
        factors[n]:=false;
        sz[n]:=Size(r);
        leafgens[n]:=f;

        sel:=[1..Length(NiceGens(r))];
        mingens[n]:=sel;
        niceranges[n]:=k;
        minranges[n]:=k{sel};

      else
        nicehom:=false;
        if IsPermGroup(Grp(r)) then
          nicehom:=IdentityMapping(Grp(r));
        elif IsBound(r!.stabilizerchain) then
          #use ActionOnOrbit?
          stbc:=BaseStabilizerChain(r!.stabilizerchain);
          act:=stbc.ops[1];
          if ForAll(stbc.ops,x->x=act) then
            # all the same action -- union
            orbit:=[];
            for j in [1..Length(stbc.points)] do
              if not stbc.points[j] in orbit then
                orbit:=Union(orbit,Orbit(Grp(r),stbc.points[j],stbc.ops[j]));
              fi;
            od;
            if Length(orbit)>1 and Length(orbit)<50000 then
              nicehom:=ActionHomomorphism(Grp(r),orbit,stbc.ops[1],"surjective");
              Info(InfoFFMat,2,"Got degree ",Length(orbit));
            fi;
          fi;
        fi;


        if nicehom=false then
# Hasfhmethsel, success method: StabilizerChain

          if IsBound(r!.projective) and r!.projective then
            g:=Grp(r);
            v:=OnLines(g.1[1],One(g));
            dom:=ShallowCopy(Orbit(g,v,OnLines));
            nicehom:=ActionHomomorphism(g,dom,OnLines,"surjective");
            while Size(Image(nicehom))<Size(r) do
              cnt:=0;
              repeat
                v:=OnLines(Random(DefaultFieldOfMatrixGroup(g)^Length(One(g))),
                           One(g));
                cnt:=cnt+1;
                if cnt>1000 then Error("no vector found");fi;
              until not v in dom;
              Append(dom,Orbit(g,v,OnLines));
              nicehom:=ActionHomomorphism(g,dom,OnLines,"surjective");
            od;
          else
            nicehom:=IsomorphismPermGroup(Grp(r));
          fi;
        fi;
        g:=Image(nicehom);
        if Size(g)<>Size(r) then
          Error("some discrepancy happened");
        fi;

        #test!!
        # .isknownsimple, .isknownalmostsimple
        if IsBound(r!.comment) and r!.comment[1]<>'_' then
          SetIsSimpleGroup(g,true);
          gd:=g;
        elif not IsSolvableGroup(g) then
          gd:=PerfectResiduum(g);
        else
          gd:=fail;
        fi;

        if gd<>fail then
          if NrMovedPoints(g)> SufficientlySmallDegreeSimpleGroupOrder(Size(gd))
            then
              nicehom:=nicehom*SmallerDegreePermutationRepresentation(g);
              gd:=Image(nicehom);
              Info(InfoFFMat,2,"Improved degree ",NrMovedPoints(g),"->",
                NrMovedPoints(gd));
              if HasIsSimpleGroup(g) and IsSimpleGroup(g) then
                SetIsSimpleGroup(gd,true);
                SetIsSolvableGroup(gd,false);
              fi;
              g:=gd;
          fi;
        fi;

        ## indicate unsolvable
        #nsol:=IsSimpleGroup(g) or
        #  (Length(Set(Factors(Size(r))))>2 and not IsSolvableGroup(g));


        if IsSolvableGroup(g) then
          Info(InfoFFMat,2,"Ignoring solvable factor of order ",Size(g));
            n:=n+Length(Factors(Size(g)));
        elif not IsSimpleGroup(g) then
          Info(InfoFFMat,2,"doing size",Size(g)," ",IsSimpleGroup(g));

          # not simple -- split
          cs:=CompositionSeries(g);
          ng:=List(NiceGens(r),x->ImagesRepresentative(nicehom,x));
          for i in [2..Length(cs)] do
            extra:=[];
            n:=n+1;
            # test generators
            u:=cs[i];
            sel:=[];
            for j in [1..f] do
              x:=ng[j];
              if x in cs[i-1] and not x in u then
                Add(sel,j);
                u:=ClosureSubgroup(u,x);
              fi;
            od;

            nicegp:=Group(ng);
            map:=EpimorphismFromFreeGroup(nicegp);

            if Size(u)<Size(cs[i-1]) then
              Info(InfoFFMat,2,"cannot compatibilize generators -- add extras");
              while Size(u)<Size(cs[i-1]) do
                x:=First(GeneratorsOfGroup(cs[i-1]),x->not x in u);
                u:=ClosureSubgroup(u,x);

                # decompose into word
                xf:=Factorization(nicegp,x);
                if xf=fail then Error("factorization error");fi;
                x:=MappedWord(xf,MappingGeneratorsImages(map)[1],allgens{k});

                extragens[extranum]:=x;
                Add(extra,extranum);
                extranum:=extranum+1;
              od;
            fi;
            hom:=NaturalHomomorphismByNormalSubgroup(cs[i-1],cs[i]);
            if not IsAbelian(Image(hom)) then
              permap[n]:=IsomorphismPermGroup(Image(hom));
            fi;

            fn:=fn+1;
            if Size(cs[i])=1 then
              furtherhom[fn]:=nicehom;
            else
              furtherhom[fn]:=nicehom*hom;
            fi;
            leafs[n]:=false;
            homs[n]:=Concatenation(h,[fn]);
            factors[n]:=Image(hom);
            sz[n]:=Size(Image(hom));
            leafgens[n]:=false;
            niceranges[n]:=false;
            minranges[n]:=Concatenation(k{sel},extra);
          od;

        else
          # found a simple case
          n:=n+1;

          permap[n]:=nicehom;
          # get smaller generating set
          sel:=[];
          u:=TrivialSubgroup(Image(nicehom));
          for i in [1..Length(NiceGens(r))] do
            ngi:=ImagesRepresentative(nicehom,NiceGens(r)[i]);
            if not ngi in u then
              u:=ClosureGroup(u,ngi);
              Add(sel,i);
            fi;
          od;

          leafs[n]:=r;
          homs[n]:=h;
          factors[n]:=g;
          sz[n]:=Size(r);
          leafgens[n]:=f;

          mingens[n]:=sel;
          niceranges[n]:=k;
          minranges[n]:=k{sel};

        fi;

      fi;
      genum:=genum+f;
    else
      #hom:=Homom(r);
      f:=ImageRecogNode(r);
      treerecurse(f,Concatenation(h,[1]));
      k:=KernelRecogNode(r);
      if k<>fail then
        treerecurse(k,Concatenation(h,[0]));
      fi;
    fi;

  end;

  extragens:=[];
  furtherhom:=[];
  fn:=2;
  n:=0;
  genum:=0;
  factors:=[];
  homs:=[];
  leafgens:=[];
  leafs:=[];
  niceranges:=[];
  minranges:=[];
  mingens:=[];
  sz:=[];
  permap:=[];
  allgens:=NiceGens(recog);

  extranum:=Length(allgens)+1;
  treerecurse(recog,[]);

  return rec(
    group:=Grp(recog),
    n:=n,
    genum:=Length(allgens),
    recog:=recog,
    leafs:=leafs,
    factors:=factors,
    furtherhom:=furtherhom,
    sz:=sz,
    homs:=homs,
    leafgens:=leafgens,
    minranges:=minranges,
    mingens:=mingens,
    extragens:=extragens,
    permap:=permap,
    niceranges:=niceranges);
end;

# assume that hom i can be applied
CSIImageHomNr:=function(csi,n,x)
local
  r,  # the recognition node, walked down the tree along the path
  h,  # the homomorphism path csi.homs[n] (a list of 0/1/other step codes)
  i;  # loop index over the path entries of `h`
  r:=csi.recog;
  h:=csi.homs[n];
  for i in h do
    if i=0 then
      r:=KernelRecogNode(r);
    elif i=1 then
      x:=ImageElm(Homom(r),x);
      r:=ImageRecogNode(r);
    else
      x:=ImageElm(csi.furtherhom[i],x);
    fi;
  od;
  return x;
end;

# since nice ranges can contain negatives, we cannot simply sublist
CSINiceGens:=function(csi,a)
  if IsList(a) then
    return List(a,x->CSINiceGens(csi,x));
  elif a<=csi.genum then
    return NiceGens(csi.recog)[a];
  else
    return csi.extragens[a];
  fi;
end;

CSIDepthElm:=function(csi,x)
local
  i;  # index of the layer currently being tested; the first non-trivial layer is returned as the depth
  i:=1;
  while i<=csi.n and
   (not IsBound(csi.leafs[i]) or # early skipped multiple solvable factor
   (
   (IsBool(csi.leafs[i]) or not (IsBound(csi.leafs[i]!.projective) and csi.leafs[i]!.projective)
     or csi.leafs[i]!.projective=false)
   and Order(CSIImageHomNr(csi,i,x))=1)
     or
    (not IsBool(csi.leafs[i]) and (IsBound(csi.leafs[i]!.projective) and csi.leafs[i]!.projective) and
    Order(ImagesRepresentative(csi.permap[i],CSIImageHomNr(csi,i,x)))=1)) do
    i:=i+1;
  od;
  return i;
end;

# to test whether an element act trivially projectively, it is not enough to
# test on base, but also need to have sum of base
# elm must be an element that is already image
CSIProjectiveBases:=function(csi,i,elm)
local
  m;  # the projective test base: a copy of the identity rows together with their sum, cached in csi.projbases[i]
  if not IsBound(csi.projbases) then
    csi.projbases:=[];
  fi;
  if not IsBound(csi.projbases[i]) then
    m:=ShallowCopy(One(elm));
    Add(m,Sum(m));
    csi.projbases[i]:=m;
  fi;
  return csi.projbases[i];
end;

CSIDepthElm:=function(csi,x)
local
  i,     # index of the layer currently being tested; the first non-trivial layer is returned as the depth
  ximg;  # the image of `x` under the layer homomorphism (used for the projective triviality test)
  i:=1;
  while i<=csi.n do
    if not IsBound(csi.leafs[i]) or
      ((IsBool(csi.leafs[i]) or not (IsBound(csi.leafs[i]!.projective)  and csi.leafs[i]!.projective)
        or csi.leafs[i]!.projective=false)
        and Order(CSIImageHomNr(csi,i,x))=1) then
      i:=i+1;
    elif (not IsBool(csi.leafs[i]) and IsBound(csi.leafs[i]!.projective) and
        csi.leafs[i]!.projective) then
      ximg:=CSIImageHomNr(csi,i,x);
      if ForAll(CSIProjectiveBases(csi,i,ximg),z->OnLines(z,ximg)=z) then
        i:=i+1;
      else
        return i;
      fi;
    else
      return i;
    fi;
  od;
  return i;
end;

CSIAelement:=function(a,localgens,l)
local
  limg,  # images of the tuple `l` under the automorphism-group map a[2]
  rep;   # element of a[1] mapping `localgens` to `limg` by conjugation (RepresentativeAction result)
  limg:=List(l,x->Image(a[2],x));
  rep:= RepresentativeAction(a[1],localgens,limg,OnTuples);
  if rep=fail then
    Error("no rep");
  fi;
  return rep;
end;

FindAct:=function(csi)
local
  # -- setup / inputs --
  csinice,      # the nice generators of the recognition tree, NiceGens(csi.recog)  (csi is the parameter)
  dci,          # the decomposition-info record being built up and returned inside the homomorphism
  n,            # the number of layers, csi.n
  # -- pool bookkeeping (grouping isomorphic non-solvable factors) --
  pools,        # list of pools; each pool is a list of layer indices sharing an isomorphic simple factor
  pool,         # the current pool being processed, pools[i]
  poolnum,      # index of the pool the current layer belongs to
  poolimggens,  # per-pool generator images inside the automorphism-representing group
  isoms,        # per-layer isomorphism from the pool's representative group to that layer's group
  isom,         # the isomorphism from the current group `gp` to the pool's normal form
  diso,         # isomorphism to the canonical form for the image at depth `d`
  conj,         # conjugator used when joining two pools
  auts,         # per-pool list of automorphisms induced by the acting generators
  doesaut,      # per-generator record of which layers that generator acts on as an automorphism
  perms,        # per-generator record of how a generator permutes the components of a pool
  # -- generators and images --
  goodgens,     # per-layer canonical generators chosen for each factor
  lgens,        # the good generators for the current layer i
  conjgens,     # the conjugates l^k of the current generators
  genims,       # generator images in the depth-`d` permutation group
  genimgs,      # per-generator, per-layer list of generator images
  gens,         # generator images of `lgens` in the permutation group `gp`
  gp,           # the permutation group Image(csi.permap[i])
  dgp,          # the image group at depth d, Image(csi.permap[d])
  hom,          # homomorphism between the factor images / candidate automorphism
  imgdepth,     # the common depth of the conjugated generator images (consistency check)
  # -- process queue --
  process,      # work-list (queue) of layer indices still to be processed
  ii,           # index into `process` (the main processing loop counter)
  kn,           # index of the conjugating nice generator
  # -- wreath / direct-product construction --
  aels,         # per-pool automorphism-representing-group data
  localgens,    # generators mapped into the automorphism-representing group
  biggens,      # the nice generators used as the "big" generators of the wreath image
  wrimages,     # wreath-product images of the big generators
  img,          # wreath-product image accumulated for one generator
  wrs,          # list of the wreath products, one per pool
  wrsoc,        # per-pool socle subgroups of the wreath products
  dfgens,       # per-pool big-generator lists (direct-factor generators)
  dfimgs,       # per-pool wreath images (direct-factor images)
  dfs,          # the direct factors of the fitting-free socle
  e,            # WARNING: two roles -- (1) the list of embeddings of a wreath product, and
                #          (2) a single direct-product embedding dci.embeddings[i]
  m,            # length of the current pool (number of direct factors)
  w,            # WARNING: two roles -- (1) a WreathProduct group, and (2) a trivial-subgroup accumulator
                #          for the socle image
  # -- heavily reused scratch variables --
  a,            # WARNING: multi-role -- a pool index (integer), the group generators, an automorphism-
                #          representing-group datum, an image group, a matrix dimension, ... (reassigned throughout)
  b,            # WARNING: multi-role -- a list of identity images, a symmetric group, a trivial-subgroup
                #          accumulator, and finally the image subgroup
  c,            # a conjugate element l^k, and elsewhere a minranges range list
  d,            # WARNING: multi-role -- a layer depth (integer), a DirectProduct group, a wreath-image element,
                #          and a permutation (reassigned throughout)
  x,            # loop variable / element
  # -- loop counters --
  i,            # loop index
  j,            # loop index (over earlier layers, over generators, ...)
  k,            # the acting nice generator element, CSINiceGens(csi,kn)
  l;            # loop variable (over lgens and over pool elements)


  # get a list of generator images and construct the appropriate element of
  # a[1] inducing these generator images by conjugation
  aels:=[];

  # dci is decomposition information needed for restriction to proper
  # subgroups
  dci:=rec(csi:=csi);


  # isoms[i] An isomorphism from the first group in its pool to group i

  goodgens:=[];
  poolimggens:=[];
  genimgs:=List([1..Maximum(Union(csi.minranges))],x->[]);
  n:=csi.n;
  csinice:=NiceGens(csi.recog);
  # Unused: act:=[];
  auts:=[];
  pools:=[];
  perms:=[];
  isoms:=[];
  n:=csi.n;

  doesaut:=[];

  process:=[n];
  ii:=1;
  while ii<=n do
    # do we need to feed in another layer?
    if Length(process)<ii then
      Add(process,First([n,n-1..1],x->not x in process));
    fi;
    i:=process[ii];

    if IsBound(csi.permap[i]) and not IsSolvableGroup(Image(csi.permap[i])) then

      if not IsBound(goodgens[i]) then
        # no generators set yet -- just take new ones
        #lgens:=csinice{csi.minranges[i]};
        lgens:=CSINiceGens(csi,csi.minranges[i]);
        goodgens[i]:=lgens;
      else
        lgens:=goodgens[i];
      fi;

      gp:=Image(csi.permap[i]);
      # Unused: act[i]:=[];
      gens:=List(lgens,x->ImagesRepresentative(csi.permap[i],CSIImageHomNr(csi,i,x)));

#      for x in [1..Length(gens)] do
#       j:=CSIImageHomNr(csi,i,csinice[csi.minranges[i][x]]);
#       k:=Image(csi.permap[i],j);
#       if k<>gens[x] then
#         Error("GENS");
#       fi;
#      od;

      # nonsolvable factor -- Is it in pools?
      poolnum:=First([1..Length(pools)],x->i in pools[x]);
      # pools will always be joined TO the current one, so isom will not
      # have to change on the way.
      if poolnum=fail then
        Add(pools,[i]);
        poolnum:=Length(pools);
        auts[poolnum]:=[];
        isoms[i]:=false; # representative
        isom:=IdentityMapping(gp);
        Add(poolimggens,gens);
      elif i<>pools[poolnum][1] then
        isom:=isoms[i]; # iso
      else
        isom:=IdentityMapping(gp); # first one -- iso is trivial
      fi;

      #find out what the previous factors do on it
      for j in Filtered([1..i-1],x->IsBound(csi.minranges[x])) do
        Info(InfoFFMat,2,"Try ",i," mapped by ",j);

        for kn in csi.minranges[j] do
          #k:=csinice[kn];
          k:=CSINiceGens(csi,kn);
          if not IsBound(perms[kn]) then
            perms[kn]:=[];
            genimgs[kn]:=[];
            doesaut[kn]:=[];
          fi;

          conjgens:=[];
          genims:=[];
          imgdepth:=fail;
          # calculate images
          for l in lgens do
            c:=l^k; # conjugate
            Add(conjgens,c);
            d:=CSIDepthElm(csi,c);
            if imgdepth=fail then
              # new image -- isomorphism
              imgdepth:=d;
              perms[kn][i]:=d; # record permutations
            elif imgdepth<>d then
              # this cannot happen
              Error("incompatible depths");
            fi;
            Add(genims,ImagesRepresentative(csi.permap[d],CSIImageHomNr(csi,d,c)));
          od;

          dgp:=Image(csi.permap[d]);
          if AssertionLevel()>2 then
            hom:=GroupHomomorphismByImages(gp,dgp,gens,genims);
          else
            hom:=GroupHomomorphismByImagesNC(gp,dgp,gens,genims);
          fi;
          if hom=fail then Error("should not happen");fi;

          if d=i then
            # pull generator images back in original group
            genimgs[kn][i]:=List(genims,x->PreImagesRepresentativeNC(isom,x));

#if ForAny(Flat(genimgs),x->not x in Image(csi.permap[pools[poolnum][1]])) then
#  Error("imerrD");
#fi;

            # component is not moved we just need to conjugate with isom
            hom:=isom*hom*InverseGeneralMapping(isom);
          SetIsBijective(hom,true);
            if not IsInnerAutomorphism(hom) then
              Add(auts[poolnum],hom);
              AddSet(doesaut[kn],i);
            elif not IsOne(hom) then
              # not automorphism, but still acting
              AddSet(doesaut[kn],i);
            fi;
          elif d in pools[poolnum] then
            # we know already that the groups are isomorphic -- just store
            # the new automorphism

            # isomorphism is always to the first element in the pool
            if d=pools[poolnum][1] then
              # the group is the normal form -- we do not need to translate
              diso:=IdentityMapping(dgp);
            else
              # isomorphism to canonical form
              diso:=InverseGeneralMapping(isoms[d]);
            fi;
            genimgs[kn][i]:=List(genims,x->ImagesRepresentative(diso,x));

# if ForAny(Flat(genimgs),x->not x in Image(csi.permap[pools[poolnum][1]])) then
#   Error("imerrC");
# fi;
            hom:=isom*hom*diso;
          SetIsBijective(hom,true);
            if not IsInnerAutomorphism(hom) then
              Add(auts[poolnum],hom);
              AddSet(doesaut[kn],i);
            fi;
          else
            # isomorphism wasn't known yet -- we need to join pools
            # hom is the joiner, isom*hom*isoms[d]^-1 the conjugator of the
            # canonical map

            # the pool we add. We always join the image pool to the current
            # one.
            a:=First([1..Length(pools)],x->d in pools[x]);
            if a=fail then
              # we did not yet know layer d -- add it
              Info(InfoFFMat,2,"Found Layer ",d);
              Add(pools[poolnum],d);
              isoms[d]:=isom*hom;
              Assert(3,IsBijective(isoms[d]));
              # newly included component -- the generators *are* simply the
              # images
              genimgs[kn][i]:=List(genims,
                x->PreImagesRepresentativeNC(isoms[d],x));

#if ForAny(Flat(genimgs),x->not x in Image(csi.permap[pools[poolnum][1]])) then
  #Error("imerrB");
#fi;
              goodgens[d]:=conjgens;
              Add(process,d);
            else
              Error("This cannot happen??");
              # we knew the layer and join pools
              Info(InfoFFMat,2,"Join ",pools[a]," to ",pools[poolnum]);
              conj:=isom*hom;
              if not d=pools[a][1] then
                # translate to normal form
                conj:=conj*InverseGeneralMapping(isoms[d]);
              fi;

              # join pools
              for x in pools[a] do
                Add(pools[poolnum],x);
                isoms[x]:=conj*isoms[x]; # reroot isomorphism
              od;

              # translate automorphisms
              for x in auts[a] do
                Add(auts[poolnum],conj*x*InverseGeneralMapping(conj));
              od;

              # delete old
              pools[a]:=[];
              auts[a]:=[];
            fi;
          fi;
        od;
      od;
    fi;
    ii:=ii+1;
  od;

  wrs:=[];
  wrsoc:=[];
  dfgens:=[];
  dfimgs:=[];

#if ForAny(Flat(genimgs),x->not x in Image(csi.permap[pools[poolnum][1]])) then
#  Error("imerrA");
#fi;

  if Length(pools)=0 then
    # solvable
    d:=Group(());
    a:=GeneratorsOfGroup(csi.group);
    b:=List(a,x->One(d));
    if RUN_IN_GGMBI=false then
    RUN_IN_GGMBI:=true; # hack to skip Nice treatment
      hom:=GroupHomomorphismByImagesNC(csi.group,d,a,b:);
    RUN_IN_GGMBI:=false;
    else
      hom:=GroupHomomorphismByImagesNC(csi.group,d,a,b:Run_In_GGMBI:= true);
    fi;

    dci:=rec(isTrivial:=true);
    SetRecogDecompinfoHomomorphism(hom,dci);
    SetImagesSource(hom,Group(()));
    return hom;
  fi;

  dci.pools:=pools;
  dci.wreathemb:=[];
  dci.goodgens:=goodgens;
  dci.isoms:=isoms;
  dci.poollocalgens:=[];

  # now build each group
  for i in [1..Length(pools)] do
    pool:=pools[i];
    m:=Length(pools[i]);
    a:=Image(csi.permap[pools[i][1]]);
    a:=[a,Concatenation(auts[i],
            List(GeneratorsOfGroup(a),x->InnerAutomorphism(a,x)))];
    a:=AutomorphismRepresentingGroup(a[1],a[2]);
    aels[i]:=a;
    localgens:=List(poolimggens[i],x->Image(a[2],x));
    dci.poollocalgens[i]:=localgens;
    b:=SymmetricGroup(m);
    w:=WreathProduct(a[1],b);
    e:=List([1..m+1],x->Embedding(w,x));
    dci.wreathemb[i]:=e;
    biggens:=[];
    wrimages:=[];

    for j in Filtered([1..n],x->IsBound(csi.minranges[x])) do
      c:=csi.minranges[j];
    Info(InfoFFMat,2,"Images for ",j," ",c);

      # does any of the elements actually do something -- if so, what?

#      if true or ForAny(c,x->IsBound(doesaut[x])
#        and ForAny(pools[i],y->y in doesaut[x])) then
        Info(InfoFFMat,2,c," acts ",pools[i],j);
        for kn in c do
          #Add(biggens,csinice[kn]);
          Add(biggens,CSINiceGens(csi,kn));
          img:=One(w);
          for l in [1..Length(pools[i])] do
            if IsBound(genimgs[kn][pools[i][l]]) then
              d:=Image(e[l],CSIAelement(a,localgens,genimgs[kn][pools[i][l]]));
            else
              if pools[i][l]=j then
                d:=List(goodgens[pools[i][l]],
                        x->ImagesRepresentative(csi.permap[pools[i][l]],
                          CSIImageHomNr(csi,pools[i][l],x^CSINiceGens(csi,kn))));
                if isoms[j]<>false then
                  d:=List(d,x->PreImagesRepresentativeNC(isoms[j],x));
                fi;
                d:=Image(e[l],CSIAelement(a,localgens,d));
              else
                # component is not acted on as it lies higher
                d:=One(a[1]);
              fi;
            fi;
            img:=img*d;
          od;
          # is it a nontrivial permutation
          if IsBound(perms[kn]) and ForAny([1..Length(perms[kn])],
            x->IsBound(perms[kn][x]) and perms[kn][x]<>x) then
            # fill up fixed points
            for d in pool do
              if not IsBound(perms[kn][d]) then perms[kn][d]:=d;fi;
            od;

            # permuting part
            d:=PermList(List(perms[kn]{pool},x->Position(pool,x)));
            img:=img*Image(e[Length(pool)+1],d);
          fi;

          Add(wrimages,img);
        od;

    od;
    Add(wrs,w);
    a:=List([1..m],x->PerfectResiduum(Image(Embedding(w,x))));
    b:=TrivialSubgroup(w);
    for l in a do
      b:=ClosureGroup(b,l);
    od;
    SetDirectFactorsFittingFreeSocle(w,a);
    Add(wrsoc,b); # socle
    Add(dfgens,biggens);
    Add(dfimgs,wrimages);
  od;

  d:=DirectProduct(wrs);
  dci.dirprod:=d;
  dci.embeddings:=List([1..Length(pools)],x->Embedding(d,x));
  dci.aels:=aels;
  a:=List([1..Length(pools)],x->List(DirectFactorsFittingFreeSocle(wrs[x]),
                                y->Image(dci.embeddings[x],y)));
  dfs:=Concatenation(a);
  SetDirectFactorsFittingFreeSocle(d,dfs);
  a:=dfgens[1];
  b:=List(a,x->One(d));
  w:=TrivialSubgroup(d);
  for i in [1..Length(pools)] do
    e:=dci.embeddings[i];
    w:=ClosureGroup(w,Image(e,GeneratorsOfGroup(wrsoc[i])));
    for j in [1..Length(a)] do
      Assert(1,dfgens[i][j]=a[j]);

      b[j]:=b[j]*Image(e,dfimgs[i][j]);
    od;
  od;

  if IsPermGroup(csi.group) and AssertionLevel()>2 then
    hom:=GroupHomomorphismByImages(csi.group,d,a,b);
  else
    if RUN_IN_GGMBI=false then
    RUN_IN_GGMBI:=true; # hack to skip Nice treatment
    hom:=GroupHomomorphismByImagesNC(csi.group,d,a,b);
    RUN_IN_GGMBI:=false;
    else
      hom:=GroupHomomorphismByImagesNC(csi.group,d,a,b:Run_In_GGMBI:= true);
    fi;
  fi;
  if hom=fail then
    Error("fail!");
  fi;
  b:=Subgroup(d,b);
  #SetSocle(b,w);
  dci.socle:=w;
  dfs:=Filtered(dfs,x->IsSubset(b,x));
  SetDirectFactorsFittingFreeSocle(w,dfs);
  SetDirectFactorsFittingFreeSocle(b,dfs);
  SetImagesSource(hom,b);
  SetRecogDecompinfoHomomorphism(hom,dci);
  return hom;

end;


CSIElementAct:=function(dci,elm)
local
  # -- setup --
  csi,     # the recognition tree info, dci.csi
  pools,   # the pools, dci.pools
  result,  # the accumulated image of `elm` in the direct product (return value)
  # -- per-pool work --
  e,       # the wreath embeddings for the current pool, dci.wreathemb[i]
  a,       # the automorphism-representing-group data for the pool, dci.aels[i]
  img,     # the wreath-product image accumulated for the current pool
  perm,    # the permutation of the pool's components induced by `elm`
  b,       # the conjugated good generators (goodgens ^ elm) at the current component
  dp,      # depth of the conjugated element (which layer it lands in)
  d,       # WARNING: multi-role -- a list of generator images, then a wreath-image element, then a permutation
  # -- loop counters --
  i,       # loop index over pools
  l;       # loop index over the components of a pool
  csi:=dci.csi;
  pools:=dci.pools;
  result:=One(dci.dirprod);
  for i in [1..Length(pools)] do
    e:=dci.wreathemb[i];
    a:=dci.aels[i];
    img:=One(Image(e[1]));
    perm:=[];
    for l in [1..Length(pools[i])] do
      b:=List(dci.goodgens[pools[i][l]],x->x^elm);
      dp:=CSIDepthElm(csi,b[1]);
      Assert(2,ForAll(b,x->CSIDepthElm(csi,x)=dp));
      perm[l]:=Position(pools[i],dp);

      d:=List(dci.goodgens[pools[i][l]],
              x->ImagesRepresentative(csi.permap[dp],
                CSIImageHomNr(csi,dp,x^elm)));
      if dci.isoms[dp]<>false then
        d:=List(d,x->PreImagesRepresentativeNC(dci.isoms[dp],x));
      fi;
      d:=Image(e[l],CSIAelement(a,dci.poollocalgens[i],d));
      img:=img*d;
    od;

    # permuting part
    d:=PermList(perm);
    img:=img*Image(e[Length(pools[i])+1],d);
    result:=result*Image(dci.embeddings[i],img);

  od;
  return result;
end;

#############################################################################
##
#M  ImagesRepresentative(<hom>,<x>)
##
InstallMethod(ImagesRepresentative,"for recognition mappings",
        FamSourceEqFamElm,
  [ IsGroupGeneralMapping and
  HasRecogDecompinfoHomomorphism,IsMultiplicativeElementWithInverse], 0,
function(hom, x)
local
  d;  # the decomposition-info record of the homomorphism, RecogDecompinfoHomomorphism(hom)
  d:=RecogDecompinfoHomomorphism(hom);
  if IsBound(d.isTrivial) then return ();fi;
  if IsBound(d.fct) then
    return d.fct(x);
  else
    return CSIElementAct(RecogDecompinfoHomomorphism(hom),x);
  fi;
end);

#############################################################################
##
#M  PreImagesSetNC(<hom>,<x>)
##
InstallMethod(PreImagesSetNC,"for recognition mappings", CollFamRangeEqFamElms,
  [ IsGroupGeneralMapping and
  HasRecogDecompinfoHomomorphism,IsGroup], 0,
function(hom, U)
local
  gens,  # a small generating set of the subgroup U
  pre;   # preimages under `hom` of those generators
  gens:=SmallGeneratingSet(U);
  pre:=List(gens,x->PreImagesRepresentativeNC(hom,x));
  U:=RecogDecompinfoHomomorphism(hom).LiftSetup;
  U:=SubgroupByFittingFreeData(Source(hom),pre,gens,U.pcgs);
  return U;
end);


#TODO: Detect Nonsolvable permuters (via perms) and leave out of pool

BasePointsActionsOrbitLengthsStabilizerChain:=function(c)
local
  l,  # accumulated list of [base point, action, orbit length] triples (return value)
  o;  # the orbit record at the current level of the stabilizer chain
  l:=[];
  while c<>false do
    o:=c!.orb;
    Add(l,[o!.orbit[1],o!.op,Length(o!.orbit)]);
    c:=c!.stab;
  od;
  return l;
end;


FactorspaceActfun:=function(field,bas)
local
  heads;  # the heads info of the semi-echelonized basis `bas`, captured for use by the returned action function
  bas:=TriangulizedMat(bas);
  heads:=HeadsInfoOfSemiEchelonizedMat(bas,Length(bas[1]));
  return function(v,g)
    local
      c;  # position of the first nonzero entry of the sifted vector, used to normalize (scale) it
    v:=v*g;
    v:=SiftedVectorForGaussianRowSpace(field,bas, heads, v );
    c:=PositionNonZero(v);
    if c<=Length(v) then
      v:=Inverse(v[c])*v;
    fi;
    MakeImmutable(v);
    return v;
  end;
end;

INVTRANSPCACHE:=[];
AsInverseTranspose:=function(x,g)
local
  a,  # the transposed inverse of `g`: looked up in the cache INVTRANSPCACHE, or computed and then cached
  p;  # WARNING: two roles -- (1) the loop index scanning the cache, and
      #          (2) a reordering list built to move a frequently-hit entry towards the front of the cache
  a:=fail;
  p:=0;
  while a=fail and p<Length(INVTRANSPCACHE) do
    p:=p+1;
    if IsIdenticalObj(INVTRANSPCACHE[p][1],g) then
      a:=INVTRANSPCACHE[p][2];
      if p>50 then
        p:=Concatenation([p],[1..p-1],[p+1..Length(INVTRANSPCACHE)]);
        INVTRANSPCACHE:=INVTRANSPCACHE{p};
      fi;
    fi;
  od;
  if a=fail then
    a:=TransposedMat(Inverse(g));
    if Length(INVTRANSPCACHE)>100 then
      # clean out old
      INVTRANSPCACHE:=Concatenation([[g,a]],INVTRANSPCACHE{[1..90]});
    else
      INVTRANSPCACHE:=Concatenation([[g,a]],INVTRANSPCACHE);
    fi;
  fi;
  return OnRight(x,a);
end;


ModuleStructureBase:=function(mats)
local
  # -- local helper functions --
  orbtranslimit,     # helper: build an orbit with a transversal, aborting if it exceeds a length limit
  trymultipleorbits, # helper: try several seed vectors and keep the best (longest) orbit found
  siftchain,         # helper: sift an element down through the accumulated orbit/transversal chain
  # -- group / module setup --
  whole,             # the whole group Group(mats)
  f,                 # the field of the matrix list, FieldOfMatrixList(mats)
  mo,                # the GModule of the generators
  limit,             # orbit-length limit (1000) beyond which an orbit is abandoned
  gens,              # the current working generating set (shrinks and grows as the chain is (re)built)
  # -- composition-series basis --
  bas,               # basis of the space, adapted to the composition series
  basn,              # normalized (OnLines) images of the basis vectors under a generator
  dims,              # list of the dimensions of the composition-series terms
  # -- chain construction --
  use,               # the list of orbit records (one per base point) forming the constructed chain
  total,             # running product of the orbit lengths (the group order computed so far)
  fct,               # the action function currently in use (OnRight / OnLines / factor-space / dual action)
  vec,               # the seed vector(s) from which the current orbit is grown
  dict,              # dictionary for orbit-membership lookup of the current orbit
  new,               # the new stabilizer generators found by random sifting
  alt,               # counter of sifting attempts made
  sub,               # dimension of the subspace defining a factor-space action
  # -- heavily reused scratch variables --
  a,                 # WARNING: multi-role -- a basis-vector accumulator, a reversed list of dimensions,
                     #          a random group element, and a sift result (reassigned throughout)
  p,                 # WARNING: two roles -- a position/dimension index (integer) and a list of positions
  orb,               # WARNING: two roles -- an orbit (or the [dict,orb,t] triple returned by the helpers),
                     #          and a plain loop counter in the dual-action retry
  t,                 # WARNING: two roles -- a set of candidate vectors, and the transversal of an orbit
  b;                 # a saved base position index

  orbtranslimit:=function(vec,fct,limit)
  local
    dict,   # dictionary mapping each orbit point to its index in `orb`
    orb,    # the orbit (list of points) being enumerated (return component)
    t,      # transversal: for each orbit point, a group element carrying `vec` to it (return component)
    pnt,    # index of the orbit point currently being expanded
    g,      # loop variable over the generators
    img,    # image fct(orb[pnt],g) of the current point under a generator
    p,      # looked-up index of `img` in `dict` (or fail if new)
    coinc;  # tracks whether the orbit map coincided (a point was reached twice); assigned but not returned
    dict:=NewDictionary(vec,true,f^Length(vec));
    orb:=[vec];
    AddDictionary(dict,vec,1);
    t:=[One(gens[1])];
    pnt:=1;
    coinc:=true;
    while pnt<=Length(orb) do
      for g in gens do
        img:=fct(orb[pnt],g);
        p:=LookupDictionary(dict,img);
        if p=fail then
          Add(orb,img);
          if Length(orb)>limit then
            return fail;
          fi;
          AddDictionary(dict,img,Length(orb));
          Add(t,t[pnt]*g);
        elif coinc then
          coinc:=false;
        fi;
      od;
      pnt:=pnt+1;
    od;
    return [dict,orb,t];
  end;

  trymultipleorbits:=function(seeds,fct,limit)
  local
    best,   # the best orbit triple [dict,orb,t] found so far (return value)
    bestl,  # length of the best orbit found so far
    bestc,  # how many times the current best length has recurred (a stability counter)
    i,      # loop variable over the seed vectors
    orb;    # the orbit computed for the current seed
    best:=fail;
    bestl:=0;
    bestc:=0;
    for i in seeds do
      # form orbit with possible initial canonization
      orb:=orbtranslimit(fct(i,One(gens[1])),fct,limit);
      if orb<>fail and Length(orb[2])>1 then
        if Length(orb[2])>bestl then
          best:=orb;
          bestl:=Length(orb[2]);
          bestc:=1;
          if bestl>100 then
            return best;
          fi;
        else
          bestc:=bestc+1;
          # if we got 30 times the same length then use it.
          if bestc>30 then
            return best;
          fi;
        fi;
      elif bestl>0 then
        # we found at least one and others failed
        return best;
      fi;
    od;
    return best;
  end;

  siftchain:=function(use,x)
  local
    i,    # index of the chain level reached (returned as part of the failure marker)
    img,  # image of the current level's action / its looked-up orbit position
    u;    # loop variable over the chain records in `use`
    i:=1;
    for u in use do
      img:=u.fct(u.vec,x);
      img:=LookupDictionary(u.dict,img);
      if img=fail then
        return [fail,i,x];
      fi;
      x:=x/u.t[img];
      i:=i+1;
    od;
    return x;
  end;

  limit:=1000;

  whole:=Group(mats);
  f:=FieldOfMatrixList(mats);
  total:=1;
  gens:=mats;
  mo:=GModuleByMats(gens,f);
  bas:=MTX.BasesCompositionSeries(mo);
  dims:=List(bas,Length);
  if ForAll([2..Length(bas)],x->IsSubset(bas[x],bas[x-1])) then
    basn:=List([2..Length(bas)],x->Filtered(bas[x],y->not y in bas[x-1]));
    bas:=Concatenation(basn);
  else
    a:=ShallowCopy(bas[2]); # 1 is empty
    p:=3;
    while p<=Length(bas) do
      t:=Difference(bas[p],a);
      t:=Difference(t,bas[p-1]);
      repeat
        vec:=List([1..Length(bas[p])-Length(bas[p-1])],x->Random(t));
      until RankMat(Concatenation(a,vec))=Length(bas[p]);
      a:=Concatenation(a,vec);
      p:=p+1;
    od;
    bas:=a;
  fi;
  use:=[];

  while Length(gens)>0 do
    basn:=List(bas,x->OnLines(x,One(gens[1])));
    p:=PositionProperty(basn,x->ForAny(gens,y->OnLines(x,y)<>x));
    if p=fail then
      p:=PositionProperty(bas,x->ForAny(gens,y->x*y<>x));
      vec:=bas{[p..Minimum(p+10,Length(bas))]};
      fct:=OnRight;
    elif Size(f)^p<limit then
      p:=PositionProperty(dims,x->Size(f)^x>limit);
      if p=fail then
        p:=Length(dims);
      else
        p:=p-1;
      fi;
      orb:=OrbitsDomain(Group(gens),NormedRowVectors(VectorSpace(f,bas{[1..dims[p]]},Zero(bas[1]))),OnLines);
      orb:=Filtered(orb,x->Length(x)>1);
      Sort(orb,function(a,b) return Length(a)>Length(b);end);
      vec:=List(orb,x->x[1]);
      fct:=OnLines;
    else
      vec:=bas{[p..Minimum(p+10,Length(bas))]};
      fct:=OnLines;
    fi;

    # nontrivial orbit
    orb:=trymultipleorbits(vec,fct,limit);
    if orb=fail then
      b:=p;

      a:=Reversed(Filtered(dims,x->x<b));
      p:=fail;
      while p=fail and a[1]<>0 do
        sub:=a[1];
        fct:=FactorspaceActfun(f,bas{[1..sub]});
        basn:=List(bas,x->fct(x,One(gens[1])));
        p:=Filtered([a[1]+1..Minimum(a[1]+10,Length(bas))],
                            x->ForAny(gens,y->fct(basn[x],y)<>basn[x]));
        if p=[] then
          p:=Filtered([a[1]+10..Length(bas)],
                              x->ForAny(gens,y->fct(basn[x],y)<>basn[x]));
        fi;
        if p=[] then
          p:=fail;
        else
          # try to find one with orbit not just length p but also not too
          # large
          p:=Reversed(p);
          while Length(p)>1 and Size(f)^(p[1]-a[1])/(Size(f)-1)>limit do
            p:=p{[2..Length(p)]};
          od;
          p:=p[1];
        fi;
        a:=a{[2..Length(a)]};
      od;
      if p<>fail then
        orb:=trymultipleorbits(bas{[p..Minimum(p+10,Length(bas))]},fct,limit);
        if orb=fail then
          # try dual action
          fct:=AsInverseTranspose;
          p:=First([b..Length(bas)],x->ForAny(gens,y->fct(bas[x],y)<>bas[x]));
          if p<>fail then
            p:=[p];
            for orb in [1..5] do
              Add(p,p[1]+orb);
              Add(p,p[1]-orb);
            od;
            p:=Filtered(p,x->x>0 and x<Length(bas));
            orb:=trymultipleorbits(bas{p},fct,limit);
          fi;
          if orb=fail then
            Info(InfoFFMat,2,"even too long in dual");
            return rec(points:=[],ops:=[]);
          else
            Info(InfoFFMat,2,"dual ");
          fi;
        else
          Info(InfoFFMat,2,"factor space ");
        fi;
      else
        return rec(points:=[],ops:=[]);
        Error("p=fail -- should not happen");
      fi;

    fi;

    dict:=orb[1];
    t:=orb[3];
    orb:=orb[2];
    vec:=orb[1];
    Info(InfoFFMat,2,"got Length ",Length(orb),":",
      Collected(Factors(total*Length(orb))));

    Add(use,rec(vec:=vec,fct:=fct,orb:=orb,dict:=dict,t:=t,gens:=gens));

    #sift 20 random elements
    new:=[];
    p:=Maximum(Minimum(1000,Length(orb)*Length(gens)),10);
    alt:=0;
    a:=fail;
    while a=fail and Length(new)<20 and alt<p do
      a:=PseudoRandom(whole);
      a:=siftchain(use,a);
      alt:=alt+1;
      if a[1]<>fail then
        if fct(vec,a)<>vec then
          Error("not stab");
        fi;
        if not IsOne(a) and not a in new then
          Add(new,a);
          alt:=0;
        fi;
        a:=fail; # sifted through
      fi;
    od;

    if a<>fail then
      gens:=ShallowCopy(use[a[2]].gens);
      Add(gens,a[3]);
      use:=use{[1..a[2]-1]};
      total:=Product(List(use,x->Length(x.orb)));
      Info(InfoFFMat,2,"construction fail on level ",a[2]," of ",Length(use),
        " -- redo ",total);
    else
      total:=total*Length(orb);
      gens:=new;
    fi;

  od;
  return rec(points:=List(use,x->x.vec),ops:=List(use,x->x.fct),
             order:=total);
end;

MATGRP_AddGeneratorToStabilizerChain:="2bdefined";

MATGRP_StabilizerChainInner:=
  function( gens, size, layer, cand, opt, parentS )
    # Computes a stabilizer chain for the group generated by gens
    # with known size size (can be false if not known). This will be
    # layer layer in the final stabilizer chain. cand is a (shared)
    # record for base point candidates and opt the (shared) option
    # record. This is called in StabilizerChain and calls itself.
    # It also can be called if a new layer is needed.
    local
      # -- chain being built --
      S,          # the stabilizer-chain record for this layer (built and returned)
      base,       # the base of S, S!.base
      next,       # record holding the next base point and its action (point/op)
      # -- Schreier-tree shallowing --
      merk,       # saved number of orbit generators before making the Schreier tree shallower
      merk2,      # saved number of strong generators before making the Schreier tree shallower
      # -- random stabilizer generation (dead code after the early `return S`) --
      stabgens,   # random stabilizer generators collected for the next layer
      x,          # a random element generated from above
      i;          # loop counter (over random-generator attempts / verification elements)

    Info(InfoGenSS,4,"Entering MATGRP_StabilizerChainInner layer=",layer);
    next:=rec(point:=cand.points[1],op:=cand.ops[1]);
    #next := GENSS_NextBasePoint(gens,cand,opt,parentS);
    #cand := next.cand;   # This could have changed

    S := GENSS_CreateStabChainRecord(parentS,gens,size,
                                     next.point,next.op,cand,opt);
    base := S!.base;

    Info( InfoGenSS, 3, "Entering orbit enumeration layer ",layer,"..." );
    repeat
        Enumerate(S!.orb,opt.OrbitLengthLimit);
        if not(IsClosed(S!.orb)) then
            if opt.FailInsteadOfError then
                return "Orbit too long, increase opt.OrbitLengthLimit";
            else
                Error("Orbit too long, increase opt.OrbitLengthLimit");
            fi;
        fi;
    until IsClosed(S!.orb);
    Info(InfoGenSS, 2, "Layer ", layer, ": Orbit length is ", Length(S!.orb));

    if layer > 1 then
        parentS!.stab := S;   # such that from now on random element
                              # generation works!
    else
        if (Length(S!.orb) > 50 or S!.orb!.depth > 5) and
           S!.opt.OrbitsWithLog then
            Error("ARGH!");
            Info(InfoGenSS, 3, "Trying to make Schreier tree shallower (depth=",
                 S!.orb!.depth,")...");
            merk := Length(S!.orb!.gens);
            merk2 := Length(S!.stronggens);
            MakeSchreierTreeShallow(S!.orb);
            Append(S!.stronggens,S!.orb!.gens{[merk+1..Length(S!.orb!.gens)]});
            Append(S!.layergens,[merk2+1..Length(S!.stronggens)]);
            Info(InfoGenSS, 3, "Depth is now ",S!.orb!.depth);
        fi;
    fi;
    S!.orb!.gensi := List(S!.orb!.gens,x->x^-1);

    # as we add normalizing, we are done
    return S;

    # Are we done?
    if size <> false and Length(S!.orb) = size then
        S!.proof := true;
        Info(InfoGenSS,4,"Leaving MATGRP_StabilizerChainInner layer=",layer);
        return S;
    fi;

    # Now create a few random stabilizer elements:
    stabgens := EmptyPlist(opt.RandomStabGens);
    for i in [1..opt.RandomStabGens] do
        x := GENSS_RandomElementFromAbove(S,layer);
        if not(S!.IsOne(x)) then
            Add(stabgens,x);
        fi;
    od;
    Info(InfoGenSS,3,"Created ",opt.RandomStabGens,
         " random stab els, ",
         Length(stabgens)," non-trivial.");

    if Length(stabgens) > 0 then   # there is a non-trivial stabiliser
        Info(InfoGenSS,3,"Found ",Length(stabgens)," non-trivial ones.");
        if size <> false then
            S!.stab := MATGRP_StabilizerChainInner(stabgens,size/Length(S!.orb),
                                                   layer+1,cand,opt,S);
        else
            S!.stab := MATGRP_StabilizerChainInner(stabgens,false,
                                                   layer+1,cand,opt,S);
        fi;
        if IsString(S!.stab) then return S!.stab; fi;
        if opt.ImmediateVerificationElements > 0 then
            Info(InfoGenSS,2,"Doing immediate verification in layer ",
                 S!.layer," (",opt.ImmediateVerificationElements,
                 " elements)...");
            i := 0;
            while i < opt.ImmediateVerificationElements do
                i := i + 1;
                x := GENSS_RandomElementFromAbove(S,layer);
                if MATGRP_AddGeneratorToStabilizerChain(S!.topS,x) then
                    Info( InfoGenSS, 2, "Immediate verification found error ",
                          "(layer ",S!.layer,")..." );
                    i := 0;
                fi;
            od;
        fi;

        S!.proof := S!.stab!.proof;   # hand up information
    else
        # We are not sure that the next stabiliser is trivial, but we believe!
        Info(InfoGenSS,3,"Found no non-trivial ones.");
        S!.proof := false;
    fi;

    Info(InfoGenSS,4,"Leaving MATGRP_StabilizerChainInner layer=",layer);
    return S;
  end;


MATGRP_AddGeneratorToStabilizerChain:=
  function( S, el,pt,op )
    # Increases the set represented by S by the generator el.
    local
      SS,           # the stabilizer-chain node at which the new generator (or new layer) is added
      r,            # the sift result record for the element `el`
      n,            # loop variable over the component names when copying a freshly built trivial-group chain
      newstrongnr;  # index that the newly added strong generator will occupy in the strong-generator list

    if IsBound(S!.trivialgroup) and S!.trivialgroup then
        if S!.IsOne(el) then
            return false;
        fi;
        #SS := StabilizerChain(Group(el),S!.opt);
        SS:=StabilizerChain(Group(el),rec(Cand:=rec(ops:=[op],points:=[pt]),Reduced:=false,StrictlyUseCandidates:=true));

        if IsString(SS) then return SS; fi;
        for n in NamesOfComponents(SS) do
            S!.(n) := SS!.(n);
        od;
        Unbind(S!.trivialgroup);
        return true;
    fi;

    r := SiftGroupElement( S, el );
    # if this is one then el is already contained in the stabilizer chain
    if r.isone then     # already in the group!
        return false;
    fi;
    # Now there remain two cases:
    #  (1) the sift stopped somewhere and we have to add a generator there
    #  (2) the sift ran all through the chain and the element still was not
    #      the identity, then we have to prolong the chain
    if r.S <> false then   # case (1)
        SS := r.S;
        Info( InfoGenSS, 2, "Adding new generator to stab. chain ",
              "in layer ", SS!.layer, " from ",Length(SS!.stronggens) );
        Add(SS!.stronggens,r.rem);
        Add(SS!.layergens,Length(SS!.stronggens));
        AddGeneratorsToOrbit(SS!.orb,[r.rem]);
        Add(SS!.orb!.gensi,r.rem^-1);
        newstrongnr := Length(SS!.stronggens);
        Info( InfoGenSS, 4, "Entering orbit enumeration layer ",SS!.layer,
              "..." );
        repeat
            Enumerate(SS!.orb,S!.opt.OrbitLengthLimit);
            if not(IsClosed(SS!.orb)) then
                if S!.opt.FailInsteadOfError then
                    return "Orbit too long, increase S!.opt.OrbitLengthLimit";
                else
                    Error("Orbit too long, increase S!.opt.OrbitLengthLimit!");
                fi;
            fi;
        until IsClosed(SS!.orb);
        Info( InfoGenSS, 4, "Done orbit enumeration layer ",SS!.layer,"..." );
        SS!.proof := false;
    else   # case (2)
        # Note that we do not create a pr instance here for one
        # generator, this will be done later on as needed...
        SS := r.preS;
        newstrongnr := Length(SS!.stronggens)+1;  # r.rem will end up there !
        SS!.stab := MATGRP_StabilizerChainInner([r.rem],false,
                           SS!.layer+1,rec(points:=[pt],ops:=[op]), SS!.opt, SS );
        if IsString(SS!.stab) then return SS!.stab; fi;
        SS := SS!.stab;
    fi;
    # Now we have added a new generator (or a new layer) at layer SS,
    # the new gen came from layer S (we were called here, after all),
    # thus we have to check, whether all the orbits between S (inclusively)
    # and SS (exclusively) are also closed under the new generator r.rem,
    # we add it to all these orbits, thereby also making the Schreier trees
    # shallower:
    while S!.layer < SS!.layer do
        Info(InfoGenSS,2,"Adding new generator to orbit in layer ",S!.layer);
        Add(S!.layergens,newstrongnr);
        AddGeneratorsToOrbit(S!.orb,[r.rem]);
        Add(S!.orb!.gensi,r.rem^-1);
        S := S!.stab;
    od;
    # Finally, we have to add it to the product replacer!
    AddGeneratorToProductReplacer(S!.opt!.pr,r.rem);
    return true;
  end;


SolvableBSGS:=function(arg)
local
  # -- local helper functions --
  CHAINTEST,             # helper: consistency checker for a stabilizer chain (debugging assertions)
  normalizingGenerator,  # helper: add a generator to the chain S and return the new strong generators
  solvNC,                # helper: compute one layer of the pcgs (a solvable normal-closure step)
  # -- group / field setup --
  gens,                  # the current working generating set (reassigned repeatedly: input gens, random gens,
                         #   then the reversed pcgs to re-sift, ...) -- always "the generators to process"
  acter,                 # the acting generators used to close each layer under conjugation
  sz,                    # the known target group order (arg[3]), or fail if not supplied
  df,                    # the default field of the matrix group, DefaultFieldOfMatrixGroup(a)
  shortlim,              # short-orbit limit (scaled by the field size)
  bound,                 # Dixon bound on the solvable length, used as a safety limit
  opt,                   # options record passed to StabilizerChain (later reset to 1 to free memory)
  goodbase,              # the reference base (points/ops) used when rebuilding chains
  CBase,                 # the reference stabilizer chain / its base, reused to rebuild chains on the same base
  # -- stabilizer chain working state --
  S,                     # the stabilizer chain currently being built up (the main working chain)
  oldsz,                 # saved Size(S) before an addition, to detect whether the group grew
  # -- pcgs accumulation --
  pcgs,                  # the pc generating sequence being accumulated (a list, later the Pcgs object)
  U,                     # the list of pcgs elements found for the current layer (returned by solvNC)
  prime,                 # the relative order (a prime) of the current layer
  relord,                # list of the relative orders of the pcgs
  depths,                # list of the layer-boundary indices (IndicesEANormalSteps) of the pcgs
  n,                     # the number of pc generators, Length(pcgs)
  x,                     # a generator/element currently being sifted or processed
  c,                     # derived-series depth counter (guards against a non-solvable input via `bound`)
  # -- strong generators and layer vectors --
  strongs,               # the strong generators, arranged so that they coincide with the pcgs
  stronglevs,            # the layer number of each strong generator
  stronglayers,          # the layer of each strong generator (layers indexed by stronglevs)
  slvec,                 # exponent vectors of the strong generators (for decomposition)
  layerzero,             # zero-vector template for each layer
  laynum,                # the number of layers, Length(depths)-1
  layers,                # for each pcgs index, the number of the layer it belongs to
  # -- decomposition indexing --
  baspts,                # for each pcgs element, the index of the base point (level) that moves it
  levp,                  # for each base point, the list of pcgs indices decomposed at that point
  blocksz,               # for each pcgs element, its block size (orbit-position divisor)
  # -- heavily reused scratch variables --
  w,                     # WARNING: multi-role -- a ModuleStructureBase record, then a base-point list, then
                         #   the current candidate pcgs generator (a group element inside solvNC)
  a,                     # WARNING: multi-role -- the group, an integer dimension, a chain-walking node, a
                         #   boolean flag, a matrix, and a saved copy of the pcgs (reassigned throughout)
  p,                     # WARNING: multi-role -- a prime, a list of orbit positions, and a list of strong-
                         #   generator positions (reassigned throughout)
  xp,                    # WARNING: two roles -- a loop counter over `gens`, and later the index list [1..|pcgs|]
  sel,                   # a selected list of indices (of orbit-moving generators / of a strong layer)
  # -- loop counters --
  i,                     # loop index
  j,                     # loop index
  # -- vestigial (assigned but effectively unused) --
  check,                 # assigned true/false but never read
  primes,                # assigned Set(Factors(Size(CBase))) but never read afterwards
  stabs,                 # assigned [] but only populated/used in commented-out code
  slpval,                # a nested function, defined but only called from commented-out code
  vals;                  # only referenced by slpval and in commented-out code

  goodbase:=[];
  CHAINTEST:=function(X,str)
    while X<>false do
      #if IsBound(X!.stronggens) and
      #  Length(X!.layergens)>Length(Factors(Size(X))) then
      #  Error("UGH");
      #fi;
      if not X!.orb!.orbit[1] in goodbase.points then
        Error("new point!");
      fi;
      #if Length(X!.orb!.gens)<>Length(Set(X!.orb!.gens)) then
      #  Error("duplicate!");
      #fi;
      #if IsBound(X!.opt) and IsBound(X!.opt.StrictlyUseCandidates) and
      #  X!.opt.StrictlyUseCandidates=false then Error("eh5!"); fi;

      if IsBound(X!.orb!.gens) and Length(X!.orb!.gens)<>Length(Factors(Size(X))) then
        Error("length!");
      fi;

      if IsBound(X!.orb!.gens) and Length(Difference(X!.orb!.gens,strongs))>0 then
        Error("XTRA!");
      fi;

      if IsBound(X!.orb!.gensi) and List(X!.orb!.gens,Inverse)<>X!.orb!.gensi then
        Error(str,"inverse!");
      fi;
      if IsBound(X!.orb!.gensi) and ForAny(X!.orb!.schreiergen,IsInt)  and
        Maximum(Filtered(X!.orb!.schreiergen,IsInt))>Length(X!.orb!.gensi)
        then
        Error("length!");
      fi;
      X:=X!.stab;
    od;
  end;

  gens:=arg[Minimum(2,Length(arg))];
  if IsGroup(gens) then
    a:=gens;
    gens:=GeneratorsOfGroup(gens);
  else
    a:=Group(gens);
  fi;
  if false and Length(arg)>2 and Length(gens)>5+Length(Factors(arg[3])) then
    gens:=List([1..5+Length(Factors(arg[3]))],x->PseudoRandom(a));
    gens:=Set(gens);
    a:=Group(gens);
    check:=true;
    #SetSize(a,sz);
  else
    check:=false;
  fi;
  if IsBound(arg[3]) then
    sz:=arg[3];
    #SetSize(a,sz);
  else
    sz:=fail;
  fi;
  SetIsSolvableGroup(a,true);

  if Length(arg)=1 then
    acter:=gens;
  else
    acter:=arg[1];
    if IsGroup(acter) then
      acter:=GeneratorsOfGroup(acter);
    fi;
  fi;

  shortlim:=Size(DefaultFieldOfMatrixGroup(a))*3000;

  df:=DefaultFieldOfMatrixGroup(a);
  if false then
    # try vectors from big group
    w:=BasisVectorsForMatrixAction(Group(acter));
    w:=ImmutableMatrix(df,w);
    if Size(df)=2 then
      ops:=List(w,x->OnRight);
    else
      ops:=List(w,x->OnLines);
    fi;
  fi;
  w:=ModuleStructureBase(gens);
  if IsBound(arg[3]) and IsBound(w.order) and w.order<arg[3] then
    for i in gens[1] do
      Add(w.points,i);
      Add(w.ops,OnRight);
    od;
  fi;

  opt:=rec(RandomStabGens:=10,
         Cand:=rec(points:=w.points,ops:=w.ops),
      #TODO XXX Find better strategy for limit iof field is larger
         VeryShortOrbLimit:=shortlim
                             );
  FUNCSPACEHASH:=[]; # needed in base point selection code
  CBase:=StabilizerChain(a,opt);
  if sz<>fail and Size(CBase)<sz then
    Info(InfoWarning,1,"Wrong size -- redo with `doall' option");
    return fail;
  fi;

  primes:=Set(Factors(Size(CBase)));
  FUNCSPACEHASH:=[];
  w:=BasePointsActionsOrbitLengthsStabilizerChain(CBase);
  Info(InfoGenSS,2,"Suggested Base Points",List(w,x->Position(opt.Cand.points,x[1])),
                   " Lengths ",List(w,x->x[3]));

  goodbase:=BaseStabilizerChain(CBase);
  w:=1;
  opt:=1;

  # Dixon's (1968 - the solvable length of a solvable linear group) bounds
  if IsPerm(gens[1]) then
    a:=Length(MovedPoints(gens));
    bound := Int( LogInt( Maximum(1,a ^ 5), 3 ) / 2 );
  elif IsMatrixOrMatrixObj(gens[1]) then
    a:=NrRows(gens[1]); # dimension

    # we would need proper log
    #bound := Int( (8.55*Log(a+1,10)+0.36 );
    # since it does not exist, replace 8.55 by 9 and simplify rounded up.
    # this anyhow only matters for catching wrong input, so its harmless
    bound := Log((a+1)^9,10)+1;
  else
    # no bound known
    bound:=infinity;
  fi;

  normalizingGenerator:=function(g)
  local
    oldsz,      # saved Size(S) before adding `g`, to detect that the group grew
    oldstrong,  # saved copy of the strong generators before adding `g`
    pba,        # the base points actually moved by `g`
    pbp;        # the chosen base point (the first one moved by `g`)

    oldsz:=Size(S);
    oldstrong:=ShallowCopy(StrongGenerators(S));
    # work around the ommission of prior redundant base points
    pba:=Filtered([1..Length(goodbase.points)],x->goodbase.ops[x](goodbase.points[x],g)<>goodbase.points[x]);
    pbp:=pba[1];

    #Print("\n");
    #View(S);
    #Print("\n use ",pba,"\n");

    MATGRP_AddGeneratorToStabilizerChain(S,g,goodbase.points[pbp],goodbase.ops[pbp]);

    while Size(S)=oldsz do
      Error("orderr BUG !!!\n");
      MATGRP_AddGeneratorToStabilizerChain(S,g);
    od;
    return Difference(StrongGenerators(S),oldstrong);
  end;

  solvNC:=function()
  local
    N,          # a fresh copy of the chain S, used for membership tests while building the layer
    pow,        # the relative order (power) of `w` being determined
    wp,         # successive powers of `w` used to determine `pow`
    process,    # work-list of elements (the conjugates) still to process
    layergens,  # the generators found so far for the current layer
    U,          # the accumulated layer elements (return value)  [shares outer name]
    g,          # loop variable over `process` (an element)
    sift,       # the sift result of `g` through S
    h,          # loop variable over `layergens`
    comm,       # a commutator Comm(g,h) (tests whether the layer really is central/abelian enough)
    V,          # the new strong generators returned by normalizingGenerator(g)
    x;          # loop variable over `acter`

    N:=StabilizerChain(Group(StrongGenerators(S)),rec(Base:=CBase,Size:=Size(S),Reduced:=false,StrictlyUseCandidates:=true)); # really should be copy

    # determine relative order
    pow:=2;
    wp:=w*w;
    while not SiftGroupElement(N,wp).isone do
      wp:=wp*w;
      pow:=pow+1;
    od;
    if not IsPrimeInt(pow) then
      prime:=Factors(pow)[1];
      w:=w^(pow/prime);
    else
      prime:=pow;
    fi;
    Info(InfoFFMat,2,"Relative order ",pow," prime ",prime);

    process:=[w];
    layergens:=[];
    U:=[];
    for g in process do
      sift:=SiftGroupElement(S,g);
      if not sift.isone then
        Info(InfoFFMat,2,"SS=",Size(S)," ",Length(U));

        for h in layergens do
          comm:=Comm(g,h);
  #Print(Order(comm),"\n");
          if not SiftGroupElement(N,comm).isone then
            # wrong element -- get new generator for layer below
            a:=false;
            w:=comm;
            S:=N; # don't we want to forget what we added?
            return false;
          fi;
        od;

        # the sifted element will be as good as generator, but allows for
        # cleaner decomposition.
        g:=sift.rem;

        V:=normalizingGenerator(g);
        #g:=V[1];
  if g in U then Error("duplicate");fi;
        U:=Concatenation([g],U);
        Add(layergens,g);
        for x in acter do
          Add(process,g^x);
        od;

      fi;
    od;
    a:=true;
    return U;

  end;

  S:=StabilizerChain(Group(One(gens[1])),rec(Base:=CBase,Reduced:=false,StrictlyUseCandidates:=true));

  pcgs:=[];
  depths:=[];
  relord:=[];
  xp:=1;
  while xp<=Length(gens) do
    Info(InfoFFMat,2,"Processing ",xp);
    x:=gens[xp];
    # feed through derived series if necessary
    while not SiftGroupElement(S,x).isone do
      c:=0;
      w:=x;
      a:=false;
      oldsz:=Size(S);
      while not a do
        c:=c+1;
         if c>bound then
           # not solvable
           Error("not solvable");
         fi;

        U:=solvNC(); # will change a,w

      od;
Info(InfoFFMat,2,"Found Layer ",Size(S)/oldsz);
#View(S);
Info(InfoFFMat,2,"n");
      if prime^Length(U)<>Size(S)/oldsz then Error("layerlen"); fi;
      pcgs:=Concatenation(U,pcgs);
      Add(depths,Length(pcgs));
      Append(relord,ListWithIdenticalEntries(Length(U),prime));

    od;
    xp:=xp+1;
  od;
  n:=Length(pcgs);
  depths:=Reversed(List(depths,x->n+1-x));
  relord:=Reversed(relord);
  Add(depths,n+1);

  w:=BasePointsActionsOrbitLengthsStabilizerChain(S);
  Info(InfoGenSS,2,"USED Base Points",List(w,x->Position(goodbase.points,x[1])),
                   " Lengths ",List(w,x->x[3]));


  # now build the chain once more with memory so we can decompose

  blocksz:=[];

  Info(InfoFFMat,2,"NOW",Size(S)," ",Length(pcgs));
  if Length(Factors(Size(S)))<>Length(pcgs) then Error("redundancies");fi;

  # now sift the generators through so that the strong generators *ARE* a pcgs
  gens:=Reversed(pcgs); # start with low level gens again
  pcgs:=[];
  acter:=[];
  #S:=StabilizerChain(Group(One(gens[1])),rec(Base:=CBase,Reduced:=false,StrictlyUseCandidates:=true));
  #S := GENSS_CreateStabChainRecord(false,
#                                   [],1,
#                                   goodbase.points[1],
#                                   goodbase.ops[1],goodbase,
#                                   rec(Base:=CBase,Reduced:=false,StrictlyUseCandidates:=true,InitialHashSize:=10));
  S:=false; # work around special treatment for trivial group.
  oldsz:=1;
  stabs:=[];
  xp:=1;
  strongs:=[];
  stronglevs:=[];
  while xp<=Length(gens) do


    if S<>false then
      #CHAINTEST(S,"E");
      oldsz:=Size(S);
      Info(InfoFFMat,2,"ProcessiNg ",xp," ",Size(S));


#if Length(StrongGenerators(S))>Length(Factors(Size(S))) then Error("more2\n");fi;
      x:=gens[xp];
      x:=SiftGroupElement(S,x).rem;
      #SiftGroupElement(S,x);

      if not IsOne(x) then
        normalizingGenerator(x);
      fi;
    else
      x:=gens[1];
      # OrbitsWithLog will add powers as strong generators to make the tree shallower.
      # This messes with the correspondence of strong generators and pcgs elements

      S:=StabilizerChain(Group(x),rec(Base:=CBase,Reduced:=false,StrictlyUseCandidates:=true,
        OrbitsWithLog := false,RandomStabGens:=1,Size:=Order(x)));

      # Remove those ~!@#$ extra generators (powers of x) that are added
      # to make the tree shallower, but cause problems here.
      strongs:=[x];

      a:=S;

      while a<>false do
        a!.stronggens:=ShallowCopy(strongs);
        a!.layergens:=[1];
        if Length(a!.orb!.gens)>0 then
          if a!.orb!.gens<>strongs then
            a!.orb:=Orb(ShallowCopy(strongs),a!.orb[1],a!.orb!.op,rec(schreier:=true));
            a!.orb!.gensi:=List(a!.orb!.gens,Inverse);
            repeat
              Enumerate(a!.orb);
            until IsClosed(a!.orb);
            a!.orb!.depth:=Size(a!.orb)-1;
            #if Length(a!.orb)=1 then
            #  a!.orb!.schreiergen:=ShallowCopy(strongs);
            #else
            #  a!.orb!.schreiergen:=[];
            #fi;
          fi;
        fi;
        a:=a!.stab;
      od;

      #CHAINTEST(S,"F2");

      strongs:=[];

    fi;

    if Size(S)=oldsz then Error("no change!");fi;

    Add(pcgs,x);
    a:=Set(Filtered(StrongGenerators(S),x->not x in strongs));
    a:=Difference(a,List([0..Order(x)],y->x^y));
    if Length(a)>0 then
      # redo
      Info(InfoFFMat,1,"nanu-redo");
      #Error("nanu");
      return SolvableBSGS(arg[1],arg[2],arg[3]);
    fi;
    Append(strongs,[x]);

    #CHAINTEST(S,"F");


    Append(stronglevs,ListWithIdenticalEntries(1,n+1-xp));
    #if n+1-xp in depths then
    #  Add(stabs,StructuralCopy(S));
    #fi;
    xp:=xp+1;
  od;

  pcgs:=Reversed(pcgs);

  CBase:=BaseStabilizerChain(S);
  w:=BasePointsActionsOrbitLengthsStabilizerChain(S);
  Info(InfoGenSS,2,"Used Base Points",List(w,x->Position(goodbase.points,x[1])),
                   " Lengths ",List(w,x->x[3]));

  baspts:=[];
  levp:=[];
  xp:=[1..Length(pcgs)];
  a:=S;
  for i in [1..Length(CBase.points)] do
    p:=List(xp,x->Position(a!.orb!.orbit,CBase.ops[i](CBase.points[i],pcgs[x])));
    sel:=Filtered([1..Length(xp)],x->p[x]<>1);
    baspts{xp{sel}}:=ListWithIdenticalEntries(Length(sel),i);
    levp[i]:=xp{sel};
    blocksz{xp{sel}}:=p{sel}-1;
    xp:=Difference(xp,xp{sel});
    a:=a!.stab;
  od;

  slpval:=function(w)
  local
    a,  # the value accumulated while evaluating the straight-line-program line
    i;  # loop index over the operands of the line
    a:=w[2]*vals[w[1]];
    for i in [3,5..Length(w)-1] do
      a:=(a+w[i+1]*vals[w[i]]) mod p;
    od;
    return a;
  end;

  layerzero:=[];
  laynum:=Length(depths)-1;
  layers:=List([1..n],x->First([laynum,laynum-1..1],y->x>=depths[y]));
  stronglayers:=layers{stronglevs};
  # get vectors for strong generators (b/c we decompose into these)
  slvec:=[];
  for i in [1..laynum] do
    p:=relord[depths[i]];
    a:=IdentityMat(depths[i+1]-depths[i]);
    layerzero[i]:=Zero(a[1]);
    sel:=Positions(stronglayers,i);
    slvec{sel}:=a{List(strongs{sel},x->Position(pcgs,x))-depths[i]+1};

#    layervecs:=ListWithIdenticalEntries(n,fail); # this will trap use of
#    # higher gens
#    layervecs{[depths[i]..depths[i+1]-1]}:=a;
#
#    for j in [depths[i+1]..n] do
#      layervecs[j]:=Zero(a[1]);
#    od;
#    slp:=LinesOfStraightLineProgram(SLPOfElms(strongs{sel}));
#    vals:=Reversed(layervecs);
#    for j in slp do
#      if ForAll(j,IsList) then
#       #last line
#       slvec{sel}:=List(j,x->slpval(x));
#      elif IsList(j[1]) then
#       Error("reassign syntax?");
#      else
#       Add(vals,slpval(j));
#      fi;
#    od;
  od;
  ForgetMemory(strongs);
  ForgetMemory(gens);
  ForgetMemory(S);

  S!.pcgs:=pcgs;
  S!.layerzero:=layerzero;
  S!.layerprimes:=relord{depths{[1..laynum]}};
  S!.strongvecs:=slvec;
  S!.layranges:=List([1..laynum],x->[depths[x]..depths[x+1]-1]);
  S!.revranges:=List([1..laynum],x->Reversed([1..Length(S!.layranges[x])]));

  # reindexing
  a:=S;
  while a<>false do
    p:=List(a!.orb!.gensi,x->Position(strongs,x^-1));
    if fail in p then
      Error("not found!");
    fi;
    a!.gensistrongpos:=p;
    a!.gensilayers:=stronglayers{p};
    a:=a!.stab;
  od;

  #TODO: Use chain for decomposition -- at the moment it uses subgroups
  a:=pcgs;
  pcgs:=PcgsByPcSequenceCons(IsPcgsDefaultRep,
    IsPcgsMatGroupByStabChainRep and IsPcgs and IsPrimeOrdersPcgs,
    FamilyObj(gens[1]),pcgs,[]);
  pcgs!.stabilizerChain:=S;
  #pcgs!.stabs:=Reversed(stabs);
  pcgs!.laynums:=[1..laynum];
  pcgs!.layranges:=S!.layranges;
  pcgs!.baspts:=baspts;
  pcgs!.blocksz:=blocksz;
  pcgs!.levp:=levp;
  pcgs!.invpows:=
    List([1..Length(pcgs)],x->List([1..relord[x]-1],y->a[x]^(-y)));
  SetRelativeOrders(pcgs,relord);
  SetIndicesEANormalSteps(pcgs,depths);

  return rec(
    pcgs:=pcgs,
    depths:=depths,
    relord:=relord,
    baspts:=baspts,
    blocksz:=blocksz);
end;

# old ersion of exponent routines, allowing for arbitrary orbit arrangemnts
MatPcgsSift:=function( S, x,l )
local
  # -- setup --
  vecs,   # the strong-generator vectors, S!.strongvecs
  prime,  # the relative order (prime) of the requested layer `l`
  v,      # the accumulated exponent vector for layer `l` (return value)
  # -- per-level sifting --
  o,      # the orbit at the current level of the chain
  p,      # the image point o!.op(o[1],x)
  po,     # position of `p` in the orbit `o`
  r,      # the Schreier-generator index while sifting through the Schreier tree
  # -- declared but unused --
  preS;   # assigned false but never read
  v:=S!.layerzero[l];
  prime:=S!.layerprimes[l];
  vecs:=S!.strongvecs;
  preS := false;
  while S <> false do
      o := S!.orb;
      if IsObjWithMemory(x) then
        p := o!.op(o[1],x!.el);
      else
        p := o!.op(o[1],x);
      fi;
      po := Position(o,p);
      if po = fail then   # not in current stabilizer
          Error("element not in group");
          return rec( v:=v,rem := x, S := S );
      fi;
      # Now sift through Schreier tree:
      while po > 1 do
        r:=o!.schreiergen[po];
        x := x * S!.orb!.gensi[r];
        po := o!.schreierpos[po];
        if S!.gensilayers[r]=l then
          # generator on the desired layer -- add up vectors
          v:=(v+vecs[S!.gensistrongpos[r]]) mod prime;
        fi;
      od;
      S := S!.stab;
  od;
  return v;
end ;

MatPcgsExponentsOld:=function(S,laynums,x)
local
  v,  # the accumulated full exponent vector (return value)
  a,  # the exponents contributed by one layer (from MatPcgsSift)
  i,  # loop index over the requested layers `laynums`
  j;  # loop index over the entries of `a` (used to divide off the found part)
  v:=[];
  for i in [1..Length(laynums)] do
    #S:=stabs[laynums[i]];
    a:=MatPcgsSift(S,x,laynums[i]);
    Append(v,a);
    if i<Length(laynums) then
      # need to divide off what is given by the vector
      for j in [1..Length(a)] do
        if a[j]<>0 then
          x:=LeftQuotient(S!.pcgs[S!.layranges[laynums[i]][j]]^a[j],x);
        fi;
      od;
    fi;
  od;
  return v;
end;


# new routine, using orbit positions
#decomposition as product, but not in canonical order, thus in multi stages
#TODO compare cost matrix multiplication vs. collection
MatPcgsExponents:=function(arg)
local
  # -- inputs (from arg) --
  pcgs,     # arg[1]: the pcgs
  laynums,  # arg[2]: the layers whose exponents are wanted
  ox,       # arg[3]: the element to decompose (the "original x", left untouched across passes)
  deponly,  # arg[4]: flag -- only the depth (leading nonzero) is needed, not the full exponents
  # -- chain / pcgs data --
  pS,       # the stabilizer chain, pcgs!.stabilizerChain
  dep,      # the layer boundaries, IndicesEANormalSteps(pcgs)
  md,       # the lowest layer depth that is being ignored (last entry of laynums)
  z,        # zero-vector template of full pcgs length
  # -- decomposition state --
  x,        # the working element, progressively divided down
  S,        # the current stabilizer-chain node
  e,        # the accumulated exponent vector (return value)
  delta,    # the exponent contributions of the current pass (or fail to signal termination)
  layer,    # the topmost layer currently being decomposed
  prime,    # the relative order of `layer`
  curran,   # the index range of `layer`, pcgs!.layranges[layer]
  ind,      # index of the current chain level (stabilizer depth)
  sel,      # the pcgs indices decomposed at the current base point, pcgs!.levp[ind]
  bs,       # the block sizes for the current base point, pcgs!.blocksz{sel}
  seli,     # the current pcgs index, sel[i]
  o,        # the orbit at the current chain level
  p,        # the image point o!.op(o[1],x)
  po,       # position of `p` in the orbit (minus 1)
  q,        # quotient QuoInt(po, bs[i])
  rem,      # remainder po - q*bs[i]
  # -- loop counter --
  i,        # loop index
  # -- vestigial / unused --
  prd,      # assigned [] but only used in commented-out code
  pos,      # assigned [] but never read
  preS;     # assigned false but never read

  pcgs:=arg[1];
  laynums:=arg[2];
  ox:=arg[3];
  deponly:=Length(arg)>3 and arg[4];

  dep:=IndicesEANormalSteps(pcgs);
  md:=laynums[Length(laynums)]; #the layer depth which we ignore.

  z:=ListWithIdenticalEntries(dep[Length(dep)]-1,0);
  e:=ShallowCopy(z);
  pS:=pcgs!.stabilizerChain;

  repeat
    # factor the element using orbit positions  -- this is correct only for
    # the topmost layer used
    x:=ox;
    S:=pS;
    prd:=[];
    pos:=[];
    preS := false;
    ind:=1;
    layer:=md; # Maximal layer we do
    delta:=ShallowCopy(z); # we need to forget all in the previous layer
    curran:=pcgs!.layranges[layer];
    prime:=RelativeOrders(pcgs)[curran[1]];
    while S <> false do
      sel:=pcgs!.levp[ind];
      bs:=pcgs!.blocksz{sel};

      o := S!.orb;
      if IsObjWithMemory(x) then
        p := o!.op(o[1],x!.el);
      else
        p := o!.op(o[1],x);
      fi;
      po := Position(o,p)-1;

      # decompose
      for i in [1..Length(sel)] do
        seli:=sel[i];

        q:=QuoInt(po,bs[i]);
        rem:=po-q*bs[i];
        if q>0 then
          if seli<dep[layer] then
            # decrease layer in which we decompose
            layer:=layer-1;
            while seli<dep[layer] do
              layer:=layer-1;
            od;
            delta:=ShallowCopy(z); # we need to forget all in the previous layer
            prime:=RelativeOrders(pcgs)[seli];
            curran:=pcgs!.layranges[layer];
          fi;
          #Add(prd,[sel[i],q]);;

          # remember exponent entry if right layer
          if seli<dep[layer+1] then
            delta[seli]:=delta[seli]+q mod prime;
          fi;

          #x:=x/(pcgs[sel[i]]^q);
          x:=x*pcgs!.invpows[seli][q]; # use stored inverse powers

        fi;
        po:=rem;
      od;

      #if o[1]<>o!.op(o[1],x) then Error("not all off!");fi;

      S := S!.stab;
      ind:=ind+1;
    od;

    if delta<>fail then
      # now the affected (topmost) layer is correct
      e:=e+delta;

      # we don't need to correct the last layer we use
      if deponly and not IsZero(delta) then
        delta:=fail; # pop out.
      elif layer<laynums[Length(laynums)] then
        # we have d describe a the highest level. So we want x=d*newx
        for i in curran do
          if delta[i]>0 then
            ox:=pcgs!.invpows[i][delta[i]]*ox;
          fi;
        od;
      fi;

    fi;

  # stop when we've set the lowest layer or if the remainder is the identity
  until layer=laynums[Length(laynums)] or delta=fail;

  # do we only want some?
  if laynums=pcgs!.laynums then
    return e;
  else
    return e{Concatenation(pcgs!.layranges{laynums})};
  fi;

end;

InstallMethod(ExponentsOfPcElement,"matrix pcgs",IsCollsElms,
  [IsPcgsMatGroupByStabChainRep and IsPcgs and IsPrimeOrdersPcgs,
   IsMatrixOrMatrixObj],
function(pcgs,x)
  return MatPcgsExponents(pcgs,pcgs!.laynums,x);
end);

InstallOtherMethod(ExponentsOfPcElement,"matrix pcgs,indices",IsCollsElmsX,
  [IsPcgsMatGroupByStabChainRep and IsPcgs and IsPrimeOrdersPcgs,
   IsMatrixOrMatrixObj,IsList],
function(pcgs,x,inds)
local
  lay,   # the set of layers that need to be computed
  sel,   # the positions within those layers' exponents that correspond to the requested indices
  ip,    # pointer into the requested index list `inds`
  sl,    # running offset accumulated as layers are added to `sel`
  r,     # the index range of the current layer, pcgs!.layranges[i]
  use,   # boolean: whether the current layer contributes any requested index
  i,     # loop index over the layers, pcgs!.laynums
  j;     # loop index over the entries of the current layer range `r`
  # is it a layer?
  for i in pcgs!.laynums do
    if inds=pcgs!.layranges[i] then
      # only this layer
      return MatPcgsExponents(pcgs,[i],x);
    fi;
  od;
  # which layers do we need?
  if not IsSortedList(inds) then
    # perverse case -- old
    return MatPcgsExponents(pcgs,pcgs!.laynums,x){inds};
  fi;
  sel:=[];
  lay:=[];
  ip:=1;
  sl:=0;
  for i in pcgs!.laynums do
    r:=pcgs!.layranges[i];
    use:=false;
    for j in [1..Length(r)] do
      if ip<=Length(inds) and  r[j]=inds[ip] then
        Add(sel,j+sl);
        use:=true;
        ip:=ip+1;
      fi;
    od;
    if use then
      AddSet(lay,i);
      sl:=sl+Length(r);
    fi;
  od;
  return MatPcgsExponents(pcgs,lay,x){sel};
end);

InstallMethod(DepthOfPcElement,"matrix pcgs",IsCollsElms,
  [IsPcgsMatGroupByStabChainRep and IsPcgs and IsPrimeOrdersPcgs,
  IsMatrixOrMatrixObj],
function(pcgs,x)
local
  e;  # the (depth-only) exponent vector; the depth is its first nonzero position
  e:=MatPcgsExponents(pcgs,pcgs!.laynums,x,true);
  return PositionNonZero(e);
end);

BindGlobal("TFDepthLeadExp",function(pcgs,x)
local
  p;  # position of the first nonzero exponent (the depth), used to read off the leading exponent
  x:=MatPcgsExponents(pcgs,pcgs!.laynums,x,true); # first nonzero
  p:=PositionNonZero(x);
  if p>Length(x) then
    # identity
    return [p,0];
  else
    return [p,x[p]];
  fi;
end);

InstallMethod(FittingFreeLiftSetup,"fields, using recognition",true,
  [IsMatrixGroup],0,
function(G)
local
  # -- recognition / factor --
  csi,        # the recognition tree information, GetInformationFromRecog(r)
  factorhom,  # the homomorphism onto the fitting-free factor group, FindAct(csi)
  r,          # WARNING: two roles -- (1) the default field of the matrix group, then (2) the recognition
              #          result RecognizeGroup(G)
  sz,         # the order of the radical, Size(G)/Size(Image(factorhom))
  # -- radical / solvable part --
  k,          # generators of the kernel (radical) collected via the cokernel iterator
  it,         # iterator over cokernel generators of the (inverse) factor homomorphism
  x,          # an element drawn from the cokernel iterator
  stop,       # loop flag: stop collecting kernel generators
  sbs,        # the solvable-radical BSGS setup record (SolvableBSGS result), extended and returned
  rad,        # the solvable radical subgroup of G
  pc,         # the pc group built from the radical's pcgs, PcGroupWithPcgs(sbs.pcgs)
  hom,        # the isomorphism radical -> pc group
  # -- loop counter --
  i;          # loop index
  r:=DefaultFieldOfMatrixGroup(G);
  if not IsField(r) then
    TryNextMethod();
  fi;

  r:=RecognizeGroup(G);
  SetSize(G,Size(r));
  csi:=GetInformationFromRecog(r);
  G!.storedrecog:=csi;
  factorhom:=FindAct(csi);


  #TODO: Better kernel gens by random selection
  sz:=Size(G)/Size(Image(factorhom));
  if Size(Image(factorhom))=1 then
    # solvable
    k:=GeneratorsOfGroup(G);
  else
    it:=CoKernelGensIterator(InverseGeneralMapping(factorhom));
    k:=Filtered(List([1..csi.genum],x->CSINiceGens(csi,x)),x->IsOne(ImagesRepresentative(factorhom,x)));
    if sz>1 then
      stop:=true;
      repeat
        for i in [1..3*Length(Factors(sz))] do
          if not IsDoneIterator(it) then
            x:=NextIterator(it);
            if not IsOne(x) and not x in k then
              Add(k,x);
            fi;
          fi;
        od;
        if sz>1 and Length(k)<Length(Factors(sz)) then stop:=false; fi; # work around issue
        if ValueOption("doall")=true then stop:=false;fi;
      until stop or IsDoneIterator(it);
    fi;
  fi;

  Info(InfoGenSS,3,"|k|=",Length(k));
  # TODO: Proper test

  if ForAll(k,IsOne) then
    rad:=TrivialSubgroup(G);
    sbs:=rec(pcgs:=[],depths:=[1],relord:=[],pcisom:=IsomorphismPcGroup(rad),radical:=rad);
  else
    sbs:=SolvableBSGS(G,k,sz);
    while sbs=fail do
      sbs:=SolvableBSGS(G,k,sz:doall);
    od;
    rad:=SubgroupNC(G,sbs.pcgs);
    SetSize(rad,Product(RelativeOrders(sbs.pcgs)));
    sbs.radical:=rad;
    pc:=PcGroupWithPcgs(sbs.pcgs);
    if RUN_IN_GGMBI=false then
    RUN_IN_GGMBI:=true; # hack to skip Nice treatment
    hom:=GroupHomomorphismByImagesNC(rad,pc,sbs.pcgs,FamilyPcgs(pc));
    RUN_IN_GGMBI:=false;
    else
      hom:=GroupHomomorphismByImagesNC(rad,pc,sbs.pcgs,FamilyPcgs(pc):Run_In_GGMBI:= true);
    fi;
    SetIsBijective(hom,true);
    sbs.pcisom:=hom;
  fi;
  sbs.csi:=csi;
  SetKernelOfMultiplicativeGeneralMapping(factorhom,rad);
  sbs.factorhom:=factorhom;
  RecogDecompinfoHomomorphism(factorhom)!.LiftSetup:=sbs;

  return sbs;
end);

FFStats:=function(g)
local
  start,  # the runtime at the start, used to report elapsed time
  f;      # the fitting-free lift setup of `g`, FittingFreeLiftSetup(g)
  start:=Runtime();
  f:=FittingFreeLiftSetup(g);
  Print("Time:",Runtime()-start,"\n");
  Print("Factordegree ",NrMovedPoints(Range(f.factorhom)),"\n");
  Print("PcgsDegrees ",Maximum(List(
    BasePointsActionsOrbitLengthsStabilizerChain(f.pcgs!.stabilizerChain),
    x->x[3])),"\n");

end;

# work over residue class rings

MyZmodnZObj:=function(a,b)
  if IsPrimeInt(b) and b<65536 then
    return Z(b)^0*a;
  else
    return ZmodnZObj(a,b);
  fi;
end;

ReduceModM:=function(a,m)
  local
    b,  # the reduced matrix or vector being assembled (return value)
    r,  # WARNING: two roles -- (1) a single reduced row while building a matrix, and
        #          (2) the residue ring `Integers mod m` in the matrix-object branch
    i,  # loop index over rows
    j;  # loop index over entries
  if IsList(a) then
    if IsList(a[1]) then
      # matrix
      b:=[];
      for i in a do
        r:=[];
        for j in i do
          if IsZmodnZObjNonprime(j) then
            Add(r,MyZmodnZObj(j![1],m));
          else
            Add(r,MyZmodnZObj(j,m));
          fi;
        od;
        Add(b,r);
      od;
      if IsPrimeInt(m) then
        b:=ImmutableMatrix(m,b);
      fi;
      return b;
    else
      # vector
      b:=[];
      for j in a do
        if IsZmodnZObjNonprime(j) then
          Add(b,MyZmodnZObj(j![1],m));
        else
          Add(b,MyZmodnZObj(j,m));
        fi;
      od;
      return b;
    fi;
  elif IsMatrixObj(a) then
    a:=List(a,x->List(x,Int));
  else
    Error("weird object");
  fi;
  r:=Integers mod m;
  if IsPrimeInt(m) and m<=65536 then
    return ImmutableMatrix(r,a*One(r));
  else
    return Matrix(r,a*One(r));
  fi;
end;

ReduceModMFunc:=function(m)
  return a->ReduceModM(a,m);
end;

UnreduceModM:=Error;

InstallMethod(FittingFreeLiftSetup,"residue class rings",true,
  [IsMatrixGroup],0,
function(g)
local
  # -- local helper functions --
  elmimg,        # helper: map an element to its image in the direct product of the per-prime factor images
  addCleanUpper, # helper: clean an element in the induced (quotient) pcgs of factor i
  addPcElement,  # helper: add a pc element into the per-level module bases; returns the first level changed
  procrels,      # helper: process the conjugation relations in the linear layers
  # -- group / ring setup --
  r,             # WARNING: two roles -- (1) the coefficient ring DefaultFieldOfMatrixGroup(g), and
                 #          (2) the decomposition-info record assembled and returned at the end
  triv,          # a trivial group, used as the range for trivial factor homomorphisms
  idmat,         # the identity matrix over the ring (used to subtract off the identity part)
  gens,          # the generators of g, converted to compact matrix objects if necessary
  gnew,          # the group on the converted generators (equal to g if no conversion was needed)
  m,             # the modulus, Size(r)
  f,             # the list of prime-power factors of `m`
  # -- per-prime factor decomposition --
  homs,          # list of the reduction homomorphisms (mod p) for each prime power
  hom,           # WARNING: two roles -- (1) a temporary list of reduced generator images, and
                 #          (2) a genuine group homomorphism (the reduction / overall factor map)
  img,           # the image group of a reduction homomorphism, Group(hom)
  ff,            # the FittingFreeLiftSetup of a factor image
  ffp,           # list of the factor pcgs
  ffppc,         # list of the pcisom maps of the factors
  ffhoms,        # list of the factor homomorphisms (reduction composed with the factor's factorhom)
  d,             # WARNING: two roles -- (1) the list of factor images, and (2) their DirectProduct group
  # -- module-level tables (indexed by module level) --
  moli,          # the moduli (prime powers) of the successive module levels
  pli,           # the prime of each module level
  idx,           # for each level, the index of its prime in `pli`
  fac,           # for each level, the scaling factor (1/p^k) applied when extracting the residue
  # -- induced pcgs / factor state --
  ffpi,          # per factor: the induced pcgs computed anew
  ffsubs,        # per factor: the subgroup already covered by the induced pcgs
  upperpcgs,     # per factor: preimages (in g) of the induced-pcgs generators
  # -- per-level bases of the linear layers --
  bas,           # per level: an echelonized basis of residue vectors
  basrep,        # per level: group-element representatives of the basis vectors
  basrepi,       # per level: inverses of those representatives
  # -- iteration and work-lists --
  it,            # iterator over cokernel generators of the (inverse) factor homomorphism
  layerlimit,    # upper bound on the dimension of each linear layer
  stacks,        # set of elements already processed (avoids reprocessing)
  stack,         # work-list of G-conjugates still to process
  ao,            # the original element, saved before it is reduced/cleaned
  bl,            # saved base lengths per level, to detect when a base actually grew
  fertig,        # "done" flag: all layers are full, so processing can stop  (German: finished)
  # -- final pcgs assembly --
  pcgs,          # the accumulated pc sequence, later the Pcgs object
  relord,        # the relative orders of the pcgs
  levs,          # the EA-normal-step (layer) indices of the pcgs
  s,             # WARNING: two roles -- (1) EA-normal-step index lists while assembling levels, and
                 #          (2) the radical subgroup SubgroupNC(g,pcgs) at the end
  # -- heavily reused scratch variables --
  a,             # WARNING: multi-role -- a factor list, generator images, an iterator element, and a
                 #          conjugate/power element (reassigned throughout)
  p,             # WARNING: multi-role -- a prime, a running count, the pc group, and a homomorphism
  # -- loop counters --
  i,             # loop index
  j,             # loop index
  k,             # loop index / element from `stack`
  l;             # loop index

  triv:=TrivialSubgroup(CyclicGroup(2));
  r:=DefaultFieldOfMatrixGroup(g);
  if IsField(r) or
    not CategoryCollections(IsZmodnZObjNonprime)(r) then
    TryNextMethod();
  fi;
  idmat:=IdentityMat(Length(One(g)),1);

  # convert to compact type
  gens:=GeneratorsOfGroup(g);
  if not ForAll(gens,IsMatrixObj) then
    gens:=List(gens,x->Matrix(r,One(r)*List(x,r->List(r,Int))));
    gnew:=Group(gens);
    if HasSize(g) then SetSize(gnew,Size(g));fi;
  else
    gnew:=g;
  fi;

  m:=Size(r);
  # the prime power factors occurring
  f:=List(Collected(Factors(m)),x->x[1]^x[2]);
  Sort(f);
  homs:=[];
  ffp:=[];
  ffppc:=[];
  ffhoms:=[];
  moli:=[1];
  pli:=[];
  idx:=[false];
  fac:=[false];
  #idmats:=[];
  for i in f do
    a:=Factors(i);
    p:=a[1];
    if Length(a)>1 then
      Add(moli,moli[Length(moli)]*p^2);
      Add(fac,1/p);
      Add(pli,p);
      Add(idx,Length(pli));
      for j in [3..Length(a)] do
        Add(moli,moli[Length(moli)]*p);
        Add(idx,Length(pli));
        Add(fac,fac[Length(fac)]/p);
      od;
    fi;
    hom:=List(gens,x->ReduceModM(x,p));
    if ForAll(hom,IsOne) then
      hom:=GroupHomomorphismByFunction(gnew,Group(()),x->());
      Info(InfoFFMat,2,"alltrivial ",p);
      Add(ffp,[]);
      Add(homs,hom);
      Add(ffhoms,hom);
      hom:=GroupHomomorphismByFunction(Group(()),triv,x->One(triv));
      Add(ffppc,hom);
    else
      img:=Group(hom);
      hom:=GroupHomomorphismByFunction(gnew,img,ReduceModMFunc(p),false,
            x->UnreduceModM(x,m));
      SetImagesSource(hom,img);
      Add(homs,hom);
      ff:=FittingFreeLiftSetup(img);
      Add(ffp,ff.pcgs);
      Add(ffppc,ff.pcisom);
      Add(ffhoms,hom*ff.factorhom);
    fi;
  od;
  if Length(ffhoms)=0 then
    hom:=GroupHomomorphismByFunction(gnew,Group(()),x->());
  else
    d:=List(ffhoms,Image);
    d:=DirectProduct(d);
    elmimg:=function(x)
     local p,  # the accumulated product image in the direct product `d`
       i;  # loop index over the factor homomorphisms `ffhoms`
       p:=One(d);
       for i in [1..Length(ffhoms)] do
          p:=p*ImagesRepresentative(Embedding(d,i),
            ImagesRepresentative(ffhoms[i],x));
       od;
       return p;
     end;
    a:=List(gens,elmimg);
    hom:=GroupHomomorphismByFunction(gnew,SubgroupNC(d,a),elmimg);
  fi;

  # we can't rescue the existing pcgs for the factors as we will not know
  # preimages. Compute all anews.
  ffpi:=List([1..Length(ffp)],x->[]);
  #$ffpi[1]:=ffp[1]; # the first pcgs is guaranteed to be always fully there

  ffsubs:=List([1..Length(ffpi)],x->
            TrivialSubgroup(Image(ffppc[x])));

  upperpcgs:=List([1..Length(ffpi)],x->[]);

  # Do we have an upper limit for the space on each layer?
  layerlimit:=ValueOption("layerlimit");
  if layerlimit=fail then
    layerlimit:=Length(One(g))^2; # full matrix space
  fi;

  it:=CoKernelGensIterator(InverseGeneralMapping(hom));
  bas:=List(moli,x->[]);
  basrep:=List(moli,x->[]);
  basrepi:=List(moli,x->[]);

  addCleanUpper:=function(i,a)
  local
    r,  # the image of `a` in factor i, ImagesRepresentative(homs[i],a)
    p,  # the image of `r` in the pc factor (tests whether the induced pcgs must be extended)
    s,  # the extended induced pcgs together with its preimages (a [pcgs,preimages] pair)
    e;  # exponents of `r`, then the corresponding linear combination of upper-pcgs preimages
    r:=ImagesRepresentative(homs[i],a);
    p:=ImagesRepresentative(ffppc[i],r);
    if not p in ffsubs[i] then
      # need to extend induced pcgs
      s:=CanonicalPcgsByGeneratorsWithImages(ffp[i],
            Concatenation(ffpi[i],[r]),
            Concatenation(upperpcgs[i],[a]));
      ffpi[i]:=s[1];
      upperpcgs[i]:=s[2];
      ffsubs[i]:=ClosureGroup(ffsubs[i],p);
    fi;
    # divide off the nontrivial part in the quotient
    if Length(ffpi[i])>0 then
      e:=ExponentsOfPcElement(ffpi[i],r);
      e:=LinearCombinationPcgs(upperpcgs[i],e);
      a:=LeftQuotient(e,a);
    fi;
    return a;
  end;

  addPcElement:=function(a,start)
    local
      added,  # the first module level at which a new basis vector was added (return value)
      bot,    # boolean: whether we are at the bottom (last) module level
      p,      # the prime of the current level, pli[idx[i]]
      e,      # the residue vector of `a` at the current level
      s,      # solution of `e` over the current basis (fail if `e` is independent)
      i,      # loop index over the module levels
      j;      # loop index over the entries of the solution `s`

      added:=fail;
      for i in [start..Length(moli)] do
        bot:=i=Length(moli); # last step -- powers are multiples
        p:=pli[idx[i]];
        e:=List(a,x->List(x,y->y![1] mod moli[i]));
        e:=e-idmat;
        e:=fac[i]*Concatenation(e);
        e:=ImmutableVector(GF(p),e*Z(p)^0);
        if not IsZero(e) then
          if bot and IsBound(bas[i]) and Length(bas[i])=layerlimit then
            # Bottom layer is full, nothing else needed to do
            a:=fail;
          elif IsBound(bas[i]) and Length(bas[i])>0 then
            s:=SolutionMat(bas[i],e);
            if s=fail then
              Add(bas[i],e);
              Add(basrep[i],a);
              Add(basrepi[i],a^-1);
              if added=fail then added:=i;fi;
              if not bot then a:=a^p;fi; # no need to do if on bottom
            else
              s:=List(s,Int);
              # avoid inverse by imediately dividing off
              for j in [Length(s),Length(s)-1..1] do
                 if not IsZero(s[j]) then
                  if not bot then
                    a:=basrepi[i][j]^s[j]*a;
                  else
                    # else case is not needed, as we are on the bottom :-)
                    a:=fail;
#                    # multiplication is addition of p-residues
#                    tmp:=(basrepi[i][j]![1]*s[j]-s[j]*One(basrepi[i][j]![1])+a![1]);
#                    tmp:=tmp mod Characteristic(a);
#                    tmp:=MakeZmodnZMat(ElementsFamily(ElementsFamily(FamilyObj(a))),tmp);
#                    a:=tmp;
                  fi;
                fi;
              od;
            fi;
          else
            bas[i]:=[e];
            basrep[i]:=[a];
            basrepi[i]:=[a^-1];
            if added=fail then added:=i;fi;
            if not bot then a:=a^p;fi; # no need to do if on bottom
          fi;
        fi;
      od;
      return added;
    end;

  fertig:=false;
  stacks:=[];
  repeat
    a:=NextIterator(it);
    if not IsOne(a) then
      ao:=a;

      for i in [1..Length(ffpi)] do
        a:=addCleanUpper(i,a);
      od;

      bl:=List([2..Length(moli)],x->Length(bas[x]));
      addPcElement(a,2);
      if ForAny([2..Length(moli)],x->Length(bas[x])>bl[x-1]) then
        AddSet(stacks,ao);
      fi;

    fi;
    fertig:=Length(moli)>1 and ForAll([2..Length(moli)],x->Length(basrep[x])=layerlimit);
  until IsDoneIterator(it) or fertig;

  if fertig then stack:=[];fi;
  Info(InfoFFMat,2,"layerdimsc:",List(bas,Length));

  stack:=ShallowCopy(stacks); # we'll add to the list

  # G-conjugates
  for k in stack do
    for j in gens do
      a:=k^j;

      for i in [1..Length(ffpi)] do
        a:=addCleanUpper(i,a);
      od;

      if not a in stacks then

        AddSet(stacks,a);

        # old base length
        bl:=List([2..Length(moli)],x->Length(bas[x]));
        addPcElement(a,2);
        if ForAny([2..Length(moli)],x->Length(bas[x])>bl[x-1]) then
          if not a in stack then Add(stack,a); fi;
        fi;
      fi;

    od;
  od;
  Info(InfoFFMat,2,"layerdimsb:",List(bas,Length));

  Unbind(stack);
  Unbind(stacks);

  pcgs:=[];
  relord:=[];
  levs:=[];
  p:=0;
  for i in [1..Length(ffp)] do
    if Length(upperpcgs[i])>0 then

      #r:=RelativeOrders(ffpi[i]);
      if not fertig then
        j:=1;
        while j<=Length(upperpcgs[i]) do
          a:=upperpcgs[i][j]^(RelativeOrders(ffpi[i])[j]);
          for k in [i..Length(ffp)] do
            a:=addCleanUpper(k,a);
          od;
          addPcElement(a,2);
          for l in Concatenation(pcgs,upperpcgs[i]) do
            a:=upperpcgs[i][j]^l;
            for k in [i..Length(ffp)] do
              a:=addCleanUpper(k,a);
            od;
            addPcElement(a,2);
          od;
          j:=j+1;
        od;
      fi;

      Append(pcgs,upperpcgs[i]);
      if not HasIndicesEANormalSteps(ffpi[i]) then
        s:=FamilyPcgs(PcGroupWithPcgs(ffpi[i]));
        if not IsPcgsElementaryAbelianSeries(s) then Error("not EA pcgs");fi;
        s:=IndicesEANormalSteps(s);
      else
        s:=IndicesEANormalSteps(ffpi[i]);
      fi;
      s:=s{[1..Length(s)-1]}+p;
      Append(levs,s);
      Append(relord,RelativeOrders(ffpi[i]));
      p:=Length(pcgs);
    fi;
  od;
  Add(levs,Length(pcgs)+1);

  Info(InfoFFMat,2,"layerdimsa:",List(bas,Length));

  # conjugation relations in linear bits. Will rarely add new elements (so the
  # danger of running through multiple times is minimal).
  procrels:=function()
  local
    b,  # boolean return value: true means no new element was added (the process has converged)
    a,  # the result of addPcElement (the level changed, or fail if nothing was added)
    k,  # the conjugating element (a basis representative or a pcgs element)
    l,  # the element being conjugated (a basis representative)
    i,  # loop index over module levels
    j;  # loop index over module levels
    b:=true;
    for i in [2..Length(moli)] do
      # if j>=Length(moli)-1+@ the conjugate is guaranteed the same
      for j in [i..Length(moli)-i+1] do
        if IsBound(basrep[i]) and IsBound(basrep[j]) then
          for k in basrep[i] do
            for l in basrep[j] do
              a:=addPcElement(l^k,j);
              if a<>fail then b:=false;fi;
            od;
          od;
        fi;
      od;
    od;
    for j in [2..Length(moli)] do
      if IsBound(basrep[j]) then
        for k in pcgs do
          for l in basrep[j] do
            a:=addPcElement(l^k,j);
            if a<>fail then b:=false;fi;
          od;
        od;
      fi;
    od;
    return b;
  end;

  repeat until fertig or procrels();

  Info(InfoFFMat,2,"layerdims:",List(bas,Length));

  for i in [2..Length(bas)] do
    if IsBound(bas[i]) then
      Append(pcgs,basrep[i]);
      Append(relord,ListWithIdenticalEntries(Length(basrep[i]),pli[idx[i]]));
    fi;
    Add(levs,Length(pcgs)+1);
  od;

  pcgs:=PcgsByPcSequenceCons(IsPcgsDefaultRep,
    IsPcgsResidueMatGroupRep and IsPcgs and IsPrimeOrdersPcgs,
    FamilyObj(One(g)),pcgs,[]);

  # decomposition info for pcgs
  r:=rec(pcgs:=pcgs,factorhom:=hom,
          homs:=homs,ffp:=ffp,ffpi:=ffpi,ffppc:=ffppc,idmat:=idmat,
          upperpcgs:=upperpcgs,moli:=moli,pli:=pli,idx:=idx,fac:=fac,
          depths:=levs,bas:=bas,basrep:=basrep);

  pcgs!.decompInfo:=r;
  SetRelativeOrders(pcgs,relord);
  SetIndicesEANormalSteps(pcgs,levs);

  s:=SubgroupNC(g,pcgs);
  SetSize(s,Product(RelativeOrders(pcgs)));
  SetSize(g,Size(Image(hom))*Size(s));
  SetKernelOfMultiplicativeGeneralMapping(hom,s);

  r:=rec(pcgs:=pcgs,radical:=s,factorhom:=hom,depths:=levs);

  if ValueOption("pcisom")<>false then
    i:=List(pcgs,x->ExponentsOfPcElement(pcgs,x^-1));
    p:=PcGroupWithPcgs(pcgs:inversehints:=i);
    if RUN_IN_GGMBI=false then
    RUN_IN_GGMBI:=true; # hack to skip Nice treatment
    p:=GroupHomomorphismByImagesNC(s,p,pcgs,FamilyPcgs(p));
    RUN_IN_GGMBI:=false;
    else
      p:=GroupHomomorphismByImagesNC(s,p,pcgs,FamilyPcgs(p):Run_In_GGMBI:= true);
    fi;
    SetIsBijective(p,true);
    r.pcisom:=p;
  fi;

  SetRecogDecompinfoHomomorphism(hom,rec(fct:=elmimg,LiftSetup:=r));
  return r;

end);

ExponentsResiduePcgs:=function(r,elm)
local
  exp,  # the accumulated exponent vector (return value)
  e,    # the residue vector of the current linear layer
  p,    # the prime of the current layer, r.pli[r.idx[i]]
  s,    # WARNING: two roles -- (1) an exponent list (from ExponentsOfPcElement / a matrix solution), and
        #          (2) the corresponding group element (a linear combination of pcgs elements)
  i;    # loop index (over the upper pcgs factors, then over the linear layers)
  exp:=[];
  # upper pcgs part
  for i in [1..Length(r.ffpi)] do
    if Length(r.ffpi[i])>0 then
      s:=ExponentsOfPcElement(r.ffpi[i],ImagesRepresentative(r.homs[i],elm));
      Append(exp,s);
      elm:=LeftQuotient(LinearCombinationPcgs(r.upperpcgs[i],s),elm);
    fi;
  od;

  for i in [2..Length(r.moli)] do
    p:=r.pli[r.idx[i]];
    e:=List(elm,x->List(x,y->y![1] mod r.moli[i]));
    e:=e-r.idmat;
    e:=r.fac[i]*Concatenation(e);
    e:=e*Z(p)^0;
    if not IsZero(e) then
      s:=SolutionMat(r.bas[i],e);
      if s=fail then
        Error("not in span of pcgs");
      else
        s:=List(s,Int);
        Append(exp,s);
        s:=LinearCombinationPcgs(r.basrep[i],s);
        elm:=LeftQuotient(s,elm);
      fi;
    elif IsBound(r.bas[i]) then
      Append(exp,ListWithIdenticalEntries(Length(r.bas[i]),0));
    fi;
#Print(exp);
  od;

  return exp;
end;

InstallMethod(ExponentsOfPcElement,"matrix residue pcgs",IsCollsElms,
  [IsPcgsResidueMatGroupRep and IsPcgs and IsPrimeOrdersPcgs,
    IsMatrixOrMatrixObj], function(pcgs,x)
  return ExponentsResiduePcgs(pcgs!.decompInfo,x);
end);

#############################################################################
##
#M  RestrictedMapping(<hom>,<U>)
##
InstallMethod(RestrictedMapping,"for recognition mappings",
  CollFamSourceEqFamElms,[IsGroupGeneralMapping and
  HasRecogDecompinfoHomomorphism,IsGroup],
  # make this ranked higher than even some default subset methods in grp.gi
  SUM_FLAGS+100,
function(hom, U)
local
  d,     # the decomposition-info record of `hom`, RecogDecompinfoHomomorphism(hom)
  gens,  # the generators of the subgroup U
  imgs,  # the images of those generators under the restriction
  rest;  # the restricted homomorphism (return value)
  d:=RecogDecompinfoHomomorphism(hom);
  gens:=GeneratorsOfGroup(U);
  if IsBound(d.fct) then
    imgs:=List(gens,d.fct);
  elif IsBound(d.isTrivial) then
    imgs:=List(gens,x->());
  else
    imgs:=List(gens,x->CSIElementAct(d,x));
  fi;

  if IsPermGroup(U) and AssertionLevel()>2 then
    rest:=GroupHomomorphismByImages(U,Range(hom),gens,imgs);
  else
    if RUN_IN_GGMBI=false then
    RUN_IN_GGMBI:=true; # hack to skip Nice treatment
    if IsBound(d.fct) then
      rest:=GroupHomomorphismByFunction(U,Range(hom),d.fct);
    else
      rest:=GroupHomomorphismByImagesNC(U,Range(hom),gens,imgs);
    fi;
    RUN_IN_GGMBI:=false;
    else
      PushOptions( rec( Run_In_GGMBI:= true ) ); # hack to skip Nice treatment
      if IsBound(d.fct) then
        rest:=GroupHomomorphismByFunction(U,Range(hom),d.fct);
      else
        rest:=GroupHomomorphismByImagesNC(U,Range(hom),gens,imgs);
      fi;
      PopOptions();
    fi;



  fi;

  SetRecogDecompinfoHomomorphism(rest,d);
  return rest;
end);

# temporary -- missing in library
InstallMethod(MaximalSubgroupClassReps,"TF method",true,
  [IsGroup and IsFinite and HasFittingFreeLiftSetup],
  {} -> RankFilter( IsHandledByNiceMonomorphism ),
DoMaxesTF);

TFSubgroupMembership:=function(g,u,elm)
local
  ffu,  # the fitting-free setup for the subgroup u inside g, FittingFreeSubgroupSetup(g,u)
  ff,   # the parent fitting-free setup, ffu.parentffs
  x;    # the image of `elm` in the fitting-free factor
  ffu:=FittingFreeSubgroupSetup(g,u);
  ff:=ffu.parentffs;
  x:=ImagesRepresentative(ff.factorhom,elm);
  if not x in Image(ffu.rest) then
    return false;
  fi;
  elm:=elm/PreImagesRepresentativeNC(ffu.rest,x);
  elm:=ImagesRepresentative(ff.pcisom,elm);
  if not IsBound(ffu.pcsub) then
    ffu.pcsub:=Subgroup(Image(ff.pcisom),
      List(ffu.pcgs,x->ImagesRepresentative(ff.pcisom,x)));
  fi;
  return elm in ffu.pcsub;
end;

# careful -- this implicitly assumes we always stay in the parent
InstallMethod(\in,"ff subgroup",IsElmsColls,
  [IsMultiplicativeElementWithInverse,IsGroup and IsMatrixGroup],
  {} -> RankFilter( IsHandledByNiceMonomorphism ),
function(elm,u)
  if not HasParentAttr(u) or not HasFittingFreeLiftSetup(ParentAttr(u)) then
    TryNextMethod();
  fi;
  return TFSubgroupMembership(ParentAttr(u),u,elm);
end);

# generic method, avoid nice method
TFNormalClosure:=function ( G, N )
local   gensG,      # generators of the group <G>
        genG,       # one generator of the group <G>
        gensN,      # generators of the group <N>
        genN,       # one generator of the group <N>
        cnj;        # a conjugate genN^genG of a generator of
                    # <N> by a generator of <G>

  gensG := GeneratorsOfGroup( G );
  gensN := ShallowCopy( GeneratorsOfGroup( N ) );
  for genN  in gensN  do
    for genG  in gensG  do
      cnj := genN ^ genG;
      if not TFSubgroupMembership(G,N,cnj) then
        Add( gensN, cnj );
        N := SubgroupNC(G,gensN);
      fi;
    od;

  od;

  # return the normal closure
  return N;

end;

