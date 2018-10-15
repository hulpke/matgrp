#############################################################################
##
#W  zmodmat.g                   matgrp package               Alexander Hulpke
##
##
#Y  Copyright (C)  2018, Alexander Hulpke
##
##  Implement a data type for matrices over residue class rings
##

#############################################################################
##
#R  IsZmodnZMatRep( <obj> )
##
##  Representations for matrices and vectors
##
DeclareRepresentation( "IsZmodnZMatRep", IsPositionalObjectRep, [ 1 ] );
DeclareRepresentation( "IsZmodnZVecRep", IsPositionalObjectRep, [ 1 ] );

IsZmodnZMat:=IsZmodnZObjNonprimeCollColl and IsAssociativeElement and IsAdditivelyCommutativeElement
      and IsZmodnZMatRep and IsMultiplicativeElementWithInverse and IsList and
      IsDenseList and IsTable and IsExtAElement and IsHomogeneousList and IsMatrix;

IsZmodnZVec:=IsZmodnZObjNonprimeCollection and IsAssociativeElement
      and IsAdditivelyCommutativeElement
      and IsZmodnZVecRep and IsList and IsDenseList and IsExtAElement and IsHomogeneousList;

InstallTrueMethod(CanEasilySortElements,IsZmodnZMat);
InstallTrueMethod(CanEasilySortElements,IsZmodnZVec);

MakeZmodnZMat:=function(elmfam,mat)
  if not IsBound(elmfam!.mattype) then
    elmfam!.mattype:=NewType(CollectionsFamily(CollectionsFamily(elmfam)),IsZmodnZMat);
  fi;
  return Objectify(elmfam!.mattype,[mat]);
end;

ZmodnZMat:=function(r,mat)
local fam;
  fam:=FamilyObj(One(r));
  mat :=mat mod Characteristic(fam);
  return MakeZmodnZMat(fam,mat);
end;

InstallMethod(PrintObj,"ZmodnZMat",true,[IsZmodnZMat],0,
function(a)
  Print(a![1],"*ZmodnZObj(1,",Characteristic(ElementsFamily(ElementsFamily(FamilyObj(a)))),")");
end);

InstallOtherMethod(One,"ZmodnZMat",true,[IsZmodnZMat],0,
function(a)
  return MakeZmodnZMat(ElementsFamily(ElementsFamily(FamilyObj(a))),One(a![1]));
end);

InstallOtherMethod(Zero,"ZmodnZMat",true,[IsZmodnZMat],0,
function(a)
  return MakeZmodnZMat(ElementsFamily(ElementsFamily(FamilyObj(a))),Zero(a![1]));
end);

InstallMethod(\+,"ZmodnZMat",IsIdenticalObj,[IsZmodnZMat,IsZmodnZMat],0,
function(a,b)
local fam;
  fam:=ElementsFamily(ElementsFamily(FamilyObj(a)));
  b:=a![1]+b![1];
  b:=b mod Characteristic(fam);
  return MakeZmodnZMat(fam,b);
end);

InstallMethod(\=,"ZmodnZMat",IsIdenticalObj,[IsZmodnZMat,IsZmodnZMat],0,
function(a,b)
  return a![1]=b![1];
end);

InstallMethod(\<,"ZmodnZMat",IsIdenticalObj,[IsZmodnZMat,IsZmodnZMat],0,
function(a,b)
  return a![1]<b![1];
end);

InstallOtherMethod(Length,"ZmodnZMat",true,[IsZmodnZMat],0,
function(a)
  return Length(a![1]);
end);

InstallOtherMethod(DimensionsMat,"ZmodnZMat",true,[IsZmodnZMat],0,
function(a)
  return DimensionsMat(a![1]);
end);

InstallMethod(\*,"ZmodnZMat",IsIdenticalObj,[IsZmodnZMat,IsZmodnZMat],0,
function(a,b)
local fam;
  fam:=ElementsFamily(ElementsFamily(FamilyObj(a)));
  b:=a![1]*b![1];
  b:=b mod Characteristic(fam);
  return MakeZmodnZMat(fam,b);
end);

InstallOtherMethod(InverseOp,"ZmodnZMat",true,[IsZmodnZMat],0,
function(a)
local fam,d;
  fam:=ElementsFamily(ElementsFamily(FamilyObj(a)));
  d:=DeterminantMat(a![1]);
  a:=InverseOp(1/d*a![1]);
  a:=a mod Characteristic(fam);
  return MakeZmodnZMat(fam,a);
end);

InstallMethod(\^,"ZmodnZMat,int",true,[IsZmodnZMat,IsInt],0,
function(a,e)
local fam;
  if e<0 then a:=InverseOp(a);e:=-e;
  elif e=0 then return One(a);
  fi;

  fam:=ElementsFamily(ElementsFamily(FamilyObj(a)));
  a:=a![1]^e;
  a:=a mod Characteristic(fam);
  return MakeZmodnZMat(fam,a);
end);

InstallOtherMethod(AdditiveInverseOp,"ZmodnZMat",true,[IsZmodnZMat],0,
function(a)
local fam;
  fam:=ElementsFamily(ElementsFamily(FamilyObj(a)));
  a:=-a![1] mod Characteristic(fam);
  return MakeZmodnZMat(fam,a);
end);

# vectors

MakeZmodnZVec:=function(elmfam,vec)
  if not IsBound(elmfam!.vectype) then
    elmfam!.vectype:=NewType(CollectionsFamily(elmfam),IsZmodnZVec);
  fi;
  return Objectify(elmfam!.vectype,[vec]);
end;

ZmodnZVec:=function(r,vec)
local fam;
  fam:=FamilyObj(One(r));
  vec :=vec mod Characteristic(fam);
  return MakeZmodnZVec(fam,vec);
end;

InstallMethod(PrintObj,"ZmodnZVec",true,[IsZmodnZVec],0,
function(a)
  Print(a![1],"*ZmodnZObj(1,",Characteristic(ElementsFamily(FamilyObj(a))),")");
end);

InstallOtherMethod(Zero,"ZmodnZVec",true,[IsZmodnZVec],0,
function(a)
  return MakeZmodnZVec(ElementsFamily(FamilyObj(a)),Zero(a![1]));
end);

InstallMethod(\+,"ZmodnZVec",IsIdenticalObj,[IsZmodnZVec,IsZmodnZVec],0,
function(a,b)
local fam;
  fam:=ElementsFamily(FamilyObj(a));
  b:=a![1]+b![1];
  b:=b mod Characteristic(fam);
  return MakeZmodnZVec(fam,b);
end);

InstallOtherMethod(AdditiveInverseOp,"ZmodnZVec",true,[IsZmodnZVec],0,
function(a)
local fam;
  fam:=ElementsFamily(FamilyObj(a));
  a:=-a![1] mod Characteristic(fam);
  return MakeZmodnZVec(fam,a);
end);

InstallMethod(\=,"ZmodnZVec",IsIdenticalObj,[IsZmodnZVec,IsZmodnZVec],0,
function(a,b)
  return a![1]=b![1];
end);

InstallMethod(\<,"ZmodnZVec",IsIdenticalObj,[IsZmodnZVec,IsZmodnZVec],0,
function(a,b)
  return a![1]<b![1];
end);

InstallOtherMethod(Length,"ZmodnZVec",true,[IsZmodnZVec],0,
function(a)
  return Length(a![1]);
end);

InstallMethod(\[\],"ZmodnZVec",true,[IsZmodnZVec,IsPosInt],0,
function(a,p)
local fam;
  fam:=ElementsFamily(FamilyObj(a));
  return ZmodnZObj(fam,a![1][p]);
end);



InstallMethod(\[\],"ZmodnZMat",true,[IsZmodnZMat,IsPosInt],0,
function(a,p)
local fam;
  fam:=ElementsFamily(ElementsFamily(FamilyObj(a)));
  return MakeZmodnZVec(fam,a![1][p]);
end);

# mixed media
InstallOtherMethod(\+,"ZmodnZMat+List",IsIdenticalObj,[IsZmodnZMat,IsList],0,
function(a,b)
local fam;
  fam:=ElementsFamily(ElementsFamily(FamilyObj(a)));
  b:=a![1]+b;
  b:=List(b,r->List(r,Int));
  return MakeZmodnZMat(fam,b);
end);

InstallOtherMethod(\+,"List+ZmodnZMat",IsIdenticalObj,[IsList,IsZmodnZMat],0,
function(a,b)
local fam;
  fam:=ElementsFamily(ElementsFamily(FamilyObj(b)));
  b:=a+b![1];
  b:=List(b,r->List(r,Int));
  return MakeZmodnZMat(fam,b);
end);

InstallOtherMethod(\+,"intmat+ZmodnZMat",true,[IsList and IsTable,IsZmodnZMat],0,
function(a,b)
local fam;
  if not CollectionsFamily(FamilyObj(Integers))=FamilyObj(a) then
    TryNextMethod();
  fi;
  fam:=ElementsFamily(ElementsFamily(FamilyObj(b)));
  b:=a+b![1];
  b:=List(b,r->List(r,Int));
  return MakeZmodnZMat(fam,b);
end);

InstallOtherMethod(\*,"ZmodnZMat*List",IsIdenticalObj,[IsZmodnZMat,IsList],0,
function(a,b)
local fam;
  fam:=ElementsFamily(ElementsFamily(FamilyObj(a)));
  b:=a![1]*b;
  b:=List(b,r->List(r,Int));
  return MakeZmodnZMat(fam,b);
end);

InstallOtherMethod(\*,"List*ZmodnZMat",IsIdenticalObj,[IsList,IsZmodnZMat],0,
function(a,b)
local fam;
  fam:=ElementsFamily(ElementsFamily(FamilyObj(b)));
  b:=a*b![1];
  b:=List(b,r->List(r,Int));
  return MakeZmodnZMat(fam,b);
end);

InstallOtherMethod(\=,"ZmodnZMat=List",IsIdenticalObj,[IsZmodnZMat,IsList],0,
function(a,b)
local fam;
  b:=List(b,r->List(r,Int));
  return a![1]=b;
end);

InstallOtherMethod(\=,"List=ZmodnZMat",IsIdenticalObj,[IsList,IsZmodnZMat],0,
function(a,b)
local fam;
  a:=List(a,r->List(r,Int));
  return a=b![1];
end);

InstallOtherMethod(\*,"ZmodnZVec*ZmodnZMat",IsElmsColls,[IsZmodnZVec,IsZmodnZMat],0,
function(a,b)
local fam;
  fam:=ElementsFamily(ElementsFamily(FamilyObj(b)));
  a:=a![1]*b![1];
  a:=a mod Characteristic(fam);
  return MakeZmodnZVec(fam,a);
end);

InstallOtherMethod(\^,"ZmodnZVec^ZmodnZMat",IsElmsColls,[IsZmodnZVec,IsZmodnZMat],0,
function(a,b)
local fam;
  fam:=ElementsFamily(ElementsFamily(FamilyObj(b)));
  a:=a![1]^b![1];
  a:=a mod Characteristic(fam);
  return MakeZmodnZVec(fam,a);
end);

InstallMethod(\*,"Rat*ZmodnZVec",true,[IsRat,IsZmodnZVec],0,
function(a,b)
local fam;
  fam:=ElementsFamily(FamilyObj(b));
  b:=a*b![1];
  b:=b mod Characteristic(fam);
  return MakeZmodnZVec(fam,b);
end);
