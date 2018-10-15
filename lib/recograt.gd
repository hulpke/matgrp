#############################################################################
##
#W  recograt.gd                 matgrp package               Alexander Hulpke
##
##
#Y  Copyright (C)  2014-18, Alexander Hulpke
##
##  basic setup for matrix fitting free.
##

DeclareInfoClass("InfoFFMat");

#############################################################################
##
#R  IsPcgsMatGroupByStabChainRep   . . . . . . . . . . . .  pcgs of mat group
##
##  This is the representation for a pcgs of a matrix group which computes
##  exponents via a stabilizer chain. It may not be set for subsets (tails)
##  as this could lead to wrong exponents.
##
DeclareRepresentation( "IsPcgsMatGroupByStabChainRep",
    IsPcgsDefaultRep and IsFiniteOrdersPcgs, [ "stabilizerChain" ] );

#############################################################################
##
#R  IsModuloPcgsPermGroupRep  . . induced/modulo pcgs for layers of mat group
##
DeclareRepresentation( "IsSublayersPcgsMatGroupRep",
    IsPcgsMatGroupByStabChainRep, [ "stabilizerChain", "layers"] );

#############################################################################
##
#R  IsPcgsResidueMatGroupRep   . . . . . . . . . . . .  pcgs of mat group
##
##  This is the representation for a pcgs of a matrix group over a residue
##  class ring.
##
DeclareRepresentation( "IsPcgsResidueMatGroupRep",
    IsPcgsDefaultRep and IsFiniteOrdersPcgs, [ "decomp" ] );

