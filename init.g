#############################################################################
##
##  init.g                matgrp package                  Alexander Hulpke
##
##  Copyright 2014 by the author.
##  Licensed under GPL v2 or v3
##
#############################################################################

#I introducing globally the NC versions of PreImages...  
if not IsBound( PreImagesSetNC ) then 
    BindGlobal( "PreImagesSetNC", PreImagesSet ); 
fi; 
if not IsBound( PreImagesRepresentativeNC ) then 
    BindGlobal( "PreImagesRepresentativeNC", PreImagesRepresentative ); 
fi; 

#############################################################################

ReadPackage("matgrp","lib/recograt.gd");
