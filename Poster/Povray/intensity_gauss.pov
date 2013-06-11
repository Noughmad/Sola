#include "colors.inc"    // The include files contain
#include "stones.inc"    // pre-defined scene elements

background {
  color rgb <0.75, 0.85, 0.0>
}

camera {
  location <2, 2, 0.5>
  look_at  <0.5, 0.5, 0.5>
}

#declare F = function { 
  pattern {
    density_file df3 "intensity_0.df3"
    interpolate 1
  }
}

isosurface {
  function {
    0.0002 - F(x,y,z/2)
  }
  contained_by { box { -0, 5 } }
  texture {
    pigment { color Green }
  }
}

light_source {
  <3, 0, 5> color White
}