#include "colors.inc"    // The include files contain
#include "stones.inc"    // pre-defined scene elements

background {
  color rgb <0.75, 0.85, 0.0>
}

camera {
  location <1.5, 1, 1.2>
  look_at  <0.5, 0.5, 0.5>
}

#declare F = function { 
  pattern {
    // Include all these (or only some) in the poster
    // density_file df3 "intensity_1000.df3"
    // density_file df3 "intensity_0.df3"
    // density_file df3 "intensity_2200.df3"
    // density_file df3 "intensity_500.df3"
    density_file df3 "intensity_1400.df3"
    interpolate 1
  }
}

isosurface {
  function {
    0.05 - F(x,y,z)
  }
  contained_by { box { 0.1, 0.9 } }
  texture {
    pigment { color Green }
    finish { phong 0.3 }
  }
}

isosurface {
  function {
    (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) - 0.16
  }
  contained_by { box {0, 1 } }
  texture {
    pigment { color Red transmit 0.7 }
  }
}

isosurface {
  function {
    (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) - 0.0001
  }
  contained_by { box {0, 1 } }
  texture {
    pigment { color Black transmit 0.2 }
    finish { phong 0.7 }
  }
}

light_source {
  <3, 0, 2> color shadowless
}

light_source {
  <3, 3, 3> color White
  area_light <1, 1, 0>, <0, 0, 1>, 5, 10
  adaptive 1
}