#include "colors.inc"
#include "transforms.inc"

#declare T = texture {
  pigment { color Green }
  finish { phong 0.3 }
}

background {
  color rgb <0.5, 0.7, 0.0>
}

camera {
  location <3.5, 2, 6>
  look_at  <0.5, 0.4, 3>
  right x*3
  angle 50
  sky <-0.35, 1, 0>
}

#declare F0 = function { 
  pattern {
    density_file df3 "intensity_0.df3"
    interpolate 1
  }
}

#declare F1 = function { 
  pattern {
    density_file df3 "intensity_500.df3"
    interpolate 1
  }
}

#declare F9 = function { 
  pattern {
    density_file df3 "intensity_950.df3"
    interpolate 1
  }
}

#declare F4 = function { 
  pattern {
    density_file df3 "intensity_1400.df3"
    interpolate 1
  }
}

#declare F2 = function { 
  pattern {
    density_file df3 "intensity_2200.df3"
    interpolate 1
  }
}

#declare F3 = function { 
  pattern {
    density_file df3 "intensity_3000.df3"
    interpolate 1
  }
}

isosurface {
  function {
    0.05 - F0(x,y,z-0.3)
  }
  contained_by { box { 0.1, 0.9 } }
  texture { T }
}

isosurface {
  function {
    0.05 - F1(x,y,z-0.8)
  }
  contained_by { box { <0, 0, 0.5>, <1, 1, 1.5> } }
  texture { T }
}

isosurface {
  function {
    0.05 - F9(x,y,z-1.4)
  }
  contained_by { box { <0, 0, 1.5>, <1, 1, 2.5> } }
  texture { T }
}

isosurface {
  function {
    0.05 - F4(x,y,z-2.1)
  }
  contained_by { box { <0, 0, 2.2>, <1, 1, 3.5> } }
  texture { T }
}

isosurface {
  function {
    0.05 - F2(x,y,z-3)
  }
  contained_by { box { <0, 0, 2.8>, <1, 1, 4> } }
  texture { T }
}

isosurface {
  function {
    0.05 - F3(x,y,z-4)
  }
  contained_by { box { <0.25, 0.25, 3.5>, <0.75, 0.75, 5> } }
  texture { T }
}

isosurface {
  function {
    (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) - 0.12
  }
  contained_by { box {<0, 0, 0.019>, <1, 1, 4.781> } }
  texture {
    pigment { color Red transmit 0.8 }
  }
}

isosurface {
  function {
    (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) - 0.12
  }
  contained_by { box {<0, 0, 0>, <1, 1, 0.02> } }
  open
  texture {
    pigment { color Red transmit 0.2 }
  }
}

isosurface {
  function {
    (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) - 0.12
  }
  contained_by { box {<0, 0, 4.78>, <1, 1, 4.8> } }
  open
  texture {
    pigment { color Red transmit 0.2 }
  }
}

isosurface {
  function {
    (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) - 0.0001
  }
  contained_by { box {<0, 0, 0> <1, 1, 4.8> } }
  texture {
    pigment { color Black transmit 0.2 }
    finish { phong 0.7 }
  }
}

union {
  cylinder {
    <1, 0, 0>, <1, 0, 0.25>, 0.02
    texture {pigment {color Blue}}
  }
  cone {
    <1, 0, 0.25>, 0.04, <1, 0, 0.4>, 0
    texture {pigment {color Blue}}
  }
}

union {
  cylinder {
    <1, 0, 0>, <1, 0.25, 0>, 0.02
    texture {pigment {color Green}}
  }
  cone {
    <1, 0.25, 0>, 0.04, <1, 0.4>, 0
    texture {pigment {color Green}}
  }
}

union {
  cylinder {
    <1, 0, 0>, <0.75, 0, 0>, 0.02
    texture {pigment {color Red}}
  }
  cone {
    <0.75, 0, 0>, 0.04, <0.6, 0, 0>, 0
    texture {pigment {color Red}}
  }
}

#macro Axes(Texture)
  union {
    cylinder { <0,0.35,0>,<0,AxisLen,0>,0.05
               texture{Texture}
                       translate<0.1,0,0.1>}
             }
    cone{<0,AxisLen,0>,0.2,<0,AxisLen+0.7,0>,0
          texture{Dark_Texture}
         }
     } // end of union
#end

light_source {
  <3, 0, 2> color shadowless
}

light_source {
  <3, 3, 3> color White
  area_light <1, 1, 0>, <0, 0, 1>, 5, 10
  adaptive 1
}