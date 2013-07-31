# Amira Script
remove -all
remove physics.icol 0_pulse_field_800.raw m1_pulse_field_800.raw p1_pulse_field_800.raw p12_pulse_field_800.raw m12_pulse_field_800.raw EmptyPlane Vectors

# Create viewers
viewer setVertical 0

viewer 0 setBackgroundMode 0
viewer 0 setBackgroundColor 0.06 0.13 0.24
viewer 0 setBackgroundColor2 0.72 0.72 0.78
viewer 0 setTransparencyType 5
viewer 0 setAutoRedraw 0
viewer 0 show
mainWindow show

set hideNewModules 1
[ load ${AMIRA_ROOT}/data/colormaps/physics.icol ] setLabel physics.icol
physics.icol setIconPosition 0 0
physics.icol setNoRemoveAll 1
physics.icol fire
{physics.icol} setMinMax 0 0.5
physics.icol flags setValue 1
physics.icol shift setMinMax -1 1
physics.icol shift setButtons 0
physics.icol shift setIncrement 0.133333
physics.icol shift setValue 0
physics.icol shift setSubMinMax -1 1
physics.icol scale setMinMax 0 1
physics.icol scale setButtons 0
physics.icol scale setIncrement 0.1
physics.icol scale setValue 1
physics.icol scale setSubMinMax 0 1
physics.icol fire
physics.icol setViewerMask 16383

set hideNewModules 0
[ load /data/miha/Svetloba/sreda/0_pulse_field_800.raw ] setLabel 0_pulse_field_800.raw
0_pulse_field_800.raw setIconPosition 20 10
0_pulse_field_800.raw fire
0_pulse_field_800.raw setViewerMask 16383

set hideNewModules 0
[ load /data/miha/Svetloba/sreda/m1_pulse_field_800.raw ] setLabel m1_pulse_field_800.raw
m1_pulse_field_800.raw setIconPosition 20 40
m1_pulse_field_800.raw fire
m1_pulse_field_800.raw setViewerMask 16383

set hideNewModules 0
[ load /data/miha/Svetloba/sreda/p1_pulse_field_800.raw ] setLabel p1_pulse_field_800.raw
p1_pulse_field_800.raw setIconPosition 20 70
p1_pulse_field_800.raw fire
p1_pulse_field_800.raw setViewerMask 16383

set hideNewModules 0
[ load /data/miha/Svetloba/sreda/p12_pulse_field_800.raw ] setLabel p12_pulse_field_800.raw
p12_pulse_field_800.raw setIconPosition 19 100
p12_pulse_field_800.raw fire
p12_pulse_field_800.raw setViewerMask 16383

set hideNewModules 0
[ load /data/miha/Svetloba/sreda/m12_pulse_field_800.raw ] setLabel m12_pulse_field_800.raw
m12_pulse_field_800.raw setIconPosition 20 130
m12_pulse_field_800.raw fire
m12_pulse_field_800.raw setViewerMask 16383

set hideNewModules 0
create HxArbitraryCut {EmptyPlane}
EmptyPlane setIconPosition 509 10
EmptyPlane data connect m12_pulse_field_800.raw
EmptyPlane fire
EmptyPlane origin  setBoundingBox -1e+08 1e+08 -1e+08 1e+08 -1e+08 1e+08
EmptyPlane origin  setImmediate 0
EmptyPlane origin  setOrtho 0
EmptyPlane origin  showDragger 0
EmptyPlane origin  showPoints 0
EmptyPlane origin  setPointScale 1
EmptyPlane origin  showOptionButton 1
EmptyPlane origin  setNumPoints 1 1 1
EmptyPlane origin  setValue 0 71.5 71.5 357.051
EmptyPlane normal  setBoundingBox -1e+08 1e+08 -1e+08 1e+08 -1e+08 1e+08
EmptyPlane normal  setImmediate 0
EmptyPlane normal  setOrtho 0
EmptyPlane normal  showDragger 0
EmptyPlane normal  showPoints 0
EmptyPlane normal  setPointScale 1
EmptyPlane normal  showOptionButton 1
EmptyPlane normal  setNumPoints 1 1 1
EmptyPlane normal  setValue 0 0 0 1
EmptyPlane frameSettings setState item 0 1 item 2 1 color 3 1 0.5 0 
EmptyPlane options setValue 0 1
EmptyPlane options setToggleVisible 0 1
EmptyPlane options setValue 1 0
EmptyPlane options setToggleVisible 1 1
EmptyPlane options setValue 2 1
EmptyPlane options setToggleVisible 2 1
EmptyPlane options setValue 3 0
EmptyPlane options setToggleVisible 3 1
EmptyPlane translate setMinMax 0 100
EmptyPlane translate setButtons 1
EmptyPlane translate setIncrement 1
EmptyPlane translate setValue 68.8
EmptyPlane translate setSubMinMax 0 100
EmptyPlane orientation untouch
EmptyPlane setMinPlanePoint -1e+15 -1e+15 -1e+15
EmptyPlane setMaxPlanePoint 1e+15 1e+15 1e+15
EmptyPlane setTranslateRange 101.000000
EmptyPlane fire
EmptyPlane setViewerMask 16383
EmptyPlane setPickable 1

set hideNewModules 0
create HxVectors {Vectors}
Vectors setIconPosition 509 30
Vectors setLineWidth 2
Vectors setLogScale 0
Vectors data connect m12_pulse_field_800.raw
Vectors module connect EmptyPlane
Vectors colormap setDefaultColor 0.13 0.08 0.6
Vectors colormap setDefaultAlpha 0.500000
Vectors colormap setLocalRange 0
Vectors colormap enableAlpha 0
Vectors colormap enableAlphaToggle 1
Vectors colormap connect physics.icol
Vectors fire
Vectors resolution setMinMax 0 -3.40282346638529e+38 3.40282346638529e+38
Vectors resolution setValue 0 100
Vectors resolution setMinMax 1 -3.40282346638529e+38 3.40282346638529e+38
Vectors resolution setValue 1 100
Vectors scale setMinMax 0 5
Vectors scale setButtons 0
Vectors scale setIncrement 0.333333
Vectors scale setValue 0.16
Vectors scale setSubMinMax 0 5
Vectors scaleZ setMinMax 0 1
Vectors scaleZ setButtons 0
Vectors scaleZ setIncrement 0.1
Vectors scaleZ setValue 1
Vectors scaleZ setSubMinMax 0 1
Vectors mode setValue 0 0
Vectors mode setToggleVisible 0 1
Vectors mode setValue 1 0
Vectors mode setToggleVisible 1 1
Vectors mode setValue 2 1
Vectors mode setToggleVisible 2 1
Vectors mode setValue 3 0
Vectors mode setToggleVisible 3 1
Vectors mode setValue 4 0
Vectors mode setToggleVisible 4 1
Vectors colorMode setIndex 0 0
Vectors phase setMinMax 0 360
Vectors phase setButtons 0
Vectors phase setIncrement 24
Vectors phase setValue 0
Vectors phase setSubMinMax 0 360
Vectors fire
Vectors setViewerMask 16383
Vectors setShadowStyle 0
Vectors setPickable 1

set hideNewModules 0


viewer 0 setCameraOrientation 1 0 0 3.14159
viewer 0 setCameraPosition 71.5 71.5 255.929
viewer 0 setCameraFocalDistance 101.117
viewer 0 setCameraNearDistance 100.682
viewer 0 setCameraFarDistance 101.709
viewer 0 setCameraType orthographic
viewer 0 setCameraHeight 147.656
viewer 0 setAutoRedraw 1
viewer 0 redraw

