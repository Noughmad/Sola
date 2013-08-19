# Amira Script
remove -all
remove physics.icol p1_cont_field_2000.raw EmptyPlane Vectors EmptyPlane2 Vectors2 EmptyPlane3 Vectors3 EmptyPlane4 Vectors4 GlobalAxis

# Create viewers
viewer setVertical 0

viewer 0 setBackgroundMode 1
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
{physics.icol} setMinMax 0 1
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
[ load "/home/miha/Dropbox/Podatki - slow/p1_cont_field_2000.raw" ] setLabel p1_cont_field_2000.raw
p1_cont_field_2000.raw setIconPosition 20 10
p1_cont_field_2000.raw fire
p1_cont_field_2000.raw setViewerMask 16383

set hideNewModules 0
create HxArbitraryCut {EmptyPlane}
EmptyPlane setIconPosition 192 162
EmptyPlane data connect p1_cont_field_2000.raw
EmptyPlane fire
EmptyPlane origin  setBoundingBox -1e+08 1e+08 -1e+08 1e+08 -1e+08 1e+08
EmptyPlane origin  setImmediate 0
EmptyPlane origin  setOrtho 0
EmptyPlane origin  showDragger 0
EmptyPlane origin  showPoints 0
EmptyPlane origin  setPointScale 1
EmptyPlane origin  showOptionButton 1
EmptyPlane origin  setNumPoints 1 1 1
EmptyPlane origin  setValue 0 71.5 71.5 17.1429
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
EmptyPlane translate setValue 1.65289
EmptyPlane translate setSubMinMax 0 100
EmptyPlane orientation untouch
EmptyPlane setMinPlanePoint -1e+15 -1e+15 -1e+15
EmptyPlane setMaxPlanePoint 1e+15 1e+15 1e+15
EmptyPlane setTranslateRange 101.000000
EmptyPlane fire
EmptyPlane setViewerMask 16382
EmptyPlane setPickable 1

set hideNewModules 0
create HxVectors {Vectors}
Vectors setIconPosition 192 182
Vectors setLineWidth 2
Vectors setLogScale 0
Vectors data connect p1_cont_field_2000.raw
Vectors module connect EmptyPlane
Vectors colormap setDefaultColor 0.13 0.08 0.6
Vectors colormap setDefaultAlpha 0.500000
Vectors colormap setLocalRange 0
Vectors colormap enableAlpha 0
Vectors colormap enableAlphaToggle 1
Vectors colormap connect physics.icol
Vectors fire
Vectors resolution setMinMax 0 -3.40282346638529e+38 3.40282346638529e+38
Vectors resolution setValue 0 120
Vectors resolution setMinMax 1 -3.40282346638529e+38 3.40282346638529e+38
Vectors resolution setValue 1 120
Vectors scale setMinMax 0 5
Vectors scale setButtons 0
Vectors scale setIncrement 0.333333
Vectors scale setValue 0.3
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
create HxArbitraryCut {EmptyPlane2}
EmptyPlane2 setIconPosition 377 201
EmptyPlane2 data connect p1_cont_field_2000.raw
EmptyPlane2 fire
EmptyPlane2 origin  setBoundingBox -1e+08 1e+08 -1e+08 1e+08 -1e+08 1e+08
EmptyPlane2 origin  setImmediate 0
EmptyPlane2 origin  setOrtho 0
EmptyPlane2 origin  showDragger 0
EmptyPlane2 origin  showPoints 0
EmptyPlane2 origin  setPointScale 1
EmptyPlane2 origin  showOptionButton 1
EmptyPlane2 origin  setNumPoints 1 1 1
EmptyPlane2 origin  setValue 0 71.5 71.5 231.905
EmptyPlane2 normal  setBoundingBox -1e+08 1e+08 -1e+08 1e+08 -1e+08 1e+08
EmptyPlane2 normal  setImmediate 0
EmptyPlane2 normal  setOrtho 0
EmptyPlane2 normal  showDragger 0
EmptyPlane2 normal  showPoints 0
EmptyPlane2 normal  setPointScale 1
EmptyPlane2 normal  showOptionButton 1
EmptyPlane2 normal  setNumPoints 1 1 1
EmptyPlane2 normal  setValue 0 0 0 1
EmptyPlane2 frameSettings setState item 0 1 item 2 1 color 3 1 0.5 0 
EmptyPlane2 options setValue 0 0
EmptyPlane2 options setToggleVisible 0 1
EmptyPlane2 options setValue 1 0
EmptyPlane2 options setToggleVisible 1 1
EmptyPlane2 options setValue 2 1
EmptyPlane2 options setToggleVisible 2 1
EmptyPlane2 options setValue 3 0
EmptyPlane2 options setToggleVisible 3 1
EmptyPlane2 translate setMinMax 0 100
EmptyPlane2 translate setButtons 1
EmptyPlane2 translate setIncrement 1
EmptyPlane2 translate setValue 22.4876
EmptyPlane2 translate setSubMinMax 0 100
EmptyPlane2 orientation untouch
EmptyPlane2 setMinPlanePoint -1e+15 -1e+15 -1e+15
EmptyPlane2 setMaxPlanePoint 1e+15 1e+15 1e+15
EmptyPlane2 setTranslateRange 101.000000
EmptyPlane2 fire
EmptyPlane2 setViewerMask 16382
EmptyPlane2 select
EmptyPlane2 setPickable 1

set hideNewModules 0
create HxVectors {Vectors2}
Vectors2 setIconPosition 377 221
Vectors2 setLineWidth 3
Vectors2 setLogScale 0
Vectors2 data connect p1_cont_field_2000.raw
Vectors2 module connect EmptyPlane2
Vectors2 colormap setDefaultColor 0.13 0.08 0.6
Vectors2 colormap setDefaultAlpha 0.500000
Vectors2 colormap setLocalRange 0
Vectors2 colormap enableAlpha 0
Vectors2 colormap enableAlphaToggle 1
Vectors2 colormap connect physics.icol
Vectors2 fire
Vectors2 resolution setMinMax 0 -3.40282346638529e+38 3.40282346638529e+38
Vectors2 resolution setValue 0 150
Vectors2 resolution setMinMax 1 -3.40282346638529e+38 3.40282346638529e+38
Vectors2 resolution setValue 1 150
Vectors2 scale setMinMax 0 5
Vectors2 scale setButtons 0
Vectors2 scale setIncrement 0.333333
Vectors2 scale setValue 0.68323
Vectors2 scale setSubMinMax 0 5
Vectors2 scaleZ setMinMax 0 1
Vectors2 scaleZ setButtons 0
Vectors2 scaleZ setIncrement 0.1
Vectors2 scaleZ setValue 1
Vectors2 scaleZ setSubMinMax 0 1
Vectors2 mode setValue 0 0
Vectors2 mode setToggleVisible 0 1
Vectors2 mode setValue 1 0
Vectors2 mode setToggleVisible 1 1
Vectors2 mode setValue 2 0
Vectors2 mode setToggleVisible 2 1
Vectors2 mode setValue 3 0
Vectors2 mode setToggleVisible 3 1
Vectors2 mode setValue 4 0
Vectors2 mode setToggleVisible 4 1
Vectors2 colorMode setIndex 0 0
Vectors2 phase setMinMax 0 360
Vectors2 phase setButtons 0
Vectors2 phase setIncrement 24
Vectors2 phase setValue 0
Vectors2 phase setSubMinMax 0 360
Vectors2 fire
Vectors2 setViewerMask 16383
Vectors2 select
Vectors2 setShadowStyle 0
Vectors2 setPickable 1

set hideNewModules 0
create HxArbitraryCut {EmptyPlane3}
EmptyPlane3 setIconPosition 164 75
EmptyPlane3 data connect p1_cont_field_2000.raw
EmptyPlane3 fire
EmptyPlane3 origin  setBoundingBox -1e+08 1e+08 -1e+08 1e+08 -1e+08 1e+08
EmptyPlane3 origin  setImmediate 0
EmptyPlane3 origin  setOrtho 0
EmptyPlane3 origin  showDragger 0
EmptyPlane3 origin  showPoints 0
EmptyPlane3 origin  setPointScale 1
EmptyPlane3 origin  showOptionButton 1
EmptyPlane3 origin  setNumPoints 1 1 1
EmptyPlane3 origin  setValue 0 71.5 71.5 570.873
EmptyPlane3 normal  setBoundingBox -1e+08 1e+08 -1e+08 1e+08 -1e+08 1e+08
EmptyPlane3 normal  setImmediate 0
EmptyPlane3 normal  setOrtho 0
EmptyPlane3 normal  showDragger 0
EmptyPlane3 normal  showPoints 0
EmptyPlane3 normal  setPointScale 1
EmptyPlane3 normal  showOptionButton 1
EmptyPlane3 normal  setNumPoints 1 1 1
EmptyPlane3 normal  setValue 0 0 0 1
EmptyPlane3 frameSettings setState item 0 1 item 2 1 color 3 1 0.5 0 
EmptyPlane3 options setValue 0 0
EmptyPlane3 options setToggleVisible 0 1
EmptyPlane3 options setValue 1 0
EmptyPlane3 options setToggleVisible 1 1
EmptyPlane3 options setValue 2 1
EmptyPlane3 options setToggleVisible 2 1
EmptyPlane3 options setValue 3 0
EmptyPlane3 options setToggleVisible 3 1
EmptyPlane3 translate setMinMax 0 100
EmptyPlane3 translate setButtons 1
EmptyPlane3 translate setIncrement 1
EmptyPlane3 translate setValue 55.3719
EmptyPlane3 translate setSubMinMax 0 100
EmptyPlane3 orientation untouch
EmptyPlane3 setMinPlanePoint -1e+15 -1e+15 -1e+15
EmptyPlane3 setMaxPlanePoint 1e+15 1e+15 1e+15
EmptyPlane3 setTranslateRange 101.000000
EmptyPlane3 fire
EmptyPlane3 setViewerMask 16382
EmptyPlane3 setPickable 1

set hideNewModules 0
create HxVectors {Vectors3}
Vectors3 setIconPosition 164 95
Vectors3 setLineWidth 3
Vectors3 setLogScale 0
Vectors3 data connect p1_cont_field_2000.raw
Vectors3 module connect EmptyPlane3
Vectors3 colormap setDefaultColor 0.13 0.08 0.6
Vectors3 colormap setDefaultAlpha 0.500000
Vectors3 colormap setLocalRange 0
Vectors3 colormap enableAlpha 0
Vectors3 colormap enableAlphaToggle 1
Vectors3 colormap connect physics.icol
Vectors3 fire
Vectors3 resolution setMinMax 0 -3.40282346638529e+38 3.40282346638529e+38
Vectors3 resolution setValue 0 150
Vectors3 resolution setMinMax 1 -3.40282346638529e+38 3.40282346638529e+38
Vectors3 resolution setValue 1 150
Vectors3 scale setMinMax 0 5
Vectors3 scale setButtons 0
Vectors3 scale setIncrement 0.333333
Vectors3 scale setValue 0.3
Vectors3 scale setSubMinMax 0 5
Vectors3 scaleZ setMinMax 0 1
Vectors3 scaleZ setButtons 0
Vectors3 scaleZ setIncrement 0.1
Vectors3 scaleZ setValue 1
Vectors3 scaleZ setSubMinMax 0 1
Vectors3 mode setValue 0 0
Vectors3 mode setToggleVisible 0 1
Vectors3 mode setValue 1 0
Vectors3 mode setToggleVisible 1 1
Vectors3 mode setValue 2 1
Vectors3 mode setToggleVisible 2 1
Vectors3 mode setValue 3 0
Vectors3 mode setToggleVisible 3 1
Vectors3 mode setValue 4 0
Vectors3 mode setToggleVisible 4 1
Vectors3 colorMode setIndex 0 0
Vectors3 phase setMinMax 0 360
Vectors3 phase setButtons 0
Vectors3 phase setIncrement 24
Vectors3 phase setValue 0
Vectors3 phase setSubMinMax 0 360
Vectors3 fire
Vectors3 setViewerMask 16383
Vectors3 setShadowStyle 0
Vectors3 setPickable 1

set hideNewModules 0
create HxArbitraryCut {EmptyPlane4}
EmptyPlane4 setIconPosition 379 75
EmptyPlane4 data connect p1_cont_field_2000.raw
EmptyPlane4 fire
EmptyPlane4 origin  setBoundingBox -1e+08 1e+08 -1e+08 1e+08 -1e+08 1e+08
EmptyPlane4 origin  setImmediate 0
EmptyPlane4 origin  setOrtho 0
EmptyPlane4 origin  showDragger 0
EmptyPlane4 origin  showPoints 0
EmptyPlane4 origin  setPointScale 1
EmptyPlane4 origin  showOptionButton 1
EmptyPlane4 origin  setNumPoints 1 1 1
EmptyPlane4 origin  setValue 0 71.5 71.5 434.57
EmptyPlane4 normal  setBoundingBox -1e+08 1e+08 -1e+08 1e+08 -1e+08 1e+08
EmptyPlane4 normal  setImmediate 0
EmptyPlane4 normal  setOrtho 0
EmptyPlane4 normal  showDragger 0
EmptyPlane4 normal  showPoints 0
EmptyPlane4 normal  setPointScale 1
EmptyPlane4 normal  showOptionButton 1
EmptyPlane4 normal  setNumPoints 1 1 1
EmptyPlane4 normal  setValue 0 0 0 1
EmptyPlane4 frameSettings setState item 0 1 item 2 1 color 3 1 0.5 0 
EmptyPlane4 options setValue 0 0
EmptyPlane4 options setToggleVisible 0 1
EmptyPlane4 options setValue 1 0
EmptyPlane4 options setToggleVisible 1 1
EmptyPlane4 options setValue 2 1
EmptyPlane4 options setToggleVisible 2 1
EmptyPlane4 options setValue 3 0
EmptyPlane4 options setToggleVisible 3 1
EmptyPlane4 translate setMinMax 0 100
EmptyPlane4 translate setButtons 1
EmptyPlane4 translate setIncrement 1
EmptyPlane4 translate setValue 42.1488
EmptyPlane4 translate setSubMinMax 0 100
EmptyPlane4 orientation untouch
EmptyPlane4 setMinPlanePoint -1e+15 -1e+15 -1e+15
EmptyPlane4 setMaxPlanePoint 1e+15 1e+15 1e+15
EmptyPlane4 setTranslateRange 101.000000
EmptyPlane4 fire
EmptyPlane4 setViewerMask 16382
EmptyPlane4 setPickable 1

set hideNewModules 0
create HxVectors {Vectors4}
Vectors4 setIconPosition 379 95
Vectors4 setLineWidth 3
Vectors4 setLogScale 0
Vectors4 data connect p1_cont_field_2000.raw
Vectors4 module connect EmptyPlane4
Vectors4 colormap setDefaultColor 0.13 0.08 0.6
Vectors4 colormap setDefaultAlpha 0.500000
Vectors4 colormap setLocalRange 0
Vectors4 colormap enableAlpha 0
Vectors4 colormap enableAlphaToggle 1
Vectors4 colormap connect physics.icol
Vectors4 fire
Vectors4 resolution setMinMax 0 -3.40282346638529e+38 3.40282346638529e+38
Vectors4 resolution setValue 0 150
Vectors4 resolution setMinMax 1 -3.40282346638529e+38 3.40282346638529e+38
Vectors4 resolution setValue 1 150
Vectors4 scale setMinMax 0 5
Vectors4 scale setButtons 0
Vectors4 scale setIncrement 0.333333
Vectors4 scale setValue 0.3
Vectors4 scale setSubMinMax 0 5
Vectors4 scaleZ setMinMax 0 1
Vectors4 scaleZ setButtons 0
Vectors4 scaleZ setIncrement 0.1
Vectors4 scaleZ setValue 1
Vectors4 scaleZ setSubMinMax 0 1
Vectors4 mode setValue 0 0
Vectors4 mode setToggleVisible 0 1
Vectors4 mode setValue 1 0
Vectors4 mode setToggleVisible 1 1
Vectors4 mode setValue 2 1
Vectors4 mode setToggleVisible 2 1
Vectors4 mode setValue 3 0
Vectors4 mode setToggleVisible 3 1
Vectors4 mode setValue 4 0
Vectors4 mode setToggleVisible 4 1
Vectors4 colorMode setIndex 0 0
Vectors4 phase setMinMax 0 360
Vectors4 phase setButtons 0
Vectors4 phase setIncrement 24
Vectors4 phase setValue 0
Vectors4 phase setSubMinMax 0 360
Vectors4 fire
Vectors4 setViewerMask 16383
Vectors4 setShadowStyle 0
Vectors4 setPickable 1

set hideNewModules 0
create HxAxis {GlobalAxis}
GlobalAxis setIconPosition 171 231
GlobalAxis fire
GlobalAxis axis setValue 0 1
GlobalAxis axis setToggleVisible 0 1
GlobalAxis axis setValue 1 1
GlobalAxis axis setToggleVisible 1 1
GlobalAxis axis setValue 2 1
GlobalAxis axis setToggleVisible 2 1
GlobalAxis options setValue 0 1
GlobalAxis options setToggleVisible 0 1
GlobalAxis options setValue 1 0
GlobalAxis options setToggleVisible 1 1
GlobalAxis options setValue 2 0
GlobalAxis options setToggleVisible 2 1
GlobalAxis thickness setMinMax 1 15
GlobalAxis thickness setButtons 0
GlobalAxis thickness setIncrement 0.933333
GlobalAxis thickness setValue 15
GlobalAxis thickness setSubMinMax 1 15
GlobalAxis color setColor 0 1 0 0
GlobalAxis color setColor 1 0 1 0
GlobalAxis color setColor 2 0 0 1
GlobalAxis color setColor 3 1 0.8 0.5
GlobalAxis axisNames setState text 0 x text 1 y text 2 z 
GlobalAxis font setState name: {Helvetica} size: 10 bold: 0 italic: 0 color: 0.8 0.8 0.8
GlobalAxis fire
GlobalAxis setBoundingBox 10 30 -50 -30 130 170
{GlobalAxis} setDelta 0 
{GlobalAxis} setLocalMode 0 
{GlobalAxis} setFlip 0 0 
{GlobalAxis} setFlip 1 0 
{GlobalAxis} setFlip 2 0 
GlobalAxis fire
GlobalAxis setViewerMask 16383
GlobalAxis setPickable 1

set hideNewModules 0


viewer 0 setCameraOrientation 0.955073 -0.0159039 -0.295943 3.00324
viewer 0 setCameraPosition -99.0904 7.82669 47.5772
viewer 0 setCameraFocalDistance 302.02
viewer 0 setCameraNearDistance 24.1403
viewer 0 setCameraFarDistance 580.471
viewer 0 setCameraType orthographic
viewer 0 setCameraHeight 226.722
viewer 0 setAutoRedraw 1
viewer 0 redraw

