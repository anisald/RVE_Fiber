# RVE_Fiber
the random model
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random
from array import *
import math
import numpy
import os, glob      
import shutil   
from abaqus import * 
from abaqusConstants import *
import multiprocessing
import ctypes
import mesh

session.journalOptions.setValues(replayGeometry=COORDINATE)
mdb = Mdb()


#Input Parameters ---------------------------------------------------------------------------------------------------------------------------------------
field_inputs = (('Enter your volume fraction:','0.2'),('Enter sphere radius:','0.45'),('Enter the RVE side length:','7'),('Enter the the fiber radius:','0.45'),('Enter spring dash:','200'))
volumefraction,sphere_radius,RVE_length,fiber_radius,spring_dash = getInputs(fields=field_inputs,label='RVE Properties',dialogTitle = 'Create RVE')
volumefraction = eval(volumefraction)
sphere_radius = eval(sphere_radius)
RVE_length = eval(RVE_length)
fiber_radius= eval(fiber_radius)
spring_dash=eval(spring_dash)

#--------------------------------------------------

RVE_volume = RVE_length**3
fibers_volume = RVE_volume*volumefraction
fiber_volume = pi*(fiber_radius)**2*RVE_length
fibers_nu = int(math.floor(fibers_volume/fiber_volume))
thickness = fiber_radius/10
print(fiber_volume)
print(fibers_volume)
print(fibers_nu)

#Random Positions-----------------------------------------------------------------------------------------------------------------------------------------
def Random_positions (fibers_nu,thickness,RVE_length):
    dis=numpy.zeros(1000)
    num_incl = 0
    x_coordinate = []
    y_coordinate = []
    while (num_incl < fibers_nu):
            random_x=random.uniform((sphere_radius+thickness), RVE_length-(sphere_radius+4*thickness))
            random_y=random.uniform((sphere_radius+thickness), RVE_length-(sphere_radius+4*thickness))
            
            isPointIntersecting = False
            for j in range (0,len(x_coordinate)):
        
        
                dis[j]=sqrt((random_x-x_coordinate[j])**2+(random_y-y_coordinate[j])**2)

                    
                if dis[j] < (2.1*sphere_radius):

                    isPointIntersecting = True
                    break
                    
            if random_x >(RVE_length-sphere_radius) and random_x <(RVE_length+sphere_radius):
                isPointIntersecting = True
                break
            if random_y >(RVE_length-sphere_radius) and random_y <(RVE_length+sphere_radius):
                isPointIntersecting = True
                break 
                
            if (isPointIntersecting == False):
                x_coordinate.append(random_x)
                y_coordinate.append(random_y)
                num_incl = num_incl + 1    
    return  x_coordinate ,y_coordinate
    
x_coordinate ,y_coordinate = Random_positions (fibers_nu,thickness,RVE_length) 
 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
mdb.models.changeKey('Model-1','FirstModel')
mymodel = mdb.models['FirstModel'] 

#Matrix_part-----------------------------------------------------------------------------------------------------------------------------------------------------
mysketch1 = mymodel.ConstrainedSketch('firstsketch',500)
mypart1 = mymodel.Part('Firstpart',dimensionality = THREE_D,type=DEFORMABLE_BODY)
mysketch1.rectangle(point1=(0.0, 0.0), point2=(RVE_length, RVE_length))
mypart1.BaseSolidExtrude(sketch=mysketch1, depth=RVE_length)

#Solid Extrude---------------------------------------------------------------------------------------------------------------------------------------------------
f1, e1 = mypart1.faces, mypart1.edges
t = mypart1.MakeSketchTransform(sketchPlane=f1[4], sketchUpEdge=e1[7], 
    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 
    RVE_length))
mysketch2 = mymodel.ConstrainedSketch(name='__profile__', 
    sheetSize=RVE_length, gridSpacing=4.33, transform=t)
g, v, d, c = mysketch2.geometry, mysketch2.vertices, mysketch2.dimensions, mysketch2.constraints
mysketch2.setPrimaryObject(option=SUPERIMPOSE)
mypart1.projectReferencesOntoSketch(sketch=mysketch2, filter=COPLANAR_EDGES)
for i in range (0,fibers_nu):
     mysketch2.CircleByCenterPerimeter(center=(x_coordinate[i], y_coordinate[i]), point1=(x_coordinate[i], y_coordinate[i]+fiber_radius))
     
f, e = mypart1.faces, mypart1.edges
mypart1.SolidExtrude(sketchPlane=f[4], sketchUpEdge=e[7], sketchPlaneSide=SIDE1, 
    sketchOrientation=RIGHT, sketch=mysketch2, depth=RVE_length, flipExtrudeDirection=ON, 
    keepInternalBoundaries=ON)
mysketch2.unsetPrimaryObject()
del mymodel.sketches['__profile__']

#Matrix_Property-----------------------------------------------------------------------------------------------------------------------------------------------
mymaterial1 = mymodel.Material(name='MATRI')
mymaterial1.Hyperelastic(
    materialType=ISOTROPIC, testData=OFF, type=OGDEN, 
    volumetricResponse=VOLUMETRIC_DATA, table=((87.4, 4.49, 0.0), ))
mymaterial1.Density(table=((10.0, ), ))    
i = fibers_nu + 1
mysection1 = mymodel.HomogeneousSolidSection(name='Section-%d' %(i), material='MATRI', thickness=None)
mycells1 = mypart1.cells
myregion1 = (mycells1[-1],)
mypart1.SectionAssignment(region=myregion1, sectionName='Section-%d' %(i), offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION) 
    
#Fiber_Property----------------------------------------------------------------------------------------------------------------------------------------
for i in range(0,fibers_nu):
    mymaterial2 = mymodel.Material(name='Fiber-%d' %(i))
    mymaterial2.Hyperelastic(
        materialType=ISOTROPIC, testData=OFF, type=OGDEN,
        volumetricResponse=VOLUMETRIC_DATA, table=((1130.3, 4.91 
        0.0), ))
    mymaterial2.Density(table=((10.0, ), ))    
    mysection2 = mymodel.HomogeneousSolidSection(name='Section-%d' %(i), material='Fiber-%d' %(i), thickness=None)
    myregion2 = (mycells1[i],)
    mypart1.SectionAssignment(region=myregion2, sectionName='Section-%d' %(i), offset=0.0, 
         offsetType=MIDDLE_SURFACE, offsetField='', 
         thicknessAssignment=FROM_SECTION)      
#Matrix_Assembly-------------------------------------------------------------------------------------------------------------------------------------------------
myassembly = mymodel.rootAssembly
myassembly.DatumCsysByDefault(CARTESIAN) 
myinstance = myassembly.Instance(name='Firstpart-1', part=mypart1, dependent=ON)    from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random
from array import *
import math
import numpy
import os, glob      
import shutil   
from abaqus import * 
from abaqusConstants import *
import multiprocessing
import ctypes
import mesh

session.journalOptions.setValues(replayGeometry=COORDINATE)
mdb = Mdb()



#Input Parameters ---------------------------------------------------------------------------------------------------------------------------------------
field_inputs = (('Enter your volume fraction:','0.2'),('Enter sphere radius:','0.45'),('Enter the RVE side length:','7'),('Enter the the fiber radius:','0.45'),('Enter spring dash:','200'))
volumefraction,sphere_radius,RVE_length,fiber_radius,spring_dash = getInputs(fields=field_inputs,label='RVE Properties',dialogTitle = 'Create RVE')
volumefraction = eval(volumefraction)
sphere_radius = eval(sphere_radius)
RVE_length = eval(RVE_length)
fiber_radius= eval(fiber_radius)
spring_dash=eval(spring_dash)

#--------------------------------------------------

RVE_volume = RVE_length**3
fibers_volume = RVE_volume*volumefraction
fiber_volume = pi*(fiber_radius)**2*RVE_length
fibers_nu = int(math.floor(fibers_volume/fiber_volume))
thickness = fiber_radius/10
print(fiber_volume)
print(fibers_volume)
print(fibers_nu)

#Random Positions-----------------------------------------------------------------------------------------------------------------------------------------
def Random_positions (fibers_nu,thickness,RVE_length):
    dis=numpy.zeros(1000)
    num_incl = 0
    x_coordinate = []
    y_coordinate = []
    while (num_incl < fibers_nu):
            random_x=random.uniform((sphere_radius+thickness), RVE_length-(sphere_radius+4*thickness))
            random_y=random.uniform((sphere_radius+thickness), RVE_length-(sphere_radius+4*thickness))
            
            isPointIntersecting = False
            for j in range (0,len(x_coordinate)):
        
        
                dis[j]=sqrt((random_x-x_coordinate[j])**2+(random_y-y_coordinate[j])**2)

                    
                if dis[j] < (2.1*sphere_radius):

                    isPointIntersecting = True
                    break
                    
            if random_x >(RVE_length-sphere_radius) and random_x <(RVE_length+sphere_radius):
                isPointIntersecting = True
                break
            if random_y >(RVE_length-sphere_radius) and random_y <(RVE_length+sphere_radius):
                isPointIntersecting = True
                break 
                
            if (isPointIntersecting == False):
                x_coordinate.append(random_x)
                y_coordinate.append(random_y)
                num_incl = num_incl + 1    
    return  x_coordinate ,y_coordinate
    
x_coordinate ,y_coordinate = Random_positions (fibers_nu,thickness,RVE_length) 
 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------
mdb.models.changeKey('Model-1','FirstModel')
mymodel = mdb.models['FirstModel'] 

#Matrix_part-----------------------------------------------------------------------------------------------------------------------------------------------------
mysketch1 = mymodel.ConstrainedSketch('firstsketch',500)
mypart1 = mymodel.Part('Firstpart',dimensionality = THREE_D,type=DEFORMABLE_BODY)
mysketch1.rectangle(point1=(0.0, 0.0), point2=(RVE_length, RVE_length))
mypart1.BaseSolidExtrude(sketch=mysketch1, depth=RVE_length)

#Solid Extrude---------------------------------------------------------------------------------------------------------------------------------------------------
f1, e1 = mypart1.faces, mypart1.edges
t = mypart1.MakeSketchTransform(sketchPlane=f1[4], sketchUpEdge=e1[7], 
    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 
    RVE_length))
mysketch2 = mymodel.ConstrainedSketch(name='__profile__', 
    sheetSize=RVE_length, gridSpacing=4.33, transform=t)
g, v, d, c = mysketch2.geometry, mysketch2.vertices, mysketch2.dimensions, mysketch2.constraints
mysketch2.setPrimaryObject(option=SUPERIMPOSE)
mypart1.projectReferencesOntoSketch(sketch=mysketch2, filter=COPLANAR_EDGES)
for i in range (0,fibers_nu):
     mysketch2.CircleByCenterPerimeter(center=(x_coordinate[i], y_coordinate[i]), point1=(x_coordinate[i], y_coordinate[i]+fiber_radius))
     
f, e = mypart1.faces, mypart1.edges
mypart1.SolidExtrude(sketchPlane=f[4], sketchUpEdge=e[7], sketchPlaneSide=SIDE1, 
    sketchOrientation=RIGHT, sketch=mysketch2, depth=RVE_length, flipExtrudeDirection=ON, 
    keepInternalBoundaries=ON)
mysketch2.unsetPrimaryObject()
del mymodel.sketches['__profile__']
