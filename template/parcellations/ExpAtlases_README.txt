This file describes the experimental atlas files 

expFM800AtlasFile.mat
expFM400AtlasFile.mat
expFM200AtlasFile.mat
expConte69AtlasFile.mat

which contain various different parcellations.

RANDOM PARCELLATIONS
================================================
The files 
expFM800AtlasFile.mat
expFM400AtlasFile.mat
expFM200AtlasFile.mat

contain random parcellation of the 164K conte69 template from Caret 
in 800,400,200 parcels using the Fast marching method from the 
Fast Marching toolbox 
(https://www.ceremade.dauphine.fr/~peyre/teaching/manifold/tp3.html)

The above numbers refer to the entire cortical sheet. 
The number per hemisphere is half of these numbers.
The parcels have been created on the left hemisphere 
and then the matched voxels were used in the right hemisphere
as in the Conte69 atlas the left and right hemisphere voxels 
are supposed to be more or less matched. Thus these atlases are symmetric.

In each of these files you ll find a matlab variable
'atlasinfo' which contains a single structure.



atlasinfo{1}

ans = 

    parent: 'experimental'
      name: 'FastMarching800'
      mnem: 'FM800'
         L: [1x1 struct]
         R: [1x1 struct]

The field 'parent' mimics the Conte69 atlases, 
as there each atlas belongs to a parent family.
The 'name' is a name I have come up with and 
'mnem' is a mnemonic I have come up with.

The fields 'L','R' contain the actual atlas information PER HEMISPHERE.

 atlasinfo{1}.L

ans = 

           data: [163842x1 double]
          label: {1x400 cell}
            key: [1x400 double]
            rgb: [400x3 double]
        centpos: [400x3 double]
     parcelarea: [400x1 double]
    parcelNvoxs: [400x1 double]


Fields 
data:  contains the actual parcel index. This index is PER HEMISPHERE. 
        So for the 800 parcel random parcellation the values in this field
        take values between 1 and 400 in the Left hemisphere and 1 and 400 
        in the right hemisphere as well.
label: label of the parcels, They are labeled 'L_1','L_2'... for left and 
        'R_1','R_2'... for the right hemisphere.
key: This just contains the key values for each parcel. Here jsut between 
        1 and 400.
rgb: This contains a color code for each parcel. 
centpos: This contains the position of the central voxel of the parcel. 
         This is the parcel based on which the Voronoi distances were computed.
parcelarea: The area of each parcel.
parcelNvoxs: The number of voxels in each parcel.




CONTE69 atlases
==============================================
The file
expConte69AtlasFile.mat

contains 6 of the parcellations of the Conte69 atlas.

In this file you ll find a matlab variable
'atlasinfo' which contains 6 structures, one for each parcellation scheme.

These parcellations scheme have the following code names in the Conte69 atlas:



 atlasinfo{1}.name: 'CaretCompParc'  : Caret Composite Parcellation. Contains Visuotopic areas, some brodmann especially sensorimotor and orbitofrontal areas.
 atlasinfo{2}.name: 'CaretBrodmann'  : Brodmann areas
 atlasinfo{3}.name: 'CaretRSN7'      : 7 Resting state network parcels
 atlasinfo{4}.name: 'CaretRSN17'     : 17 Resting state network parcels
 atlasinfo{5}.name: 'CaretComCons'   : Caret Common Consensus. Parcels of the brain where researchers have assigned some kind of functionality.
 atlasinfo{6}.name: 'CaretComConsFilled'  : Same as above but with the unassigned cortex assigned.


Similarly to the random parcellations these atlases have a 'L' and a 'R' fiels which contain:

atlasinfo{1}.L

ans = 

     data: [163842x1 int32]
    label: {1x55 cell}
      key: [1x55 double]
     rgba: [55x4 double]


'data',and 'key' are similar as for the random parcellations.
the field 'name' is similat but does NOT have the prefix 'L_' and 'R_' so for example 
area V1 is called 'V1' in both hemispheres.
Here instad of the field 'rgb' there is the field 'rgba' which has an 
additional column which is supposed to represent oppacity.The reason 
this is here is because these colors are extracted from Caret files 
where this additional column exists.

Again all the atlas information is defined PER HEMISPHERE.













