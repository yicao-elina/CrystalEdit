#version 3.7;
// #include "colors.inc" // Not strictly needed if we define colors manually, but good practice
// #include "finish.inc"

global_settings { assumed_gamma 1.0 }
background { color rgb <1,1,1> }

// ---------- Finishes (from example.pov) ----------
#declare ase3 = finish {ambient 0.4 brilliance 2 diffuse 0.6 metallic specular 1.0 roughness 0.001 reflection 0.0}

// ---------- File IO ----------
#declare N = 0;
#fopen F1 "atoms.xyz" read
#read (F1, N)

#declare S = "";
#declare X = 0.0;
#declare Y = 0.0;
#declare Z = 0.0;

#declare AtomsUnion = union {
    #for (i, 1, N)
        #read (F1, S, X, Y, Z)
        #if (Z < 0.0)
            #declare Col = <0.2, 0.4, 0.8>; // Blueish
        #else
            #declare Col = <0.8, 0.2, 0.2>; // Reddish
        #end
        
        sphere {
            <X, Y, Z>, 0.8
            texture {
                pigment { color rgb Col }
                finish { ase3 } // Use nicer finish
            }
        }
    #end
}
#fclose F1

// ---------- Scene Setup ----------
camera {
    location <0, -120, 60>
    look_at <0, 0, 0>
    angle 30
}

light_source {
    <200, -200, 200>
    color rgb <1,1,1>
    // area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3 // Optional soft shadows
    // adaptive 1 jitter
}

// ---------- Twin Plane ----------
#declare Z_TWIN = 0.0;
#declare PLANE_THICK = 0.3;

box {
    <-100, -100, Z_TWIN - PLANE_THICK/2>,
    < 100,  100, Z_TWIN + PLANE_THICK/2>
    texture {
        pigment { color rgbt <0.7, 0.7, 0.7, 0.65> }
        finish { ambient 0.2 diffuse 0.2 specular 0.5 }
    }
}

// ---------- Normal ----------
#declare n233 = vnormalize(<2,3,3>);
#declare ArrowLen = 18;

cylinder {
    <0,0,Z_TWIN>,
    <0,0,Z_TWIN> + ArrowLen * n233
    0.6
    texture { pigment { color rgb <0,0,0> } finish { specular 0.5 } }
}

cone {
    <0,0,Z_TWIN> + ArrowLen * n233, 2.0,
    <0,0,Z_TWIN> + (ArrowLen + 4) * n233, 0.0
    texture { pigment { color rgb <0,0,0> } finish { specular 0.5 } }
}

// ---------- Render ----------
object { AtomsUnion }
