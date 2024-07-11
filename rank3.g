Read("roots.g");

for file in [
"G3a.g",  
"G3b.g", 
"G3c.g", 
"G3d.g", 
"affine.g",
"br3a.g", 
"br3b.g",
"g16a.g",
"g16b.g",
"g23a.g",
"g23b.g",
"g23c.g",
"g23d.g",
"ufo3a.g", 
"ufo3b.g",
"ufo3c.g",
"ufo3d.g",
"ufo3e.g",
"ufo4a.g",
"ufo4b.g",
"ufo4c.g",
"ufo4d.g",
"ufo4e.g",
"ufo4f.g",
"ufo4g.g",
"ufo4h.g",
"ufo4i.g"] do
  roots(Concatenation("rank3/", file));
  Print("--\n");
od;
