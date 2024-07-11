Read("roots.g");

for file in [
"brj25a.g", 
"brj25b.g",
"standard_b.g",
"standard_b_part.g",
"standard_g2a.g", 
"standard_g2b.g",
"standard_g2c.g",
"super_a.g",
"super_ap.g",
"super_b.g",
"ufo10a.g", 
"ufo10b.g",
"ufo10c.g",
"ufo10d.g",
"ufo11a.g",
"ufo11b.g",
"ufo11c.g",
"ufo11d.g",
"ufo12a.g",
"ufo12b.g",
"ufo7a.g", 
"ufo7b.g",
"ufo7c.g",
"ufo7d.g",
"ufo7e.g",
"ufo8a.g",
"ufo8b.g",
"ufo8c.g",
"ufo9a.g",
"ufo9b.g",
"ufo9c.g",
"ufo9d.g"
] do
  roots(Concatenation("rank2/", file));
  Print("--\n");
od;
