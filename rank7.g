Read("roots.g");

for file in [
"E7.g", 
"g86a.g",  
"g86b.g", 
"g86c.g", 
"g86d.g", 
"g86e.g", 
"g86f.g", 
"g86g.g", 
"g86h.g"] do
  roots(Concatenation("rank7/", file));
  Print("--\n");
od;
