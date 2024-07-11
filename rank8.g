Read("roots.g");

for file in ["E8.g"] do
  roots(Concatenation("rank8/", file));
  Print("--\n");
od;
