###
### This script computes the Lyndon words and the hyperwords 
### associated to a set of positive roots
### Written by I. Angiono, L. Vendramin
### 
### These are Algorithms 3.2 and 3.3 
### 
### How to use this script? 

Read("roots.g");

h := function(root)
  return Sum(root);
end;

lyndon_from_positive_roots := function(pos_roots)
  local l, x, words, decompositions, i, j, k, m, p, q, tmp;

  l := List(pos_roots, x->[Position(pos_roots, x), x, h(x)]);
  words := List([1..Size(pos_roots)], x->0);
  Sort(l, function(v, w) return v[3]<=w[3]; end );

  decompositions := List([1..Size(pos_roots)], x->0);

  for x in l do
    i := Position(l, x);
    tmp := [];
    if x[3]=1 then
      j := Position(IdentityMat(Size(pos_roots[1])), x[2]);
      words[i] := [j];
    else
      for k in [1..x[1]-1] do
        m := Position(pos_roots, x[2]-pos_roots[k]);
        if not m = fail then
          p := Position(List(l, y->y[1]), k);
          q := Position(List(l, y->y[1]), m);
          Add(tmp, [words[p], words[q]]);
        fi;
      od;
      words[i] := Maximum(List(tmp, y->Concatenation(y[1], y[2])));
      decompositions[i] := Filtered(tmp, y->Concatenation(y[1],y[2])=words[i]);
    fi;
  od;

  return rec( words := words, decompositions := decompositions, ordering := List(l, x->x[1]) );
  
end;

lyndon := function(file)
  local pos_roots;
  pos_roots := roots(file);
  return lyndon_from_positive_roots(pos_roots);
end; 

hyperwords_from_positive_roots := function(pos_roots)
  local alpha, beta, gamma, lst, t, l, words, x, decompositions;

  l := lyndon_from_positive_roots(pos_roots);
  words := l!.words;
  decompositions := l!.decompositions;

  for x in words do
    if Size(x)=1 then
      Print("x_{\\alpha_{", x[1], "}}=x_{", x[1], "}\n");
    else
      lst := List(decompositions[Position(words, x)], y->Size(y[1]));
      t := Position(lst, Minimum(lst));
      alpha := pos_roots[l!.ordering[Position(words, x)]];
      beta := pos_roots[l!.ordering[Position(words, decompositions[Position(words, x)][t][1])]];
      gamma := pos_roots[l!.ordering[Position(words, decompositions[Position(words, x)][t][2])]];
      Print("x_{", alpha, "}=[x_{", beta, "},x_{", gamma, "}]_c\n");
    fi;
  od;
       
end;

hyperwords := function(file)
  local pos_roots;
  pos_roots := roots(file);
  hyperwords_from_positive_roots(pos_roots);
end; 


