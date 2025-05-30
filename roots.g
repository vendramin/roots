###
### This script computes the root systems associated to a diagonal braiding
### Written by I. Angiono, L. Vendramin
### 
### This is Algorithm 3.1 
### 
### How to use this script? 
### 
### To compute a root system of rank two:
### gap> roots("rank2/ufo12b.g");
### File: rank2/ufo12b.g
### The braiding is:
### [ [  -E(7)^5,  -E(7)^3 ],
###   [        1,       -1 ] ]
### Longest word:
### [ 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2 ]
### Positive roots (12):
### [ [ 1, 0 ], [ 5, 1 ], [ 4, 1 ], [ 7, 2 ], [ 3, 1 ], [ 8, 3 ], [ 5, 2 ], [ 7, 3 ], [ 2, 1 ], [ 3, 2 ], [ 1, 1 ],
###   [ 0, 1 ] ]
### Cartan roots (6):
### [ [ 1, 0 ], [ 4, 1 ], [ 3, 1 ], [ 5, 2 ], [ 2, 1 ], [ 1, 1 ] ]


### Bounds (global variables)
### These bounds are useful for computing truncated infinite root systems

### This "M" is the bound for the exponents of q_ii
M := 7;   

### This "N" is the bound for the length of the longest word. 
### Use infinity for finite root systems.
N := 250;  

### Global variable
x := Indeterminate(Rationals, "x");
q := [];

### This function returns the new matrix (q_ij) with q_ij=1 for all i<j
convert := function(m)
  local new, i, j;

  new := NullMat(Size(m), Size(m));

  for i in [1..Size(m)] do
    for j in [1..Size(m)] do
      if i=j then
        new[i][j] := m[i][j];
      elif i>j then
        new[i][j] := m[i][j]*m[j][i];
      else
        new[i][j] := 1;
      fi;
    od;
  od;
  return new;
end;

### When all the entries computed are \geq -M-1, this 
### function computes the entries of the generalized 
### Cartan matrix associated to the braiding <q>. 
### Otherwise, it returns a matrix with an entry 
### equal to -M-1
a_ij := function(q)
  local a, n, i, j, m;

  n := Size(q);
  a := NullMat(n,n)-M-1;

  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        a[i][j] := 2;
      else
        ### a_ij <= M, M=7 works for finite root systems
        for m in [0..M] do
          if IsOne(q[i][i]^(m+1)) or IsOne(q[i][i]^m*q[i][j]*q[j][i]) then
            a[i][j] := -m;
            break;
          fi;
        od;
      fi;
    od;
  od;
  return a;
end;

### This functions returns the matrix of the reflection s_i
s := function(i, a)
  local n, m;
  n := Size(a);
  m := IdentityMat(n, n);
  m[i] := m[i]-a[i];
  return m;
end;

### This function computes the braiding after applying the reflection s_k
new_q := function(k, q, a)
  local n, m, i, j;
  n := Size(q);
  m := NullMat(n, n);
  for i in [1..n] do
    for j in [1..n] do
      m[i][j] := q[i][j]*q[i][k]^(-a[k][j])*q[k][j]^(-a[k][i])*q[k][k]^(a[k][j]*a[k][i]);
    od;
  od;
  return m;
end;

### This function returns true if q_ii^a_ij=q_ijq_ji
is_cartan := function(q, a, i)
  return ForAll([1..Size(q)], j->q[i][i]^(a[i][j]) = q[i][j]*q[j][i]);
end;

### Compute the roots
roots := function(file)
  local n, i, w, l, pos_roots, cartan_roots, done, a, r;

  Read(file);

  LogTo();
  LogTo(Concatenation(file, ".log"));
  
  Print("File: ", file, "\n");
  Print("The braiding is:\n");
  Display(q);
  
  n := Size(q);
  
  i := 1;
  w := IdentityMat(n, n);
  l := [1];
  
  pos_roots := [TransposedMat(w)[i]];
  cartan_roots := [];
  done := true;
  
  while done do
  
    a := a_ij(q);
    r := s(i, a);

    if ForAny(Flat(a), x->x <= -M) then
      Print("The root system is infinite.\c\n");
      return fail;
    fi;
  
    ### Is i a Cartan root?
    if is_cartan(q, a, i) then
      Add(cartan_roots, TransposedMat(w)[i]);
    fi;
  
    w := w*r;
    q := new_q(i, q, a);
  
    # It suffices to check if one entry of <w> is positive
    i := First(Concatenation([1..i-1],[i+1..n]), x->ForAll(TransposedMat(w)[x], y->y>=0));

    if Size(l) > Maximum([N, Size(q)^2]) then
      Print("The root system is infinite.\c\n");
      return fail;
    fi;  
   
    if i = fail then
      done := false;
    else 
      Add(pos_roots, TransposedMat(w)[i]);
      Add(l, i);
    fi;
    
  od;
  
  Print("Longest word:\n", l, "\n");
  Print("Positive roots (", Size(pos_roots), "):\n", pos_roots, "\n");
  Print("Cartan roots (", Size(cartan_roots), "):\n", cartan_roots, "\n");

  return pos_roots;

end;
