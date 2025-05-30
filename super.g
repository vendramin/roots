### This script computes the root systems associated to Lie superalgebras
### Written by I. Angiono, L. Vendramin
### 

### Bounds (global variables)
### These bounds are useful for computing truncated infinite root systems
b := [];
p := [];

### This "M" is the bound for the exponents of q_ii
M := 7;   


### This "N" is the bound for the length of the longest word. 
### Use infinity for finite root systems.
### TODO: This bound should depend on the rank (to be sure that the q is OK)

from_ptoZ := function(p, x)
  local n;
  n := Position([0..p-1]*Z(p)^0, x)-1;
  if n = 0 then 
    return 0;
  else
    return n-p;
  fi;
end;

normalization := function(b)
  local i, m;
  m := NullMat(Size(b), Size(b));
  for i in [1..Size(b)] do
    if not IsZero(b[i][i]) then 
      b[i] := 2*b[i]/b[i][i];
    else
      if not IsZero(b[i]) then
        b[i] := b[i]/First(b[i], x->not IsZero(x));
      fi;
    fi;
  od;
  return b;
end;

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

### This function computes the entries of the generalized Cartan matrix associated to the braiding <q>
a_ij := function(b, p)
  local a, n, i, j, m, char;

  n := Size(b);
  a := NullMat(n,n)-M-1;
  char := Characteristic(b);

  if char = 0 then
    for i in [1..n] do
      for j in [1..n] do
        if i = j then
          a[i][i] := 2;
          continue;
        fi;

        if p[i] = 1 then 
          if IsInt(b[i][j]) and b[i][j]<=0 and b[i][i] <> 0 then 
            a[i][j] := b[i][j];
          else
            a[i][j] := infinity;
          fi;
        else
          if IsInt(b[i][j]) and b[i][j]<=0 and b[i][i] <> 0 then 
            a[i][j] := b[i][j];
          elif b[i][i]=0 and b[i][j] <> 0 then
            a[i][j] := -1;
          elif b[i][i]=0 and b[i][j] = 0 then
            a[i][j] := 0;
          else
            a[i][j] := infinity;
          fi;
        fi;
      od;
    od;
  fi;

  if char > 0 then
    for i in [1..n] do
      for j in [1..n] do
        if i = j then
          a[i][i] := 2;
          continue;
        fi;

        if b[i][i] = 2*Z(char)^0 and b[i][j] in GF(char) then
          if p[i] = 1 or (p[i] = -1 and from_ptoZ(char, b[i][j]) mod 2 = 0) then 
            a[i][j] := from_ptoZ(char, b[i][j]);
          else
            a[i][j] := from_ptoZ(char, b[i][j])-char;
          fi;
        elif b[i][i] = 2*Z(char)^0 and not b[i][j] in GF(char) then
          a[i][j] := 1-(3-p[i])*char/2;
        else
          if IsZero(b[i][j]) then
            a[i][j] := 0;
          else
            if p[i] = 1 then
              a[i][j] := 1-char;
            else
              a[i][j] := -1;
            fi;
          fi;
        fi;
      od;
    od;
  fi;
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

### This function computes the matrix after applying the reflection s_k
new_b := function(i, b, p, a)
  local m, j, k;
  m := NullMat(Size(b), Size(b));
  for j in [1..Size(b)] do
    for k in [1..Size(b)] do
      if j=i then
        m[j][k] := a[i][k]*b[i][i]-b[i][k];
      else
        if IsZero(b[i][j]) then 
          m[j][k] := b[j][k];
        else
          m[j][k] := a[i][j]*a[i][k]*b[j][i]*b[i][i]-a[i][j]*b[j][i]*b[i][k]-a[i][k]*b[j][i]*b[i][j]+b[i][j]*b[j][k];
        fi;
      fi;
    od;
  od;
  return normalization(m);
end;

new_p := function(i, p, a)
  return List([1..Size(p)], j->p[j]*p[i]^a[i][j]);
end;

is_odd_nondegenerate:= function(b, p, i)
  if Characteristic(b) = 0 then 
    return b[i][i] = 2 and p[i]=-1;
  else
    return b[i][i] = 2*Z(Characteristic(b))^0 and p[i]=-1;
  fi;
end;

super := function(file)
  local n, N, i, w, l, pos_roots, odd_nondegenerate, odd, heights, firstb, firstp, firsta, a, done, r, nabla;

  Read(file);

  Print("The parity is: ", p, "\nThe normalized matrix is:\n");
  b := normalization(b);
  Display(b);
  
  n := Size(b);
  N := Maximum(250, n^2);
  i := 1;
  w := IdentityMat(n, n);
  l := [1];
  
  pos_roots := [TransposedMat(w)[i]];
  odd_nondegenerate := [];
  odd := [];
  heights := [];
  firstb := ShallowCopy(b);
  firstp := ShallowCopy(p);
  firsta := a_ij(firstb, firstp);
  a := ShallowCopy(firsta);
  
  if ForAny(Cartesian([1..n],[1..n]), x->firsta[x[1]][x[2]] = infinity) then
    Print("The root system is infinite!\n");
    return fail;
  fi;
  
  done := true;
  
  while done do
  
    if is_odd_nondegenerate(b, p, i) then
      Add(odd_nondegenerate, TransposedMat(w)[i]);
    fi;
  
    if p[i] = -1 then 
      Add(odd, TransposedMat(w)[i]);
    fi;
  
    a := a_ij(b, p);

    if ForAny(Flat(a), x->x <= -M) then
      Print("The root system is infinite.\c\n");
      return fail;
    fi;
 
    if ForAny(Cartesian([1..n],[1..n]), x->a[x[1]][x[2]] = infinity) then
      Print("The root system is infinite!\n");
      Display(a);
      break;
    fi;
  
    r := s(i, a);
  
    w := w*r;
    b := new_b(i, b, p, a);
    p := new_p(i, p, a);
   
    # It suffices to check if one entry of <w> is positive
    i := First(Concatenation([1..i-1],[i+1..n]), x->ForAll(TransposedMat(w)[x], y->y>=0));
    
    if Size(l) > Maximum([N, Size(b)^2]) then
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
  
  nabla := Union(pos_roots, 2*odd_nondegenerate);
  
  Print("Longest word:\n", l, "\n");
  Print("Positive roots (", Size(pos_roots), "):\n", pos_roots, "\n");
  Print("Odd positive roots (", Size(odd), "):\n", odd, "\n");
  Print("Odd nondegenerate (", Size(odd_nondegenerate), "):\n", odd_nondegenerate, "\n");
  Print("The set Nabla+ is:\n",  nabla, "\n");
  Print("Super-dimension of the contragradient Lie superalgebra g(B,p): (", 
  2*Size(nabla)+2*n-Rank(b)-2*Size(odd), "|", 2*Size(odd), ")\n");
  
end;
