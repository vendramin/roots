###
### This script computes the root systems associated to a diagonal braiding
### Written by I. Angiono, L. Vendramin
### 

Read("roots.g");

### Bounds (global variables)
### These bounds are useful for computing truncated infinite root systems

### This "M" is the bound for the exponents of q_ii
M := 7;   
x := Indeterminate(Rationals, "x");
q := [];

### This "N" is the bound for the length of the longest word. 
### Use infinity for finite root systems.

### This is our function to compute the order of an element
### It returns infinity if a parameter is included
order := function(elm)
  if ExtRepPolynomialRatFun(elm*x*x^(-1)) = false then
    return infinity;
  elif ExtRepPolynomialRatFun(elm*x*x^(-1))[1] = [] then
    return Order(elm);
  else
    return infinity;
  fi;
end;

relations := function(file)
  local n, N, i, w, l, pos_roots, cartan_roots, firstq, firsta, a, done, heights, r, j, k, z, b, v, scalar;

  Read(file);

  LogTo();
  LogTo(Concatenation(file, ".log"));
 
  Print("The braiding is:\n");
  Display(q);
  
  n := Size(q);
  N := Maximum(250, n^2);
  i := 1;
  w := IdentityMat(n, n);
  l := [1];
  
  pos_roots := [TransposedMat(w)[i]];
  cartan_roots := [];
  heights := [];
  firstq := ShallowCopy(q);
  firsta := a_ij(firstq);
  a := ShallowCopy(firsta);
  done := true;
  
  while done do
 
    ### Is i a Cartan root?
    if is_cartan(q, a, i) then
      Add(cartan_roots, TransposedMat(w)[i]);
    fi;
  
    a := a_ij(q);
    r := s(i, a);

    if ForAny(Flat(a), x->x <= -M) then
      Print("The root system is infinite.\c\n");
      return;
    fi;
 
  
    Add(heights, order(q[i][i]));
  
    ### ax^n is represented as [[1, n], a]
    ### ax^0 is [[ ], a]
    ### ax^-n is false
    w := w*r;
    q := new_q(i, q, a);
  
    # It suffices to check if one entry of <w> is positive
    i := First(Concatenation([1..i-1],[i+1..n]), x->ForAll(TransposedMat(w)[x], y->y>=0));
    
    if i = fail or Size(l)>N then
      if Size(l)>N then
        Print("The root system is infinite!\n");
      fi;
      done := false;
    else 
      Add(pos_roots, TransposedMat(w)[i]);
      Add(l, i);
    fi;
     
  od;
  
  if infinity in heights then
    Print("Dimension: infinity\n");
  else
    Print("Dimension: ", Product(heights), "\n");
    Print("Longest word:\n", l, "\n");
    Print("Positive roots (", Size(pos_roots), "):\n", pos_roots, "\n");
    Print("Cartan roots (", Size(cartan_roots), "):\n", cartan_roots, "\n");
    Print("Heights: ", heights, "\n");
    Print("Relations:\n");
    
    ### Theorem 3.1
    ### Relations (3.3) and (3.1) 
    for i in [1..n] do
      if not order(firstq[i][i]) = infinity then
        Print("x_{", i, "}^{", order(firstq[i][i]), "}\n");
      fi;
    od;
    for v in cartan_roots do
      if v in IdentityMat(n) then
        continue;
      fi;
      j := Position(pos_roots, v);
      if not heights[j] = infinity then
        Print("x_{", v, "}^{", heights[j], "}\n");
      fi;
    od;
  ### Relations (3.4)
  for i in [1..n] do
    for j in [i+1..n] do
      if IsOne(-firstq[i][i]) and IsOne(-firstq[j][j]) and IsOne(-firstq[i][j]*firstq[j][i]) and not IdentityMat(n)[i]+IdentityMat(n)[j] in cartan_roots then
        Print("x_{", i, "\\,", j, "}^2\n");
      fi;
    od;
  od;
  
  ### Relations (3.2)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      if not IsOne(firstq[i][i]^(1-firsta[i][j])) then
        if (firsta[i][j] = 0 and i < j) then
          Print("(\\ad_c x_{", i, "})x_{", j, "}\n");
        elif firsta[i][j] <> 0 then
          Print("(\\ad_c x_{", i, "})^{", 1-firsta[i][j], "}x_{", j, "}\n");
        fi;
      fi;
    od;
  od;
  
  ### Relations (3.5)
  ### Typo in page 17: in the line before (3.5): \tilde{q_{ij}}\ne1 
  for i in [1..n] do
    for j in [i+1..n] do
      for k in [j+1..n] do
        if IsOne(-firstq[j][j]) and IsOne(firstq[i][k]*firstq[k][i]) and IsOne(firstq[i][j]*firstq[j][i]*firstq[j][k]*firstq[k][j]) and not IsOne(firstq[i][j]*firstq[j][i]) then
          Print("[x_{", i, "\\,", j, "\\,", k, "},x_{", j, "}]_c\n");
        fi;
      od;
    od;
  od;
  
  ### Relations (3.6)
  for i in [1..n] do
    for j in [1..n] do
      if i=j then
        continue;
      fi;
      if (IsOne(-firstq[j][j]) and order(firstq[i][i]*firstq[i][j]*firstq[j][i]) = 6 and not IsOne(-firstq[i][j]*firstq[j][i])) and (order(firstq[i][i])=3 or -firsta[i][j] > 2) then
        Print("[x_{", i, "\\,", i, "\\,", j, "},x_{", i, "\\,", j, "}]_c\n");
      fi;
    od;
  od;
  
  ### Relations (3.7)
  ### Typo in page 17: Remove \tilde{q}_{jk}\ne-1 from the line before (3.7)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        if order(firstq[i][i])=3 and firstq[i][j]*firstq[j][i] in [firstq[i][i],-firstq[i][i]] and IsOne(firstq[i][k]*firstq[k][i]) then
          if (IsOne(-firstq[j][j]) and IsOne(firstq[i][j]*firstq[j][i]*firstq[j][k]*firstq[k][j])) or (IsOne(firstq[j][j]*firstq[i][j]*firstq[j][i]) and IsOne(firstq[j][j]*firstq[k][j]*firstq[j][k])) then
            Print("[x_{", i, "\\,", i, "\\,", j, "\\,", k, "},x_{", i, "\\,", j, "}]_c\n");
          fi;
        fi;
      od;
    od;
  od;
  
  ### Relations (3.8)
  for i in [1..n] do
    for j in [i+1..n] do
      for k in [j+1..n] do
        if not IsOne(firstq[i][k]*firstq[k][i]) and not IsOne(firstq[i][j]*firstq[j][i]) and not IsOne(firstq[j][k]*firstq[k][j]) then
          scalar := (1-firstq[j][k]*firstq[k][j])/(firstq[k][j]*(1-firstq[i][k]*firstq[k][i]));
          if IsOne(DenominatorOfRationalFunction(scalar)) then
            Print("x_{", i, "\\,", j, "\\,", k, "}-", scalar, "[x_{", i, "\\,", k, "},x_{", j, "}]_c-", firstq[i][j]*(1-firstq[j][k]*firstq[k][j]), "x_{", j, "}x_{", i, "\\,", k, "}\n");
          else 
            Print("x_{", i, "\\,", j, "\\,", k, "}-\\frac{", NumeratorOfRationalFunction(scalar), "}{", DenominatorOfRationalFunction(scalar), "}[x_{", i, "\\,", k, "},x_{", j, "}]_c-", firstq[i][j]*(1-firstq[j][k]*firstq[k][j]), "x_{", j, "}x_{", i, "\\,", k, "}\n");
          fi;
        fi;
      od;
    od;
  od;
  
  ### Relations (3.9)
  ### Typo in page 18, line 1: $-q_{ii}^{-3}$ should be $-q_{ii}^3$
  ### Typo in page 17, line -2: $\tilde{q}_{jk}^{2}$ should be $-\tilde{q}_{jk}=-q_{ii}^2$
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        if not IsOne(firstq[i][k]*firstq[k][i]) then
          continue;
        fi;
        if (IsOne(-firstq[i][i]) and IsOne(-firstq[j][j]) and IsOne(firstq[i][j]^2*firstq[j][i]^2*firstq[j][k]*firstq[k][j])) or 
          (IsOne(-firstq[i][j]*firstq[j][i]) and IsOne(-firstq[j][j]) and order(firstq[i][i])=3 and IsOne(-firstq[i][i]*firstq[j][k]*firstq[k][j])) or 
          (IsOne(-firstq[k][k]) and IsOne(-firstq[j][j]) and IsOne(-firstq[j][k]*firstq[k][j]) and order(firstq[i][i])=3 and IsOne(-firstq[i][i]^2*firstq[i][j]*firstq[j][i])) or 
          (IsOne(-firstq[j][j]) and IsOne(firstq[i][j]*firstq[j][i]*firstq[i][i]^2) and IsOne(-firstq[j][k]*firstq[k][j]*firstq[i][i]^(-3))) or 
          (IsOne(-firstq[i][i]) and IsOne(-firstq[j][j]) and IsOne(-firstq[k][k]) and order(firstq[j][k]*firstq[k][j])=3 and firstq[i][j]*firstq[j][i] in [firstq[j][k]*firstq[k][j], -firstq[j][k]*firstq[k][j]]) then 
          Print("[[x_{", i, "\\,", j, "},x_{", i, "\\,", j, "\\,", k, "}]_c,x_{", j, "}]_c\n");
        fi;
      od;
    od;
  od;
  
  ### Relation (3.10)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        if IsOne(-firstq[i][i]) and IsOne(-firstq[j][j]) and IsOne(firstq[j][k]*firstq[k][j]*firstq[i][j]^3*firstq[j][i]^3) and IsOne(firstq[i][k]*firstq[k][i]) then
          Print("[[x_{", i, "\\,", j, "},[x_{", i, "\\,", j, "},x_{", i, "\\,", j, "\\,", k, "}]_c]_c,x_{", j, "}]_c\n");
        fi;
      od;
    od;
  od;
  
  ### Relations (3.11)
  ### Typo in page 18: before (3.11) replace $\tilde{q}_{ij}^2$ with $\tilde{q}_{ij}=q_{jj}^2$
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        if order(firstq[j][j])=3 and IsOne(firstq[j][k]*firstq[k][j]*firstq[j][j]^2) and IsOne(firstq[i][j]*firstq[j][i]*firstq[j][j]) and IsOne(firstq[i][k]*firstq[k][i]) then
          Print("[[x_{", i, "\\,", j, "\\,", k, "},x_{", j, "}]_c,x_{", j, "}]_c\n");
        fi;
      od;
    od;
  od;
   
  ### Relations (3.12)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        if order(firstq[k][k])=9 and firstq[j][j]=firstq[k][k] and IsOne(firstq[j][j]*firstq[i][j]*firstq[j][i]) and IsOne(firstq[j][j]*firstq[j][k]*firstq[k][j]) and IsOne(firstq[i][k]*firstq[k][i]) and IsOne(firstq[i][i]*firstq[k][k]^3) then
          Print("[[x_{", i, "\\,", i, "\\,", j, "},x_{", i, "\\,", i, "\\,", j, "\\,", k, "}]_c, x_{", i, "\\,", j, "}]_c\n");
        fi;
      od;
    od;
  od;
   
  ### Relations (3.13)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        if order(firstq[i][i])=9 and IsOne(firstq[i][i]*firstq[i][j]*firstq[j][i]) and IsOne(firstq[j][j]*firstq[j][k]*firstq[k][j]) and IsOne(firstq[j][j]*firstq[i][i]^4) and IsOne(firstq[i][k]*firstq[k][i]) and IsOne(firstq[k][k]*firstq[i][i]^3) then
          Print("[[x_{", i, "\\,", j, "\\,", k, "},x_{", j, "}]_c,x_{", k, "}]_c-", firstq[j][k]*Inverse(1+firstq[j][k]*firstq[k][j]), "[[x_{", i, "\\,", j, "\\,", k, "},x_{", k, "}]_c,x_{", j, "}]_c\n");
        fi;
      od;
     od;
  od;
   
  ### Relations (3.14)
  ### Typo: $\tilde{q}_{ij}^3$ should be $\tilde{q}_{ij}=q_{jj}^3$
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        if order(firstq[j][j])=4 and IsOne(firstq[j][j]*firstq[i][j]*firstq[j][i]) and IsOne(firstq[j][j]^3*firstq[j][k]*firstq[k][j]) and IsOne(firstq[i][k]*firstq[k][i]) then
          Print("[[[x_{", i, "\\,", j, "\\,", k, "},x_{", j, "}]_c,x_{", j, "}]_c,x_{", j, "}]_c\n");
        fi;
      od;
    od;
  od;
   
  ### Relations (3.15)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        if IsOne(-firstq[i][i]) and IsOne(-firstq[i][j]*firstq[j][i]) and IsOne(firstq[i][k]*firstq[k][i]) and IsOne(firstq[j][j]*firstq[j][k]*firstq[k][j]) and not IsOne(-firstq[j][j]) then
          Print("[x_{", i, "\\,", j, "},x_{", i, "\\,", j, "\\,", k, "}]_c\n");
        fi;
      od;
    od;
  od;
  
  ### Relations (3.16)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        if IsOne(-firstq[i][i]) and IsOne(-firstq[k][k]) and IsOne(firstq[i][k]*firstq[k][i]) and order(firstq[i][j]*firstq[j][i])=3 
          and firstq[j][j] in [firstq[i][j]*firstq[j][i], -firstq[i][j]*firstq[j][i]] and firstq[j][j] = -firstq[j][k]*firstq[k][j] then
          Print("[x_{", i, "},x_{", j, "\\,", j, "\\,", k, "}]_c-(", (1+firstq[j][j]^2)*firstq[k][j]^(-1), ")[x_{", i, "\\,", j, "\\,", k, "},x_{", j, "}]_c-(", (1+firstq[j][j]^2)*(1+firstq[j][j])*firstq[i][j], ")x_{", j, "}x_{", i, "\\,", j, "\\,", k, "}\n");
        fi;
      od;
    od;
  od;
  
  ### Relations (3.17)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        for l in [1..n] do
          if l in [i,j,k] then 
            continue;
          fi;
          if IsOne(firstq[j][j]*firstq[i][j]*firstq[j][i]) and IsOne(firstq[j][j]*firstq[j][k]*firstq[k][j]) and 
            IsOne(-firstq[k][k]) and IsOne(firstq[i][k]*firstq[k][i]) and IsOne(firstq[i][l]*firstq[l][i]) and
            IsOne(firstq[j][l]*firstq[l][j]) and IsOne(firstq[l][l]*firstq[l][k]*firstq[k][l]) and IsOne(firstq[j][k]^2*firstq[k][j]^2*firstq[l][k]*firstq[k][l]) then
            Print("[[[x_{", i, "\\,", j, "\\,", k, "\\,", l, "},x_{", k, "}]_c,x_{", j, "}]_c,x_{", k, "}]_c\n");
          fi;
        od;
      od;
    od;
  od;
   
  ### Relations (3.18)
  ### Typo in page 18: $\tilde{q}_{lk}=\tilde{q}_{jk}^3$ should be $\tilde{q}_{lk}=\tilde{q}_{jk}$
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        for l in [1..n] do
          if l in [i,j,k] then 
            continue;
          fi;
          if order(firstq[j][j]) in [4, 6] and IsOne(firstq[j][j]*firstq[i][j]*firstq[j][i]) and IsOne(firstq[j][j]*firstq[k][j]*firstq[j][k]) and 
            IsOne(-firstq[i][i]) and IsOne(-firstq[k][k]) and IsOne(firstq[i][k]*firstq[k][i]) and IsOne(firstq[i][l]*firstq[l][i]) and 
            IsOne(firstq[j][l]*firstq[l][j]) and IsOne(firstq[l][k]*firstq[k][l]*firstq[j][k]^(-1)*firstq[k][j]^(-1)) then
            Print("[[x_{", i, "\\,", j, "\\,", k, "},[x_{", i, "\\,", j, "\\,", k, "\\,", l, "},x_{", k, "}]_c]_c,x_{", j, "\\,", k, "}]_c\n");
          fi;
        od;
      od;
    od;
  od;
   
  ### Relations (3.19)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        for l in [1..n] do
          if l in [i,j,k] then 
            continue;
          fi;
          if IsOne(firstq[i][k]*firstq[k][i]) and IsOne(firstq[i][l]*firstq[l][i]) and 
            IsOne(firstq[j][l]*firstq[l][j]) and IsOne(firstq[l][l]*firstq[l][k]*firstq[k][l]) and
            IsOne(firstq[k][k]*firstq[l][k]*firstq[k][l]) and IsOne(firstq[i][i]*firstq[i][j]*firstq[j][i]) and
            IsOne(-firstq[j][j]) and IsOne(firstq[l][l]^3*firstq[i][i]^2) then
            Print("[[[x_{", i, "\\,", j, "\\,", k, "},x_{", j, "}]_c,[x_{", i, "\\,", j, "\\,", k, "\\,", l, "},x_{", j, "}]_c]_c,x_{", j, "\\,", k, "}]_c\n");
          fi;
        od;
      od;
    od;
  od;
  
  ### Relations (3.20)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        for l in [1..n] do
          if l in [i,j,k] then 
            continue;
          fi;
           if IsOne(firstq[i][k]*firstq[k][i]) and IsOne(firstq[i][l]*firstq[l][i]) and 
            IsOne(firstq[j][l]*firstq[l][j]) and IsOne(firstq[i][i]*firstq[i][j]*firstq[j][i]) and IsOne(-firstq[k][k]) then
            if (IsOne(firstq[j][j]^2*firstq[i][j]*firstq[j][i]) and IsOne(firstq[j][j]^3*firstq[l][l]) 
              and IsOne(firstq[l][l]*firstq[k][l]*firstq[l][k]) and IsOne(firstq[j][j]*firstq[j][k]*firstq[k][j])) or 
              (IsOne(-firstq[l][l]*firstq[i][i]) and IsOne(firstq[l][l]*firstq[l][k]*firstq[k][l]) and
              IsOne(-firstq[j][j]) and IsOne(-firstq[j][k]*firstq[k][j])) then
              Print("[[x_{", i, "\\,", j, "\\,", k, "\\,", l, "},x_{", j, "}]_c,x_{", k, "}]c-(", firstq[j][k]*(firstq[i][j]^(-1)*firstq[j][i]^(-1)-firstq[j][j]), ")[[x_{", i, "\\,", j, "\\,", k, "\\,", l, "},x_{", k, "}]_c,x_{", j, "}]c\n");
            fi;
          fi;
        od;
      od;
    od;
  od;
  
  ### Relations (3.21)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        if IsOne(firstq[j][k]*firstq[k][j]) and order(firstq[i][i]) = 3 and IsOne(firstq[i][j]*firstq[j][i]*firstq[i][i]^2) 
          and IsOne(-firstq[i][k]*firstq[k][i]*firstq[i][i]^2) then
          Print("[x_{", i, "},[x_{", i, "\\,", j, "},x_{", i, "\\,", k, "}]_c]_c+(", firstq[j][k]*firstq[i][k]*firstq[j][i], ")[x_{", i, "\\,", i, "\\,", k, "},x_{", i, "\\,", j, "}]_c+(", firstq[i][j], ")x_{", i, "\\,", j, "}x_{", i, "\\,", i, "\\,", k, "}\n");
        fi;
      od;
    od;
  od;
  
  ### Relations (3.22)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      for k in [1..n] do
        if k in [i,j] then
          continue;
        fi;
        if IsOne(-firstq[j][j]) and IsOne(-firstq[k][k]) and IsOne(-firstq[j][k]*firstq[k][j]) and order(firstq[i][i])=3 
          and IsOne(-firstq[i][j]*firstq[j][i]*firstq[i][i]^2) and IsOne(firstq[i][k]*firstq[k][i]) then
          Print("[x_{", i, "\\,", i, "\\,", j, "\\,", k, "},x_{", i, "\\,", j, "\\,", k, "}]_c\n");
        fi;
      od;
    od;
  od;
  
  ### Relations (3.23)
  for i in [1..n] do
    for j in [i+1..n] do
      if firsta[i][j] <= -2 and firsta[j][i] <= -2 then
        Print("(", (1-firstq[i][j]*firstq[j][i])*firstq[j][j]*firstq[j][i], ")[x_{", i, "},[x_{", i, "\\,", j, "},x_{", j, "}]_c]_c-(", (1+firstq[j][j])*(1-firstq[j][j]*firstq[i][j]*firstq[j][i]), ")x_{", i, "\\,", j, "}^2\n");
      fi;
    od;
  od;
  
  ### Relations (3.24)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      if (-firsta[i][j] in [4,5]) or (firsta[i][j] = -3 and order(firstq[i][i])=4 and IsOne(-firstq[j][j])) then 
        Print("[x_{", i, "},x_{3\\alpha_{", i, "}+2\\alpha_{", j, "}}]_c-(", (1-firstq[i][i]*firstq[i][j]*firstq[j][i]-firstq[i][i]^2*firstq[i][j]^2*firstq[j][i]^2*firstq[j][j])/(firstq[j][i]*(1-firstq[i][i]*firstq[i][j]*firstq[j][i])), ")x_{", i, "\\,", i, "\\,", j, "}^2\n");
      fi;
    od;
  od;
  
  ### Relations (3.25)
  for i in [1..n] do 
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      if not 4*IdentityMat(n)[i]+3*IdentityMat(n)[j] in pos_roots and (IsOne(-firstq[j][j]) or firsta[j][i] <= -2) 
        and (firsta[i][j] <= -3 or (firsta[i][j] = -2 and order(firstq[i][i])=3)) and 3*IdentityMat(n)[i]+2*IdentityMat(n)[j] in pos_roots then
        Print("[x_{3\\alpha_{", i, "}+2\\alpha_{", j, "}},x_{", i, "\\,", j, "}]_c\n");
      fi;
    od;
  od;
  
  ### Relations (3.26)
  for i in [1..n] do 
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      if not IsOne(firstq[i][i]^3*firstq[i][j]*firstq[j][i]) and not IsOne(firstq[i][i]^4*firstq[i][j]*firstq[j][i]) and 
        3*IdentityMat(n)[i]+2*IdentityMat(n)[j] in pos_roots and not 5*IdentityMat(n)[i]+3*IdentityMat(n)[j] in pos_roots then
        Print("[x_{", i, "\\,", i, "\\,", j, "},x_{3\\alpha_{", i, "}+2\\alpha_{", j, "}}]_c\n");
      fi;
    od;
  od;
  
  ### Relations (3.27)
  for i in [1..n] do 
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      if 4*IdentityMat(n)[i]+3*IdentityMat(n)[j] in pos_roots and not 5*IdentityMat(n)[i]+4*IdentityMat(n)[j] in pos_roots then
        Print("[x_{4\\alpha_{", i, "}+3\\alpha_{", j, "}},x_{", i, "\\,", j, "}]_c\n");
      fi;
    od;
  od;
  
  ### Relations (3.28) 
  for i in [1..n] do 
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      if 5*IdentityMat(n)[i]+2*IdentityMat(n)[j] in pos_roots and not 7*IdentityMat(n)[i]+3*IdentityMat(n)[j] in pos_roots then
        Print("[[x_{", i, "\\,", i, "\\,", i, "\\,", j, "},x_{", i, "\\,", i, "\\,", j, "}]_c,x_{", i, "\\,", i, "\\,", j, "}]_c\n");
      fi;
    od;
  od;
  
  ### Relations (3.29)
  for i in [1..n] do
    for j in [1..n] do
      if i = j then
        continue;
      fi;
      if IsOne(-firstq[j][j]) and 5*IdentityMat(n)[i]+4*IdentityMat(n)[j] in pos_roots then
        z := firstq[i][j]*firstq[j][i];
        a := (1-z)*(1-firstq[i][i]^4*z^3)-(1-firstq[i][i]*z)*(1+firstq[i][i])*firstq[i][i]*z;
        b := (1-z)*(1-firstq[i][i]^6*z^5)-a*firstq[i][i]*z;
        Print("[x_{", i, "\\,", i, "\\,", j, "}, x_{4\\alpha_{",i,"}+3\\alpha_{",j, "}}]_c-(", (b-(1+firstq[i][i])*(1-firstq[i][i]*z)*(1+z+firstq[i][i]*z^2)*firstq[i][i]^6*z^4)/(a*firstq[i][i]^3*z^2*firstq[j][i]), ")x_{3\\alpha_{", i, "}+2\\alpha_{", j, "}}^2\n");
      fi;
    od;
  od;
  fi;
end; 
