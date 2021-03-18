default(parisizemax, 70m);

g = Mod(6, 682492462409094395392022581537473179285250139967739310024802121913471471);
A = 245036439927702828116237663546936021015004354074422410966568949608523157;

\\ L'ordre du groupe est (682492462409094395392022581537473179285250139967739310024802121913471471 - 1)
\\ qui est composé, on va employer Pohlig-Hellman pour décomposer en sous-problèmes.


\\ g: un générateur du groupe cyclique G=<g> d'ordre n
\\ A: élément du groupe G
\\ Retourne l'unique élément x du groupe tel que g^x = A
\\ Note : l'introduction avec mapput des puissances de g est très longue pour n grand.
\\ À employer avec n premier de préférence et évidemment n "petit" (testé avec n de l'ordre
\\ de 9 chiffres).
Shanks_DLP(g, A, n) = {
  my(B, a, lookup_table, i, j, fac);
  B = ceil(sqrt(n));
  lookup_table = Map(); \\ table de hachage
  a = Mod(1, n);
  for(i = 0, B, mapput(lookup_table, lift(a), i); a = a * g;);
  a = A;
  j = 0;
  fac = g^(-B);
  for(i = 0, B,
    if(mapisdefined(lookup_table, lift(a), &j), return (i * B + j));
    a = a * fac;
  ); 
}

\\ Retourne un tuple [x, a, b] où x est le terme
next_term(x, a, b, g, A) = {
  my(class_x);
  class_x = lift(x) % 3;
print(x);
  if(class_x == 0, a *= 2; b *= 2; x = x^2);
  if(class_x == 1, a += 1; x = g * x);
  if(class_x == 2, b += 1; x = A * x);
  [x, a, b];
} 

RhoPollard_DLP(g, A, p) = {
  my(xk, ak, bk, x2k, a2k, b2k, t);
  q = (p - 1);  \\ sub group
  xk = Mod(g * A, p);
  ak = Mod(1, p);
  bk = Mod(1, p);
  x2k = xk;
  a2k = ak;
  b2k = bk;
  \\ recherche d'une collision x_k=x_2k
  while(1,
    \\ pas de tortue:
    [xk, ak, bk] = next_term(xk, ak, bk, g, A);
    \\ pas de lièvre
    [x2k, a2k, b2k] = next_term(x2k, a2k, b2k, g, A);
    [x2k, a2k, b2k] = next_term(x2k, a2k, b2k, g, A);
    if(xk == x2k,
      \\ On obtient g^(a2k-ak) = A^(bk-b2k)
      \\ dont on tire: A=g^a avec a = (a2k-a)/(bk-b2k)
      lift((a2k - ak) / (bk - b2k)););
  );
}

Pohlig_Hellman(n) = {
  my(i);
  prime_factors = factor(n);
  r = matsize(prime_factors)[1];
  for(i = 1, r,
    pi = prime_factors[r, 1];
    mi = prime_factors[r, 2];
    gi = g^(n/pi^mi);
    hi = h^(n/pi^mi);
);
    

}

\\print(DLP_solve(g, A, 682492462409094395392022581537473179285250139967739310024802121913471471));
n = 604604729;
g= Mod(7894352216, n);
A=355407489;
print(RhoPollard_DLP(g, A, n));
