g = Mod(6, 682492462409094395392022581537473179285250139967739310024802121913471471);
A = 245036439927702828116237663546936021015004354074422410966568949608523157;

/* L'ordre du groupe multiplicatif est (682492462409094395392022581537473179285250139967739310024802121913471471 - 1)
 * qui est composé, on va employer Pohlig-Hellman pour décomposer en sous-problèmes de calcul de log discret
 * à résoudre via Baby Step Giant Step.
 */

/* Algo Baby Step-Giant Step:
 * g: un générateur du groupe cyclique G=<g> d'ordre p
 * A: élément du groupe G
 * Retourne l'unique élément x du groupe tel que g^x = A

 * Note : l'introduction avec mapput des puissances de g est dure longtemps pour p grand (et consomme beaucoup de RAM).
 * À employer avec p premier "petit" (testé avec p de l'ordre de 9 chiffres).
 */
Shanks_DLP(g, A, p) = {
  my(B, a, lookup_table, i, j, fac);
  B = ceil(sqrt(p));
  lookup_table = Map(); \\ table de hachage
  a = Mod(1, g.mod);
  for(i = 0, B, mapput(lookup_table, lift(a), i); a = a * g;);
  a = A;
  fac = g^(-B);
  for(i = 0, B,
    if(mapisdefined(lookup_table, lift(a), &j), return (i * B + j));
    a = a * fac;
  ); 
}

/*
\\ Pollard Rho DLP Step function
step(x, a, b, g, A) = {
  my(class_x);
  class_x = lift(x) % 3; \\ Le groupe ambiant est partitionné en 3 sous-ensembles, on récupère la classe de l'élément x
  if(class_x == 0, a *= 2; b *= 2; x = x^2);
  if(class_x == 1, a += 1; x = g * x);
  if(class_x == 2, b += 1; x = A * x);
  [x, a, b];
} 

RhoPollard_DLP(g, A, n) = {
  my(xk, ak, bk, x2k, a2k, b2k, q);
  \\p = znorder(g);
  q = (n - 1) / 2;
  ak = 1;
  bk = 1;
  ak = random(n);
  bk = random(n);
  xk = Mod(g^ak * A^bk, n);
  x2k = xk;
  a2k = ak;
  b2k = bk;
  \\ recherche d'une collision x_k=x_2k
  while(1,
    \\ pas de tortue
    [xk, ak, bk] = step(xk, ak, bk, g, A);
    \\ pas de lièvre
    [x2k, a2k, b2k] = step(x2k, a2k, b2k, g, A);
    [x2k, a2k, b2k] = step(x2k, a2k, b2k, g, A);
    if(xk == x2k,
      \\ On obtient g^(a2k-ak) = A^(bk-b2k)
      \\ dont on tire: A=g^a avec a = (ak-a2k)/(b2k-bk)
      return (lift(Mod((ak - a2k) / (b2k - bk), n))););
  );
}
*/

/* Pohlig-Hellman dans le cas où l'ordre du groupe est une puissance d'un nombre premier (ici p^e)
 * Retourne l'unique entier satisfaisant g^x = A mod p^e
 */
PH_prime_order(g, A, p, e) = {
  my(x, k, d, Ak, gamma_);
  x = 0; \\ = x0 + x1 p + x2 p^2 + ... + x_{e-1} p ^{e-1}
  gamma_ = g^(p^(e - 1)); \\ d'ordre p via Lagrange
  for(k = 0, e - 1,
    Ak = (g^(-x) * A)^(p^(e - 1 - k));
    d = Shanks_DLP(gamma_, Ak, p);
    x += p^k * d;
  );
  x;
}

/* Décomposition de Pohlig-Hellman pour un groupe d'ordre non premier.
 * Fonctionne sur le corps multiplicatif (Z/nZ)* uniquement.
 * Cours + https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm#The_general_algorithm

 * Entrée g générateur d'un groupe G d'ordre n, A élément du groupe.
 * Retourne a=log_g(A), tel que g^a=A dans G.
 */
Pohlig_Hellman(g, A) = {
  my(i, r, n, group_order, xi, pi, ei, ni, gi, Ai, X);
  n = g.mod;
  A = Mod(A, n);
  group_order = n - 1; 
  prime_factors = factor(group_order);
  X = List();
  r = matsize(prime_factors)[1];
  for(i = 1, r,
    pi = prime_factors[i, 1]; \\ facteur premier
    ei = prime_factors[i, 2]; \\ multiplicité de pi
    ni = group_order / pi^ei;
    gi = g^ni;
    Ai = A^ni;
    \\ On décompose le problème en log discret sur un groupe d'ordre plus petit
    \\ xi vérifie gi^xi = Ai mod pi^ei
    xi = PH_prime_order(gi, Ai, pi, ei);
    listput(X, Mod(xi, pi^ei));
  );
  lift(Mod(lift(chinese(X)), n)); \\ On applique le CRT sur les paires de congruences déterminées au dessus
}

/*
//tests:
n = 604604729;
g= Mod(7894352216, n);
A=Mod(355407489,n);
print(Shanks_DLP(g, A));
print(RhoPollard_DLP(g, A));
*/
print(Pohlig_Hellman(g, A));
