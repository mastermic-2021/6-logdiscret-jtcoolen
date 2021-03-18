# Logarithme discret

## Problème

En utilisant les méthodes du cours, calculer le logarithme
discret de A modulo g, où
```
A = 245036439927702828116237663546936021015004354074422410966568949608523157;
g = Mod(6, 682492462409094395392022581537473179285250139967739310024802121913471471);
```

On s'interdit donc l'utilisation de la fonction ``znlog`` de Pari.


## Format

Vous devez écrire un programme `main.gp` en Pari/GP dont l'exécution via
```
gp -fq < main.gp
```
affiche (uniquement) la solution cherchée.

taper `make check` permet de vérifier que votre solution est bonne.

Veillez à mettre des commentaires sur votre démarche et sa validité
dans le fichier ``main.gp``.

De manière alternative, vous pouvez écrire un programme `main.c` qui
fasse la même chose.

## Validation

Il reste à faire un commit et à envoyer votre solution pour l'enregistrer
```
git add main.gp
git commit -m 'ma solution'
git push
```

