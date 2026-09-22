# Maquette Quarto — à tester, ne remplace rien pour l'instant

Ce dossier est une proposition complète pour basculer ton site de
`postcards`/`distill` vers **Quarto**. Il ne touche à rien du site actuel
(`index.Rmd`, `docs/`, etc.) : c'est un dossier à part que tu peux essayer,
comparer, puis adopter (ou pas).

## Pour le tester (2 minutes)

1. Installe Quarto (une seule fois) : https://quarto.org/docs/get-started/
   — un simple installeur, pas un package R.
2. Dans RStudio, ouvre ce dossier `quarto-preview/` comme dossier de travail.
3. Soit tu cliques sur "Render" (RStudio détecte automatiquement le projet
   Quarto), soit dans un terminal, à l'intérieur de ce dossier :
   ```
   quarto preview
   ```
   Ça ouvre le site dans ton navigateur, avec rechargement automatique à
   chaque modification.

## Ce qui est repris à l'identique

* Le texte de toutes les pages (accueil, About me, Research, Teaching) —
  copié tel quel depuis `index.Rmd`, `curriculum.Rmd`, `research.Rmd`,
  `teaching.Rmd`.
* La structure de navigation (mêmes 3 liens : About me / Research / Teaching).
* La page d'accueil utilise le template Quarto **"jolla"** — le même nom, la
  même disposition (photo ronde + nom + liens) que `postcards::jolla` que tu
  utilises aujourd'hui. C'est une fonctionnalité native de Quarto, pas une
  extension tierce.
* Les couleurs et polices de `postcards.css` (vert sauge / bleu marine,
  polices Amiri/Bitter/DM Mono) sont reprises dans `styles.css` — adaptées
  aux classes HTML de Quarto (Bootstrap), qui sont différentes de celles de
  Distill. Le rendu exact peut donc légèrement différer par endroits ; à
  ajuster à l'œil une fois que tu l'auras sous les yeux.

## Ce qui reste à faire si tu adoptes cette version pour de bon

* Les liens vers `slides/Theme2_Slides.html`, `slides/Theme4_Slides.html`
  et `phdtopic/Slides.html` pointent vers les mêmes chemins relatifs
  qu'aujourd'hui. Ils ne fonctionnent pas depuis ce dossier `quarto-preview/`
  isolé (les dossiers `slides/` et `phdtopic/` sont à la racine du dépôt),
  mais fonctionneront normalement une fois ces fichiers déplacés à la racine
  du dépôt, aux côtés de ces mêmes dossiers.
* Les icônes "Google Scholar" et "ORCID" sont pour l'instant en texte seul
  (pas d'icône Bootstrap officielle pour ces deux-là) ; possible d'ajouter
  l'extension Quarto "academicons" plus tard si tu veux des icônes dédiées.
* Le fichier `quarto-publish.workflow.yml` est prêt pour automatiser le
  rendu via GitHub Actions (à chaque `push`, comme discuté) — à copier
  toi-même vers `.github/workflows/quarto-publish.yml`, je n'ai pas le droit
  d'écrire directement dans ce dossier depuis cette session.
* Le `renv.lock` du site actuel n'a pas d'équivalent obligatoire côté Quarto :
  Quarto n'est pas un package R (c'est un outil en ligne de commande), donc
  aucune dépendance R n'est nécessaire pour ces pages qui ne contiennent pas
  de code R. Si tu ajoutes des chunks R plus tard, un `renv.lock` séparé
  pourra être utile.

## Pourquoi je n'ai pas pu te montrer un rendu déjà fait

Je n'ai pas pu installer Quarto dans cet environnement cloud pour te générer
un aperçu déjà rendu : l'accès réseau y est restreint (CRAN, quarto.org et
GitHub bloqués par la politique réseau de cette session). Plutôt que de te
proposer une image que j'aurais reconstituée à la main — donc pas fiable à
100 % — je préfère te donner les vraies sources, à rendre toi-même en une
poignée de secondes une fois Quarto installé, pour voir le résultat exact.
