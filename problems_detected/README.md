# Problèmes détectés

Ce dossier centralise les problèmes détectés automatiquement par l'agent.

## Fichiers

- `latest-report.md` : dernier rapport généré automatiquement.
- `history/` : historique des rapports horodatés.

## Fonctionnement

Le workflow GitHub Actions `Problem Agent` (fichier `.github/workflows/problem-agent.yml`) :

1. Vérifie les liens du site.
2. Génère un rapport Markdown.
3. Met à jour `problems_detected/latest-report.md`.
4. Archive une copie dans `problems_detected/history/`.
5. Ouvre/actualise une issue GitHub intitulée **[Problem Agent] Rapport automatique** si des problèmes sont trouvés.
