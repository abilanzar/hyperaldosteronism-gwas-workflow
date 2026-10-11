# Public-release boundary

## Safe to publish

- Generic workflow code with portable configuration.
- Synthetic or formally public example data.
- Tests and expected outputs generated from synthetic data.
- High-level descriptions of methods already public in a paper, preprint, poster, or approved abstract.

## Do not publish

- Participant-level genotype, phenotype, covariate, or EHR data.
- Controlled-access or redistribution-restricted summary statistics.
- Credentials, tokens, cluster host details, private paths, or internal project IDs.
- Unpublished locus tables, effect estimates, manuscript drafts, collaborator comments, or reviewer correspondence without approval.
- Internal logs that expose usernames, paths, job metadata, or protected dataset names.
- Code copied from collaborators or institutional pipelines without permission and license review.

## Before creating a public GitHub remote

1. Confirm the manuscript/preprint status and collaborator approval.
2. Review the repository history, not only the current files.
3. Run a secret and large-file scan.
4. Verify that every example is synthetic or redistributable.
5. Remove local absolute paths and institutional identifiers.
6. Add an explicit license only after ownership is clear.

