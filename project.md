# Project outline


## Estimates based on gene orthology

__Prerequisites__:

- given two species A and B
- given a table of gene expression estimates for one or more biol. samples for species A and B
- given a table of gene orthologs between species A and B
- given a table of analogy information between biol. samples of species A and B

__Restrict__:

- select subset of orthologous genes between species A and B 

__Procedure__:

- map gene expression estimates from A to B and vice versa
- evaluate agreement between true expression values and mapped expression values from orthologous genes

__Produce__:

- distinguish between *positive* and *negative* mappings (analogous biol. samples)
- distinguish between *all* orthologs and *active* subset (active: TPM >= 1)
- compute the following metrics: accuracy, r2 score, kendall tau correlation
- Accuracy: weighted accuracy for "classification" task (gene off == TPM < 1 == class 0; gene on == TPM >= 1 == class 1)
- R2 score: weighted R2 score for all orthologs, unweighted for active subset
- Kendall Tau: consistency of TPM-based ranking for all orthologs and for active subset
