# ACE Cell-Cell Communication Analysis

Cell-cell communication inference using CellChat to identify altered
ligand-receptor interactions between brain cell types in ACE-exposed
individuals.

## Method

Individuals are split into ACE-high and ACE-low groups (median split on
the specified phenotype, within each sex). CellChat infers communication
probability for all ligand-receptor pairs in CellChatDB.human. Differential
interaction analysis compares signaling between groups, with focus on
microglia-neuron, astrocyte-neuron, and oligodendrocyte-neuron axes.

## Cohort Directories

| Cohort | Entry Point |
|--------|-------------|
| Tsai | `Tsai/aceCellChatT.sh` |
| DeJager | `DeJager/aceCellChatDJ.sh` |

## Key Interaction Axes

| Sender | Receiver | Pathways of Interest |
|--------|----------|---------------------|
| Mic | In-PV_Basket | Complement (C3-C3AR1), CX3C (CX3CL1-CX3CR1) |
| Mic | Exc | TNF, IL1, TGFB |
| Ast | In-PV_Basket | Glutamate, GABA signaling |
| Oli | Exc | NRG (NRG1-ERBB), MAG, MOG |
