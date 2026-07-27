---
layout: post
title: "Simulations for Biological Inference"
date: 2026-06-05 10:00:00
description: "Using simulations to find the biology hidden in experimental data."
tags: statistics simulation mitochondria
categories: research
thumbnail: assets/img/posts/mtDNA-smFISH-simulations/mtDNA-smFISH-thumbnail.webp
giscus_comments: true
toc: true
---

Modern biological experiments can generate an astonishing amount of detail about living systems. Single-cell technologies allow us to profile thousands of individual cells at once, imaging approaches let us visualize molecules inside tissues, and perturbation experiments make it possible to ask how biological systems respond when specific genes or pathways are disrupted. This is especially powerful for mitochondrial biology, where mutations and their effects can vary between cells, tissues, developmental stages, and disease contexts.

Yet more detailed measurements do not automatically make interpretation easier. Observations can show us a pattern, but they cannot tell us on their own which part of that pattern is biology, which part is noise, and which part is a consequence of the measurement process itself.

This is where simulations become powerful. To interpret an observation, we often need to compare it against a counterfactual: what would we have observed if only the simplest explanation were operating? That comparison rarely comes from the experiment itself. We have to build it, and simulations enable us to do that by making our assumptions explicit, generating the data those assumptions predict, and then asking whether what we actually observed is surprising.

> **Experiments show us what happened. Simulations ask what else could have happened.**

### A deliberately simplified world

---

The key step in simulation-based inference is building what I think of as a deliberately simplified world, also known as a null model. Here, we decide which processes are allowed to operate and which ones are left out. We might allow random sampling, measurement error, or differences in detection sensitivity, while leaving out the more interesting mechanisms we are trying to test. This lets us isolate the different factors that can influence our measurements and understand the kinds of patterns they can produce on their own.

That matters because biological data can be surprisingly good at producing convincing-looking patterns by chance alone. A difference between groups might just reflect how deeply each was measured, or a shift in a distribution might reflect technical variation rather than biology.

Simulations let us make these expectations explicit: generate the data we'd expect under the simplified model, then compare it to what we actually observed. If the observed pattern falls inside that expectation, the simulation protects us from over-interpreting a result that chance alone could explain. If it falls outside, the simulation tells us where the simple explanation breaks down and what's left to account for.

### The quantitative nature of mtDNA heteroplasmy

---

Mitochondria carry their own genome: mitochondrial DNA (mtDNA) which exists in hundreds to thousands of copies per cell. The population of mtDNA in a cell can be heterogeneous, with copies carrying different mutant variants, a phenomenon termed heteroplasmy. In the clinical context, heteroplasmy is the presence of a disease-causing variant alongside the normal, wild-type mtDNA and the proportion of the mutant copies in the cell determines whether and how severely a cell is affected. The same mutation can be silent or damaging in a cell depending on the heteroplasmy level (proportion of mutant mtDNA) of the cell.

This makes mitochondrial genetics unusually quantitative. It is not enough to ask whether a mutation is present or absent. We often need to ask how much of it is present, how the proportion varies between cells, and how it changes as cells divide, differentiate or respond to stress. Because cells contain finite numbers of mitochondrial genomes, and experiments measure only a subset of them, heteroplasmy patterns can shift for stochastic, technical, and biological reasons.

### When the simpler explanation is enough

---

The mitochondrial genome is compact and encodes only a small set of proteins, all of which are components of the oxidative phosphorylation system. Nuclear genes play a central role in mitochondrial genome biology. They regulate how mtDNA is replicated, packaged, and maintained, but establishing direct causal links between specific genes and heteroplasmy dynamics at the level of individual cells has been difficult.

Our recently published study in Nature Structural and Molecular Biology describes MitoPerturb-Seq, a technique designed to investigate these links directly. It perturbs a set of nuclear genes, then measures mitochondrial DNA copy number, heteroplasmy, chromatin accessibility, and gene expression in thousands of individual cells.

The surprising result was that mean heteroplasmy was unaffected across all these perturbations but cells with _Tfam_, _Opa1_, and _Polg_ knockdown exhibited a significant reduction in mtDNA copy number, along with increased heteroplasmy variance. At first glance, this observation can suggest that disrupting these genes makes mitochondrial genetic variation less stable across cells.

However, the increase in heteroplasmy variance can also be confounded by the reduced mtDNA copy number we observe in cells where these three genes are perturbed. When fewer mtDNA molecules are present, or captured by the assay, heteroplasmy estimates become noisier. Cells can then appear more variable even if the underlying distribution of mutant and wild-type mtDNA has not changed in an unexpected way. The observed pattern can have two possible explanations: altered heteroplasmy dynamics, or the expected consequence of measuring fewer mitochondrial genomes.

This was the ambiguity the simulations were designed to test. I built a null model that asked what heteroplasmy variance would look like if mtDNA depth alone was responsible for producing the patterns. Starting from the control heteroplasmy distribution, the simulation repeatedly resampled cells at the reduced depths observed in _Tfam_, _Opa1_, and _Polg_ knockdown cells. This created the counterfactual we needed: a simplified world where mtDNA depth dropped, but no additional biological process was added to increase heteroplasmy variance.

When we compared the observed variances to that simulated expectation, the apparent increase in heteroplasmy variance was not surprising. The depleted perturbations sat within the range expected from reduced mtDNA depth alone, no additional mechanism required. That doesn't mean the perturbations did nothing: they clearly depleted mtDNA and altered the cellular response to mitochondrial stress. But the increase in heteroplasmy variance wasn't a separate effect needing its own explanation. Once mtDNA depth was accounted for, the simpler model already covered it.

### When the simpler explanation falls short

---

The MitoPerturb-Seq study measured heteroplasmy indirectly, through sequencing reads from thousands of dissociated cells. The mtDNA-smFISH study instead visualizes mitochondrial DNA molecules inside intact cells and tissues. Using sequence-specific fluorescent probes, mtDNA-smFISH allows wild-type and mutant mitochondrial genomes to be counted directly in whole tissues, while preserving spatial and cellular context. This enabled visualization of heteroplasmy distributions in actively dividing cells such as the neural stem cells.

Applying this technique to developing neural stem cells (NSCs) in a heteroplasmic _Drosophila_ line enabled the visualization and quantification of heteroplasmy dynamics as the NSCs divided to give rise to progeny cells. We observed that the progeny cells contained 12-14x fewer mitochondrial genomes compared to their parent NSCs. This drop in mtDNA copy number suggested that neurogenesis might create a somatic mtDNA bottleneck, a sharp reduction in the number of mitochondrial genomes transmitted from a dividing cell to its progeny, akin to the well-documented germline mtDNA bottleneck during oogenesis.

From the imaging data, we observed a consistent downward shift in mean heteroplasmy levels, and increased cell-to-cell variability between sibling progeny compared to their parent NSCs. However, the imaging couldn't resolve if both the observations were the consequence of the same underlying process.

A neutral bottleneck where NSC mtDNA molecules are randomly partitioned into smaller progeny pools can lead to the increased variability in the heteroplasmy levels between progeny cells. But was random inheritance sufficient to explain the observed shift in mean heteroplasmy, or was something else shaping which mitochondrial genomes were retained in the progeny cells?

To test this, we simulated random mtDNA segregation during cell division. For each NSC lineage, using the experimentally observed starting mtDNA copy number and heteroplasmy level, we simulated what the progeny heteroplasmy would look like if the mtDNA molecules were inherited at random. This gave us a lineage-specific null model against which we could compare our experimentally observed metrics.

The tests revealed that while random segregation was sufficient to explain the observed variance in progeny cells, it could not explain the lower proportion of mutant mtDNA in the progeny. This suggested that the heteroplasmy dynamics during neurogenesis were shaped by another process in addition to the developmental bottleneck.

To test this formally, we used Approximate Bayesian Computation (ABC), a simulation-based framework for comparing models and estimating biological parameters when the biological process can be simulated directly but the likelihood is difficult to write down. We simulated progeny heteroplasmy distributions under two models: random segregation alone, and random segregation combined with purifying selection against the mutant mtDNA variant.

The selection model consistently outperformed the neutral model across all lineages. ABC gave us more than a model choice: it allowed us to estimate the strength of selection during this developmental transition. The inferred selection coefficient suggested that mutant mtDNA was roughly half as likely as wild-type mtDNA to be represented in progeny cells.

Simulations in this study therefore did not just rule out the simpler explanation. They converted a qualitative observation, that progeny cells systematically inherited a lower proportion of mutant mtDNA, into a quantitative biological claim about selection operating during neurogenesis.

### From measurement to inference

---

Bioinformatics often gets described as a pipeline: data goes in, results come out. But that data is just a measurement, and measurements can show us that a pattern exists, but not necessarily why or what process produced it.

That is where simulations become part of the scientific argument. They convert assumptions into testable logic by asking what the data should look like if only the simplest processes were operating. The value of simplifying the system is not that the simplified model is complete, but that it tells us what the simplest explanation can and cannot account for.

Together, these studies show why simulations are not just a technical add-on to biological experiments, but an essential part of biological inference. In one case, the simpler explanation was enough; in the other, it was incomplete. That contrast is the point. Simulations helped decide how much biological explanation the data required, making the final claims more cautious in one study and stronger in the other. While experiments reveal the pattern, simulations reveal what kind of process could have produced it.

### Papers and analysis reports

---

- [Analysis Report](https://jvdalab.github.io/mitoPerturb-Seq/MitoPerturbSeq_heteroplasmy_analysis.html) for the MitoPerturb-Seq study (Published in [Nature Structural and Molecular Biology](https://www.nature.com/articles/s41594-026-01779-7)).
- [Analysis Report](https://jvdalab.github.io/mtDNA-smFISH/smFISH_heteroplasmy_analysis.html) for the mtDNA-smFISH study (Published in [Developmental Cell](<https://www.cell.com/developmental-cell/fulltext/S1534-5807(26)00123-1>)).
