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

Modern biological experiments can generate an astonishing amount of detail about living systems. Single-cell technologies allow us to profile thousands of individual cells at once, imaging approaches let us visualize molecules inside tissues, and perturbation experiments make it possible to ask how biological systems respond when specific genes or pathways are disrupted. This is especially powerful for mitochondrial biology, where mutations and their effects can vary between cells, tissues, developmental stages, and even diseases.

The extraordinary resolution and detailed measurements afforded by these technologies however, hasn't made their interpretation any easier. The deeper we probe, the more levels of variation we find, and that variation can be difficult to interpret. The observations can show us a pattern, but cannot automatically tell us which part of the pattern is biology, which part is noise, and which part is simply a consequence of the measurement process itself.

This is where simulations become powerful. To interpret an observation, we need to know what we would have observed if nothing interesting had happened: the counterfactual. That reference point rarely comes from the experiment itself. We have to build it and simulations enable us to do that. They make our assumptions explicit, generate the data those assumptions predict, and then ask whether what we actually observed is surprising.

**Experiments show us what happened. Simulations ask what else could have happened.**

### Building the null world

---

The key step in simulation-based inference is building what I think of as a deliberately simplified world, also sometimes known as the null model. In this world, we decide which processes are allowed to operate and which ones are left out. We might allow random sampling, measurement error, or differences in detection sensitivity, while leaving out the more interesting mechanisms we are trying to test. This lets us isolate the different factors that can influence our measurements and understand the kinds of patterns they can produce on their own.

That matters because biological data can be surprisingly good at producing convincing-looking patterns by chance. A difference between two groups might appear because one group was measured more deeply than the other. A rare event might look enriched simply because the sample size was small. A shift in a distribution might reflect random inheritance, uneven detection, or technical variation rather than an active biological process.

Simulations let us make these expectations explicit. We can generate the data we would expect under the simpler world, then compare our observations to it. If the observed pattern falls inside that expectation, the simulation protects us from over-interpreting the result. If the observed pattern falls outside it, the simulation tells us where the simple explanation breaks down.

The comparison is the useful part. As the two studies below show, simulations can sharpen interpretation in either direction. Sometimes, the observed pattern sits comfortably inside the null expectation, and the simulation saves us from turning a plausible result into a wrong conclusion. Other times, the observed pattern falls outside that expectation, and the simulation is what makes the stronger biological claim possible.

### Mitochondrial DNA and heteroplasmy

---

Mitochondria carry their own genome: mitochondrial DNA (mtDNA) which exists in hundreds to thousands of copies per cell. The population of mtDNA in a cell can be heterogeneous, with copies carrying different mutant variants, a phenomenon termed heteroplasmy. In the clinical context, heteroplasmy is the presence of a disease-causing variant alongside the normal, wild-type mtDNA and the proportion of the mutant copies in the cell determines whether and how severely a cell is affected. The same mutation can be silent or damaging in a cell depending on the heteroplasmy level (proportion of mutant mtDNA) of the cell.

This makes mitochondrial genetics unusually quantitative. It is not enough to ask whether a mutation is present or absent. We often need to ask how much of it is present, how the proportion varies between cells, and how it changes as cells divide, differentiate or respond to stress. Because cells contain finite numbers of mitochondrial genomes, and experiments measure only a subset of them, heteroplasmy patterns can shift or vary for reasons that are partly stochastic, partly technical, and partly biological. Simulations are useful because they let us ask whether random sampling, measurement depth, or mtDNA copy number changes are sufficient to explain the observed data before invoking a stronger biological mechanism.

### When the null is the answer

---

The mitochondrial genome is sparse, with the subunits of the electron transport chain, being the only protein-coding genes it encodes. Nuclear genes play a central role in mitochondrial genome biology. They regulate how mtDNA is replicated, packaged, and maintained, but establishing direct causal links between specific genes and heteroplasmy dynamics at the level of individual cells has been difficult.

Our recently published study in Nature Structural and Molecular Biology describes MitoPerturb-Seq, a technique designed to investigate these links directly. It perturbs a set of nuclear genes, then measures mitochondrial DNA copy number, heteroplasmy, chromatin accessibility, and gene expression in thousands of individual cells.

The surprising result was that mean heteroplasmy was unaffected across all these perturbations but cells with _Tfam_, _Opa1_, and _Polg_ knockdown exhibited a significant reduction in mtDNA copy number, along with increased heteroplasmy variance. At first glance, this observation can suggest that disrupting these genes makes mitochondrial genetic variation less stable across cells.

However, the increase in heteroplasmy variance can also be confounded by the reduced mtDNA copy number we observe in cells where these three genes are perturbed. When fewer mtDNA molecules are present, or captured by the assay, heteroplasmy estimates become noisier. Cells can then appear more variable even if the underlying distribution of mutant and wild-type mtDNA has not changed in an unexpected way. The observed pattern can have two possible explanations: altered heteroplasmy dynamics, or the expected consequence of measuring fewer mitochondrial genomes.

This is where the simulations helped assess the pattern.  I built a null model that asked what heteroplasmy variance would look like if the main thing changing was mtDNA depth or what would the heteroplasmy in the control cells look like if they had mtDNA depths similar to the _Tfam_, _Opa1_, and _Polg_ perturbed cells. Starting from the control heteroplasmy distribution, the simulation repeatedly resampled cells at the reduced depths observed in the depleted perturbations. This created the counterfactual we needed: a simplified world where mtDNA depth dropped, but no additional biological process was added to increase heteroplasmy variance.

When we compared the observed variances to that simulated expectation, the apparent increase in heteroplasmy variance was not surprising. The depleted perturbations sat within the range expected from reduced mtDNA depth alone. The simulation therefore showed that the increased variants did not require an additional mechanism beyond mtDNA depletion and the sampling noise that comes with measuring fewer mitochondrial genomes.

The simulation made the interpretation more precise. The perturbations clearly did something biologically meaningful: they depleted mtDNA and changed the cellular response to mitochondrial stress. But the increase in heteroplasmy variance did not need to be interpreted as a separate effect. Once we accounted for reduced mtDNA depth, the pattern was something the simpler model could already explain.

### When the null world isn't enough

---

The MitoPerturb-Seq study measured heteroplasmy indirectly, through sequencing reads from thousands of dissociated cells. The mtDNA-smFISH study instead visualizes mitochondrial DNA molecules inside intact cells and tissues. Using sequence-specific fluorescent probes, mtDNA-smFISH allows wild-type and mutant mitochondrial genomes to be counted directly in whole tissues, while preserving spatial and cellular context. This enabled visualization of heteroplasmy distributions in actively dividing cells such as the neural stem cells.

Applying this technique on developing neural stem cells (NSCs) in a heteroplasmic _Drosophila_ line enabled the visualization and quantification of heteroplasmy dynamics as the NSCs divided to give rise to progeny cells. We observed that the progeny cells contained 12-14x fewer mitochondrial genomes compared to their parent NSCs. This drop in mtDNA copy number suggested that neurogenesis might create a somatic mtDNA bottleneck, a sharp reduction in the number of mitochondrial genomes transmitted from a diving cell to its progeny, akin to the well-documented germline mtDNA bottleneck during oogenesis.

From the imaging data, we observed a consistent downward shift in mean heteroplasmy levels, and increased cell-to-cell variability between sibling progeny compared to their parent NSCs. However, the imaging couldn't resolve if both the observations were the consequence of the same underlying process.

A neutral bottleneck where NSC mtDNA molecules are randomly partitioned into smaller progeny pools can lead to the increased variability in the heteroplasmy levels between progeny cells. But was random inheritance sufficient to explain the observed shift in mean heteroplasmy, or was something else shaping which mitochondrial genomes were retained in the progeny cells?

To test this, we simulated random mtDNA segregation during cell division. For each NSC lineage, using the experimentally observed starting mtDNA copy number and heteroplasmy level, we simulate what the progeny heteroplasmy would look like if the mtDNA molecules were inherited at random. This gives us a lineage-specific null model against which we could compare our experimentally observed metrics.

The tests revealed that while random segregation was sufficient to explain the observed variance in progeny cells, it could not explain the lower proportion of mutant mtDNA in the progeny. This suggested that the heteroplasmy dynamics during neurogenesis were shaped by another process in addition to the developmental bottleneck.

To assess this formally, we extended the simulation model to include a selection parameter, representing the likelihood of the wild-type mtDNA variant to be retained over the mutant variant. We then used Approximate Bayesian Computation (ABC) to compare the two models: one where heteroplasmy changed via random segregation alone, and another where random selection was accompanied by purifying selection against the mutant mtDNA.

Approximate Bayesian Computation (ABC) utilizes a simulation framework to formally compare competing models and estimate parameters that explain the observed data. The ABC tests corroborated the result that the random segregation accompanied by purifying selection explained the observed data better than random segregation alone and also enabled us to quantify the strength of selection against the mutant mtDNA variant.

The simulations thus helped provide quantitative support for the underlying biological processes

The mtDNA-smFISH study was published in Developmental Cell and the full analysis report is available online at - https://jvdalab.github.io/mtDNA-smFISH/smFISH_heteroplasmy_analysis.html
