---
layout: distill
title: "Version Control for Bioinformatics"
date: 2025-05-15 21:01:00
description: "Tracking changes, collaborating and ensuring reproducibility"
tags: git bioinformatics
categories: bioinformatics
thumbnail:
giscus_comments: true
#disqus_comments: true
tabs: true
disable_animation: true
toc: true
authors:
  - name: Abhilesh Dhawanjewar
    url: "https://abhilesh.github.io"
    affiliations:
      name: University of Cambridge
      url: "https://www.cam.ac.uk/"
bibliography: 2025-02-27-docker-for-bioinformatics.bib
---

<div class="row justify-content-center mt-3">
    <div class="col-12 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/version-control-for-bioinformatics/version-control-for-bioinformatics-cover.jpg" class="img-fluid rounded z-depth-1" zoomable=false %}
    </div>
</div>

Bioinformatics workflows are often iterative, involving multiple steps of data processing, analysis, and visualization. As these workflows evolve, it becomes crucial to maintain a record of changes, collaborate with others, and ensure reproducibility. Version control systems (VCS) like Git provide a robust solution to these challenges.

## Motivation

In computational biology, our analyses are driven by two key components: the **code** (scripts, pipelines) and the **data** (FASTQ files, reference genomes, alignment maps, statistical outputs). Achieving full reproducibility—a cornerstone of scientific inquiry—requires us to version _both_.

For code, the solution is established: **Git** is the industry standard for tracking changes to our scripts, R markdowns, and pipeline definitions. However, Git was not designed to handle the multi-gigabyte files (BAMs, CRAMs, or Hi-C contact maps) that are routine in our field. Committing a 50GB alignment file to a Git repository is not just impractical; it's impossible.

This leads to a common, fragmented workflow: code is tracked in Git, while large data files are managed manually on a server, an S3 bucket, or an external hard drive, often with filenames like `analysis_v2_final_REVISED.csv`. This "data gap" breaks the link between the code and the data it processed, making true reproducibility a significant challenge.

How can we unite our code and data under a single, cohesive versioning system? This post outlines my workflow using **Git** for code and **DVC (Data Version Control)** for data.

---

### The Base Layer: Git for Code

Git is non-negotiable for tracking our analytical code. It provides a complete history of our project, allowing us to:

- Track every change to our analysis scripts (Python, R, shell).
- Collaborate effectively with colleagues by merging changes.
- Create branches to test new ideas (e.g., trying a different aligner) without breaking the main analysis.
- Tag specific versions (commits) of our code to correspond with a paper submission or a specific figure.

A typical Git workflow is simple:

```bash
# Initialize a new repository
git init
# Add a script to be tracked
git add scripts/my_analysis.R
# Save a snapshot (commit) of the project
git commit -m "Add initial analysis script for differential expression"
# Send it to a remote backup (like GitHub or GitLab)
git push
```

This is perfect for text-based files. The problem begins when we introduce our data.

### The Solution: DVC for Large Data

DVC (Data Version Control) is an open-source tool designed to handle large files, models, and datasets. It works alongside Git, allowing Git to do what it does best (versioning code) while DVC handles what it does best (versioning data).

DVC does not store the data in the Git repository. Instead, it creates small "pointer" files (or metafiles) that end in .dvc. These .dvc files are lightweight, text-based, and contain the information needed to retrieve the correct version of the data, such as an MD5 hash.

Here’s the key concept:

You add your large data file (e.g., data/raw/sample1.fastq.gz) to DVC.

DVC moves the file to a central cache (usually outside your project directory) and creates a small pointer file: data/raw/sample1.fastq.gz.dvc.

You add this small .dvc pointer file to Git.

You configure DVC to push the actual data (from its cache) to a remote storage location, such as an S3 bucket, Google Cloud Storage, or an SSH server.

Your Git repository now tracks the lightweight .dvc files, which serve as a manifest of the data required for the analysis. Your colleague can git pull your code, and then run a single dvc pull command to download all the correct data versions referenced by the .dvc files.

### A Practical Bioinformatics Workflow: Git + DVC

Here is a complete, simplified workflow for a hypothetical project.

1. Setup (One-time only):

```bash
# Create the project
mkdir diff-exp-project
cd diff-exp-project

# Initialize Git (for code)
git init

# Initialize DVC (for data)
dvc init
# Git now starts tracking DVC's internal config files
git commit -m "Initialize DVC"

# Tell DVC where to store the actual data
# (This could be S3, GCP, or a server you SSH into)
dvc remote add -d my-remote /path/to/my/server-storage
```

2. Adding Data and Code: Let's say we have a large reference genome and our analysis script.

```bash
# Add the large genome file to DVC
dvc add data/reference.fasta

# Add our analysis script to Git
git add scripts/run_alignment.sh

# Now, look at what Git sees:
git status
```

Git will not see data/reference.fasta. It will only see data/reference.fasta.dvc (the small pointer file) and scripts/run_alignment.sh.

3. Committing the Project State: We commit this project snapshot, which now includes both the code and the pointers to the data.

```bash
# Add the .dvc file and our script to Git
git add data/reference.fasta.dvc scripts/run_alignment.sh
git commit -m "Add reference genome and alignment script"

# Push our code and data (in two steps)
# 1. Push code (and pointers) to GitHub
git push origin main

# 2. Push data (via DVC) to our server storage
dvc push
```

Our project is now fully versioned and backed up.

4. Collaborating (The Payoff): A new lab member wants to reproduce your analysis.

```bash
# They clone the Git repository (gets code + pointers)
git clone [https://github.com/my-lab/diff-exp-project.git](https://github.com/my-lab/diff-exp-project.git)
cd diff-exp-project

# They pull the data using DVC
# DVC reads the .dvc files and fetches the
# correct data versions from your server storage
dvc pull

# That's it. They now have the exact code and data
# from your commit and can run the analysis.
bash scripts/run_alignment.sh
```

### Conclusion: Towards Fully Reproducible Science

By combining Git and DVC, we create a single, unified workflow that versions every component of our research. The Git repository becomes the complete "source of truth," tracking our code and providing a manifest of all associated data.

This approach provides a robust framework for:

- Reproducibility: Anyone (including your future self) can check out a specific branch or commit and retrieve the exact code and data used to generate a result.

- Collaboration: Onboarding new team members is as simple as git pull and dvc pull.

- Efficiency: We stop manually tracking large files and let our tools handle the versioning, allowing us to focus on the analysis itself.
