---
layout: post
title: "A Bioinformatician's Toolkit"
date: 2025-05-20 10:01:00
description: "Essential tools for a bioinformatician"
tags: bioinformatics
categories: bioinformatics
thumbnail:
giscus_comments: true
---

<style>
.tool-icon {
  height: 45px;
  width: 45px;
  vertical-align: text-bottom;
  margin-right: 8px;
}
.tool {
  display: inline-flex;
  align-items: center;
  gap: 8px;
  font-weight: 600;
  font-size: 1.3rem;
}
/* Invert specific logos in dark mode */
.tool-icon[alt="GitHub logo"] {
  filter: var(--github-logo-filter, none);
}

/* Apply inversion in dark mode */
html[data-theme="dark"] .tool-icon[alt="GitHub logo"] {
  --github-logo-filter: invert(1);
}
</style>

A bioinformatician's daily workflow often involves wrangling large datasets, writing and debugging analysis scripts, managing development environments, and documenting analytical decisions and pipelines for reproducibility.

## Integrated Development Environments (IDEs)

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/vscode.svg" class="tool-icon" alt="VScode logo">
    <a href="https://code.visualstudio.com/" target="_blank" rel="noopener">Visual Studio Code</a>
</span>

Having tried simple text editors like Sublime Text and Atom before, the switch to a dedicated IDE like Visual Studio Code greatly enhanced my daily workflows. The integrated terminal, git support and extensive plugin system make it a powerhouse for routine bioinformatics tasks.

## Project Documentation

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/quarto.svg" class="tool-icon" alt="Quarto logo">
    <a href="https://quarto.org/" target="_blank" rel="noopener">Quarto</a>
</span>

I loved the flexibility of R Markdown and Jupyter Notebooks to document my analyses and produce beautiful reproducible reports. However, Quarto has taken this to the next level by allowing seamless integration of multiple programming languages and output formats. It’s now my go-to tool for project documentation.

## Version Control

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/git.svg" class="tool-icon" alt="Git logo">
    <a href="https://git-scm.com/" target="_blank" rel="noopener">Git</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/github.svg" class="tool-icon" alt="GitHub logo">
    <a href="https://github.com/" target="_blank" rel="noopener">GitHub</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/dvc.svg" class="tool-icon" alt="Data Version Control logo">
    <a href="https://dvc.org/" target="_blank" rel="noopener">Data Version Control (DVC)</a>
</span>

Git, GitHub and Data Version Control (DVC) have been life-savers.

## Lab Notebook (and data syncing)

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/logseq.svg" class="tool-icon" alt="Logseq logo">
    <a href="https://logseq.com/" target="_blank" rel="noopener">Logseq</a>
</span>

To keep my notes synced across my devices (3 laptops and 1 phone), I use Syncthing.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/syncthing.svg" class="tool-icon" alt="Syncthing logo">
    <a href="https://syncthing.net/" target="_blank" rel="noopener">Syncthing</a>
</span>

Lorem

## Package Management

"Dependency hell" is a term that will resonate with any bioinformatician piping together multiple tools for their analysis. Troubleshooting package conflicts and managing different software versions can become a major timesink.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/homebrew.svg" class="tool-icon" alt="Homebrew logo">
    <a href="https://brew.sh/" target="_blank" rel="noopener">Homebrew</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/anaconda.svg" class="tool-icon" alt="Anaconda logo">
    <a href="https://anaconda.org/" target="_blank" rel="noopener">Anaconda</a>
</span>

> **💡 Pro Tip**
>
> Create a separate Conda environment for each project. This keeps dependencies isolated and makes it easier to reproduce analyses months or years later. Export your environment with <code>conda env export > environment.yml</code> and commit it to version control.

## Containerization

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/docker.svg" class="tool-icon" alt="Docker logo">
    <a href="https://www.docker.com/" target="_blank" rel="noopener">Docker</a>
</span>

## Workflow Management

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/nextflow.svg" class="tool-icon" alt="Nextflow logo">
    <a href="https://www.nextflow.io/" target="_blank" rel="noopener">Nextflow</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/snakemake.svg" class="tool-icon" alt="Snakemake logo">
    <a href="https://snakemake.github.io/" target="_blank" rel="noopener">Snakemake</a>
</span>

## Reference Management

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/zotero.svg" class="tool-icon" alt="Zotero logo">
    <a href="https://www.zotero.org/" target="_blank" rel="noopener">Zotero</a>
</span>
