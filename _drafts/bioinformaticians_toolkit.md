---
layout: post
title: "A Bioinformatician's Toolkit"
date: 2025-07-20 10:01:00
description: "Essential tools for a bioinformatician"
tags: bioinformatics software-tools
categories: bioinformatics
thumbnail: assets/img/posts/bioinformaticians-toolkit/bioinformaticians-toolkit-thumbnail.webp
giscus_comments: true
disable_animation: true
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

.tool-icon[alt="GitHub logo"] {
  filter: none;
}

html[data-theme="dark"] .tool-icon[alt="GitHub logo"] {
  filter: invert(1);
}
</style>

<div class="row justify-content-center mt-3">
    <div class="col-12 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/bioinformaticians-toolkit/bioinformaticians-toolkit-cover.webp" class="img-fluid rounded z-depth-1" zoomable=false %}
    </div>
</div>

A bioinformatician's daily workflow often involves wrangling large datasets, writing and debugging analysis scripts, managing development environments, and documenting analytical decisions and pipelines for reproducibility.

### Integrated Development Environments (IDEs)

---

Writing scripts and building pipelines for data analysis is the bread and butter of bioinformatics and the humble text editor is often the first tool we reach for. While the combination of a simple text editor and a terminal is all that a bioinformatician needs for their work, the complex nature of today's analyses greatly benefits from the features provided by modern Integrated Development Environments (IDEs).

If you've started your programming journey with R or Python, you may already have come across RStudio or PyCharm. These IDEs provide a rich set of features but are often limited to a single programming language. That's where Visual Studio Code (VSCode) shines.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/vscode.svg" class="tool-icon" alt="VScode logo">
    <a href="https://code.visualstudio.com/" target="_blank" rel="noopener">Visual Studio Code</a>
</span>

Visual Studio Code is the swiss army knife of IDEs. Its rich ecosystem of extensions supports multiple programming languages, integrates with version control systems and terminals and offers powerful code editing and navigation features. The modularity of the extension system also means that you can tailor the IDE to your specific needs.

> **💡 Pro Tip**
>
> Some of my most used VSCode extensions include:
>
> **[Remote Server](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh)**: Connect to HPC clusters via SSH and work directly on remote files. Makes project management on remote systems seamless.
>
> **[R Extension](https://marketplace.visualstudio.com/items?itemName=REditorSupport.r)**: Rich R language support including syntax highlighting, code snippets and R Markdown support.
>
> **[Python Extension](https://marketplace.visualstudio.com/items?itemName=ms-python.python)**: Full Python language support with linting, debugging and Jupyter Notebook integration.
>
> **[Quarto](https://marketplace.visualstudio.com/items?itemName=quarto.quarto)**: Quarto language support for authoring documents and reports.
>
> **[Claude Code for VSCode](https://marketplace.visualstudio.com/items?itemName=anthropic.claude-code)**: Integrate the Claude AI assistant directly into VSCode for code generation and assistance.

### Version Control

---

Bioinformatics analyses always start with data exploration and involve a lot of iteration to refine methods and parameters. Keeping track of changes to code and data files as the analysis evolves is critical.

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

### Lab Notebook (and data syncing)

---

Documenting ideas, observations and critically the decisions made during analysis is an essential part of bioinformatic workflows.

My daily notes are a mix of text blocks, to-do lists, code snippets and links to papers. Working on multiple projects simultaneously adds to the complexity of managing documentation.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/logseq.svg" class="tool-icon" alt="Logseq logo">
    <a href="https://logseq.com/" target="_blank" rel="noopener">Logseq</a>
</span>

Logseq is a powerful knowledge management and note-taking tool that uses a local folder of plain text Markdown files. It supports bi-directional linking, making it easy to connect related notes and ideas. Its outliner format is perfect for organizing thoughts hierarchically, and the ability to embed code blocks allows me to keep track of analysis snippets directly within my notes.

> Logseq has been developing a more powerful graph database backend called Logseq DB. It's still in testing and drops some conveniences associated with the markdown file-based approach, but also enables better performance and new features.

To keep my notes synced across my devices (3 laptops and 1 phone), I use Syncthing.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/syncthing.svg" class="tool-icon" alt="Syncthing logo">
    <a href="https://syncthing.net/" target="_blank" rel="noopener">Syncthing</a>
</span>

### Project Documentation

---

Bioinformatics analyses are often dynamic and iterative, involving multiple tests and adjustments where decisions need to be documented for future reference. Clear documentation is essential for reproducibility and collaboration.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/quarto.svg" class="tool-icon" alt="Quarto logo">
    <a href="https://quarto.org/" target="_blank" rel="noopener">Quarto</a>
</span>

I loved the flexibility of R Markdown and Jupyter Notebooks to document my analyses and produce beautiful reproducible reports. However, Quarto has taken this to the next level by allowing seamless integration of multiple programming languages and output formats. It’s now my go-to tool for project documentation.

### Package Management

---

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

### Containerization

---

If you've ever uttered the phrase "but it works on my computer" when sharing code with a colleague, then you know the pain of environment inconsistencies. Containerization tools allow you to package your scripts along with all their dependencies into a single, portable unit that can run consistently across different systems.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/docker.svg" class="tool-icon" alt="Docker logo">
    <a href="https://www.docker.com/" target="_blank" rel="noopener">Docker</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/apptainer.svg" class="tool-icon" alt="Apptainer logo">
    <a href="https://apptainer.org/docs/user/main/introduction.html" target="blank" rel="noopener">Apptainer</a>
</span>

### Workflow Management

---

A bioinformatic pipeline used to be a bash script that chained together a series of commands. While this works for simple tasks, complex analyses with multiple steps, dependencies, and branching logic benefit from workflow management systems.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/nextflow.svg" class="tool-icon" alt="Nextflow logo">
    <a href="https://www.nextflow.io/" target="_blank" rel="noopener">Nextflow</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/snakemake.svg" class="tool-icon" alt="Snakemake logo">
    <a href="https://snakemake.github.io/" target="_blank" rel="noopener">Snakemake</a>
</span>

### Reference Management

---

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/zotero.svg" class="tool-icon" alt="Zotero logo">
    <a href="https://www.zotero.org/" target="_blank" rel="noopener">Zotero</a>
</span>

### Why These Tools?

---

While each of these tools is a powerful solution on its own, they truly shine when integrated into a cohesive workflow. The IDE along with its extensions streamlines code development and debugging. Version control with Git, GitHub, and DVC ensures that every change to code and data is tracked and recoverable. Logseq and Syncthing keep my notes organized and accessible across devices. Quarto allows me to document my analyses in a reproducible manner. Homebrew and Anaconda simplify package management, while Docker and Apptainer ensure consistent environments. Finally, Nextflow and Snakemake help me build robust, scalable pipelines.

These tools have become indispensable in my bioinformatics workflow because they address the unique challenges of managing complex analyses, large datasets, and collaborative projects. They enhance reproducibility, streamline development, and facilitate effective communication of results.

A firm believer in the spirit of open source, wherever possible I try to incorporate open-source tools into my workflow. This not only supports the broader scientific community but also ensures that my work remains accessible and reproducible.

### What's in Your Toolkit?

---

What tools do you rely on in your bioinformatics work? Share your favorites in the comments!
