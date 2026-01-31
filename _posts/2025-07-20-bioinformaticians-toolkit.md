---
layout: post
title: "A Bioinformatician's Toolkit"
date: 2025-07-20 10:01:00
description: "Tools that power my daily bioinformatics workflow"
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

A bioinformatician's daily workflow is rarely a straight line. It's a messy loop that often involves wrangling large datasets, writing and debugging analysis scripts, managing development environments, and documenting analytical decisions and pipelines for reproducibility.

Over the years, I've assembled a toolkit of software that has dramatically improved both my productivity and the reproducibility of my work. These tools have been an essential infrastructure of my daily workflow, helping keep things organized and efficient.

### Integrated Development Environments (IDEs)

---

Writing scripts and building pipelines for data analysis is the bread and butter of bioinformatics and the humble text editor is often the first tool we reach for. While the combination of a simple text editor and a terminal is all that a bioinformatician needs for their work, the complex nature of today's analyses greatly benefits from the features provided by modern Integrated Development Environments (IDEs).

If you've started your programming journey with R or Python, you may already have come across [RStudio](https://posit.co/download/rstudio-desktop/) or [PyCharm](https://www.jetbrains.com/pycharm/). These IDEs provide a rich set of features but are often limited to a single programming language. That's where Visual Studio Code (VSCode) shines.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/vscode.svg" class="tool-icon" alt="VScode logo">
    <a href="https://code.visualstudio.com/" target="_blank" rel="noopener">Visual Studio Code</a>
</span>

Visual Studio Code truly is the swiss army knife of IDEs. Its rich ecosystem of extensions supports multiple programming languages, integrates with version control systems and terminals and offers powerful code editing and navigation features. The modularity of the extension system also means that you can tailor the IDE to your specific needs.

> **💡 Pro Tip** <br><br>Some of my most used VSCode extensions include:<br><br>
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

### AI-Assisted Programming

---

The transformation that AI-assisted tools have had in our lives needs no introduction. Agentic Large Language Models (LLMs) have significantly altered the way we code, debug and learn new programming concepts.

<span class="tool">
<img src="/assets/img/posts/bioinformaticians-toolkit/claude.svg" class="tool-icon" alt="Claude logo">
<a href="https://claude.ai/" target="_blank" rel="noopener">Claude</a>
</span>

Anthropic's Claude Code excels at working with large and messy codebases. I use it to understand old code, refactor code and generate documentation. It's like having a patient code reviewer when debugging or revisiting code I wrote months ago.

<span class="tool">
<img src="/assets/img/posts/bioinformaticians-toolkit/gemini-color.svg" class="tool-icon" alt="Google Gemini logo">
<a href="https://gemini.google.com/" target="_blank" rel="noopener">Google Gemini</a>
</span>

Google's Gemini is great at exploring alternative ideas, translating code between languages and general prototyping of analyses.

Relying on these tools has replaced the satisfaction of solving a complex logic puzzle with the more passive role of a code reviewer, but I've found that they can be powerful allies in streamlining the mundane bits. Although the accuracy and reliability of these models has improved significantly, they still make occasional missteps that need careful review. I use them primarily to find code snippets across large projects, add documentation to existing codebases and suggest improvements to my code. As this landscape is rapidly evolving, I often combine multiple LLMs to cross-validate suggestions and ensure accuracy.

> **💡 Pro Tip** <br><br> Use the <b>"Project Instructions"</b> section (typically defined via a markdown file) to specify the desired behavior from the LLMs. Incorporating specific instructions like <i>"Be critical and objective; Refactor the code minimally preserving logic; provide short answers unless specifically asked"</i> help refine the LLM's outputs<br><br> Integrating these tools directly with VSCode via extensions ([Claude](https://marketplace.visualstudio.com/items?itemName=anthropic.claude-code) and [Gemini](https://marketplace.visualstudio.com/items?itemName=Google.geminicodeassist)), makes it extremely powerful and convenient due to the ability to share the codebase and get context-aware suggestions without leaving the IDE.

### Version Control

---

Bioinformatics analyses always start with data exploration and involve a lot of iteration to refine methods and parameters. Keeping track of changes to code and data files as the analysis evolves is critical and inadvertently saves the day when revisiting old analyses. Version control systems are a powerful time machine that enable me to live by the wise adage: _"Be kind to your future self"_. If your scripts are named `final_analysis_v2_final_really_final.R`, you need version control.

<span class="tool">
<img src="/assets/img/posts/bioinformaticians-toolkit/git.svg" class="tool-icon" alt="Git logo">
<a href="https://git-scm.com/" target="_blank" rel="noopener">Git</a>
</span>

Git is the backbone of modern version control. It allows you to track changes to your codebase and maintain multiple simultaneous versions (branches) for experimentation. Think of it like a sophisticated "undo" button that works across your entire history.

I now initialize a Git repository for every project I start and commit changes frequently with meaningful messages. This encourages me to follow a deliberate rhythm of coding, testing and documenting changes while supercharging my ability to "time travel" and trace back to previous versions of my code when needed.

> **💡 Pro Tip** <br><br> `git` is natively supported in VSCode, but the [Git Graph](https://marketplace.visualstudio.com/items?itemName=mhutchie.git-graph) extension is an easy way to familiarize yourself with version control concepts by providing a visual interface to manage commits, branches and merges.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/github.svg" class="tool-icon" alt="GitHub logo">
    <a href="https://github.com/" target="_blank" rel="noopener">GitHub</a>
</span>

GitHub extends Git's capabilities by providing a cloud-based collaboration platform making the coding process more social and manageable. Beyond just storing your code remotely, it comes integrated with a suite of features that make managing collaborative projects easier.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/dvc.svg" class="tool-icon" alt="Data Version Control logo">
    <a href="https://dvc.org/" target="_blank" rel="noopener">Data Version Control (DVC)</a>
</span>

While Git excels at tracking code, bioinformatics workflows involve large data and intermediate files that can quickly balloon repository sizes into the gigabytes or terabytes. We often need to version control not just our scripts, but also datasets, model files, and results. Data Version Control (DVC) is a powerful extension to Git that allows the tracking of large files, datasets, and analysis versions. DVC works by storing metadata about large files in Git enabling the precise versioning of data alongside code.

### Lab Notebook (and syncing notes)

---

Documenting ideas, observations and critically the decisions made during analysis is an essential part of bioinformatic workflows. Add to that the notes from meetings, literature summaries and to-do lists, and a mountain of unorganized notes accumulates in no time. While a traditional pen-and-paper note taking system has its charms, I realized early on that effective recall was the bottleneck in my note-taking process. I still see so many colleagues using MS Word to take notes, a tool ill-suited for this purpose, while several modern alternatives exist.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/logseq.svg" class="tool-icon" alt="Logseq logo">
    <a href="https://logseq.com/" target="_blank" rel="noopener">Logseq</a>
</span>

Logseq is a powerful knowledge management and note-taking tool that uses a local folder of plain text Markdown files. It supports bi-directional linking, making it easy to connect related notes and ideas. Its outliner format is perfect for organizing thoughts hierarchically, and the ability to embed code blocks allows me to keep track of analysis snippets directly within my notes.

> **💡 Pro Tip** <br><br>Logseq is an outliner at its core, meaning that every note is structured as a series of nested bullet points. When you first open Logseq, you are greated with a _daily journal_ page. For me, this reduces the friction of having to think where to put my notes. Everything is jotted on the daily page with the relevant projects and topics linked as needed. Over time, a web of interconnected notes builds up that is easily searchable. <br><br>If you prefer a folder-based organization system that plays well with long-form writing, [Obsidian](https://obsidian.md/) is another excellent alternative with similar features.

Often times, I am working across devices: laptop, home server and a mobile device. Having access to my notes wherever I am has made data recall a seamless process.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/syncthing.svg" class="tool-icon" alt="Syncthing logo">
    <a href="https://syncthing.net/" target="_blank" rel="noopener">Syncthing</a>
</span>

Syncthing is an amazing continuous file synchronization program that keeps my files synced across all my devices without relying on third-party cloud services and paying costly subscription fees. Paired with Logseq, it keeps all my notes up to date no matter which device I'm using.

### Project Documentation

---

Bioinformatic analyses are not just a series of steps to be executed, but rather a chain of decisions that need to be clearly documented. When sharing results with collaborators or preparing a manuscript, having a well-documented analysis workflow is a lifesaver. "Literate programming", a term coined by Donald Knuth, emphasizes the importance of writing code that is understandable by humans, by interweaving code with narrative text.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/quarto.svg" class="tool-icon" alt="Quarto logo">
    <a href="https://quarto.org/" target="_blank" rel="noopener">Quarto</a>
</span>

I loved the flexibility of R Markdown and Jupyter Notebooks to document my analyses and produce beautiful reproducible reports. However, Quarto has taken this to the next level by allowing seamless integration of multiple programming languages and output formats. It’s now my go-to tool for project documentation.

### Package Management

---

"Dependency hell" is a rite of passage for any bioinformatician piping together multiple tools for their analysis. Troubleshooting package conflicts and managing different software versions usually becomes a massive timesink.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/homebrew.svg" class="tool-icon" alt="Homebrew logo">
    <a href="https://brew.sh/" target="_blank" rel="noopener">Homebrew</a>
</span>

Homebrew is a package manager for macOS and Linux, operating systems that I've used predominantly. It's become my go-to package manager ever since I switched to macOS over a decade ago. It makes installing and managing command-line tools and libraries delightfully simple

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/anaconda.svg" class="tool-icon" alt="Anaconda logo">
    <a href="https://anaconda.org/" target="_blank" rel="noopener">Anaconda</a>
</span>

Anaconda is a powerful environment manager for Python and R that allows the creation of isolated environments with specific package versions. This is particularly useful when working on multiple projects that often require specific library versions.

> **💡 Pro Tip** <br><br> For environments where system-level access is restricted such as HPC clusters, Anaconda is a really useful tool to manage project-specific dependencies enabling the use of up-to-date versions as needed.

Anaconda has often been blamed for bloated environments and slower performance. For those concerned with Anaconda's footprint, [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) offers a lighter base install. Furthermore, other package managers like [`mamba`](https://github.com/mamba-org/mamba) and [`pixi`](https://pixi.sh/latest/) can be used as drop-in replacements for `conda` to speed up package installations. [Pip](https://pypi.org/project/pip/) and [`venv`](https://docs.python.org/3/library/venv.html) are alternative Python-only package and environment management tools.

> Homebrew is a **system-level** package manager, while Anaconda is a **project-level** package manager. While Anaconda is primarily a Python distribution, it also supports R and other languages through its conda package manager.

### Containerization

---

If you've ever uttered the phrase "but it works on my computer" when sharing code with a colleague, then you know the pain of environment inconsistencies. Containerization tools allow you to package your scripts along with all their dependencies into a single, portable unit that can run consistently across different systems.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/docker.svg" class="tool-icon" alt="Docker logo">
    <a href="https://www.docker.com/" target="_blank" rel="noopener">Docker</a>
</span>

Docker quickly became the de facto standard for containerization in bioinformatics. It allows the creation of lightweight, portable containers that ship analysis code along with the entire software environment, enabling reproducibility across different computing environments.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/apptainer.svg" class="tool-icon" alt="Apptainer logo">
    <a href="https://apptainer.org/docs/user/main/introduction.html" target="blank" rel="noopener">Apptainer</a>
</span>

Apptainer is the preferred containerization solution for high-performance computing (HPC) environments where Docker's requirement for root privileges is a limitation. It allows users to run containers in user space, making it suitable for shared computing clusters.

### Workflow Management

---

A bioinformatic pipeline used to be a bash script that chained together a series of commands. While this works for simple tasks, complex analyses with multiple steps, dependencies, and branching logic benefit from workflow management systems.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/snakemake.svg" class="tool-icon" alt="Snakemake logo">
    <a href="https://snakemake.github.io/" target="_blank" rel="noopener">Snakemake</a>
</span>

Snakemake is a popular workflow management system that uses a Python-based syntax to define rules for each step of the analysis. It automatically handles dependencies, parallelization, and resource management, making it easy to scale analyses from a single machine to a cluster or cloud environment.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/nextflow.svg" class="tool-icon" alt="Nextflow logo">
    <a href="https://www.nextflow.io/" target="_blank" rel="noopener">Nextflow</a>
</span>

Nextflow is a modern workflow manager that emphasizes scalability and reproducibility. It uses a dataflow programming model, allowing for complex workflows with dynamic branching and parallel execution. Nextflow integrates seamlessly with containerization tools and cloud platforms.

The real benefit comes when you need to run the same analysis on dozens or hundreds of samples—these tools handle the parallelization and resource management automatically, turning what could be days of manual work into a single command.

### Reference Management

---

Keeping track of the dozens of papers you read for a single project can quickly become overwhelming. Between literature reviews, methods sections, and supporting your analysis decisions with references, a good reference manager becomes essential.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/zotero.svg" class="tool-icon" alt="Zotero logo">
    <a href="https://www.zotero.org/" target="_blank" rel="noopener">Zotero</a>
</span>

Zotero is my reference manager of choice. It works seamlessly across different operating systems, and has a browser extension that makes capturing references from the internet effortless. It also integrates well with word processors for easy citation insertion.

### Why These Tools?

---

While each of these tools is a powerful solution in isolation, they truly shine when woven into a cohesive workflow.

VSCode serves as the central interface for my work. Its Remote SSH integrations bridge the gap between local development and remote HPC clusters, allowing me to edit code on a supercomputer as if it were on my laptop. Simultaneously, Git ensures every change is reliably versioned and backed up. When I encounter an interesting paper or need to document a critical decision, Logseq becomes its digital home. Its bi-directional linking keeps projects organized, and the Zotero integration allows me to annotate papers and link them directly to my notes without breaking context.

Once the analysis is complete, I rely on Quarto to transform raw scripts and messy exploration into polished, reproducible reports that collaborators can actually use. Supporting all of this, package managers and containerization solutions act as the silent infrastructure in the background. They eliminate the "dependency hell" that used to consume days of my time.

The tools listed here might not be the "best" in their class and many alternatives do exist. However, I have found this specific combination to offer the best balance of robust documentation and active community support. These traits make troubleshooting significantly easier when issues inevitably arise. I'm also an ardent supporter of open-source software, and wherever possible I try to incorporate open-source tools into my workflow.

### What's in Your Toolkit?

---

What tools do you rely on in your bioinformatics work? Share your favorites in the comments!
