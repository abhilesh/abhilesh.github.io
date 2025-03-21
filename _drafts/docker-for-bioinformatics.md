---
layout: post
title: "Docker for Bioinformatics"
date: 2025-02-27 21:01:00
description: "Portable, Scalable and Reproducible Bioinformatics workflows"
tags: docker bioinformatics
categories: bioinformatics
thumbnail: assets/img/posts/docker-for-bioinformatics/docker_for_bioinformatics_thumbnail.png
giscus_comments: true
tabs: true
toc:
  sidebar: left
---

<div class="row justify-content-center mt-3">
    <div class="col-12 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/docker-for-bioinformatics/docker_for_bioinformatics_cover.png" class="img-fluid rounded z-depth-1" zoomable=false %}
    </div>
</div>

Bioinformatics analysis often involves complex pipelines with rapidly evolving software tools, each with their own set of dependencies. System compatibility, version mismatches and dependency conflict issues can often be a nightmare, making running and sharing bioinformatic pipelines a challenging task. These challenges not only waste valuable research time but also contribute to irreproducible workflows, where results depend as much on the computing environment as on the analysis itself. Docker offers a powerful solution by packaging software and its dependencies into portable, reproducible containers—ensuring that your bioinformatics pipelines run consistently, whether on your local machine, an HPC cluster, or the cloud.

### What is Docker?

Imagine you're baking a cake, but every time you try, your kitchen is missing key ingredients or uses a different oven that bakes at the wrong temperature. Docker is like a self-contained baking kit that comes with all the right ingredients, tools, and even its own portable oven, ensuring your cake turns out exactly the same no matter where you bake it. In bioinformatics, Docker does the same for software—packaging tools, dependencies, and environments so your analyses run reliably, whether on your laptop, an HPC cluster, or the cloud.

### How Docker can help?

The **FAIR** (Findable, Accessible, Interoperable, and Reusable) principles provide practical guidelines for maximizing the value of research data by ensuring it is discoverable, accessible, compatible across systems, and usable in future studies. Docker plays a key role in aligning bioinformatics workflows with these principles:

**_Benefits of Using Docker in Bioinformatics (FAIR-aligned)_**

- **Findability**: Docker images can be easily found on registries like Docker Hub, with tags and descriptions facilitating discovery. This allows researchers to locate the specific software versions they need.

- **Accessibility**: Docker images are generally accessible to anyone with an internet connection and Docker installed, promoting open science. Tools and their dependencies are packaged, reducing access barriers caused by complex installations.

- **Interoperability**: Docker promotes interoperability by providing a consistent environment across different systems. Containers encapsulate all dependencies, minimizing conflicts and ensuring that tools function as expected in various computing environments.

- **Reusability**: Docker images promote reusability of bioinformatics tools and workflows. Researchers can easily reuse existing images, knowing that the software environment is well-defined and consistent. This fosters collaboration and reduces redundant effort.

### Getting Started

To begin, start by installing Docker on your system. Docker is available for all major operating systems and the installers can be downloaded from the [official website](https://www.docker.com/get-started/). The Docker Desktop application is recommended for Windows and macOS users, while Linux users can install Docker using the command line.

{% tabs docker-os-install %}

{% tab docker-os-install MacOS %}

- Download the [Docker Desktop installer for MacOS](https://docs.docker.com/desktop/setup/install/mac-install/)
- Follow the setup instructions to install Docker.
- Docker Desktop can now be launched by clicking the application icon.

{% endtab %}

{% tab docker-os-install Windows %}

- Download the [Docker Desktop installer for Windows](https://docs.docker.com/desktop/setup/install/windows-install/)
- Follow the setup instructions to install Docker.
- Docker Desktop can now be launched by clicking the application icon.

{% endtab %}

{% tab docker-os-install Linux %}

Assuming a Debian-based Linux distribution (e.g., Ubuntu):

- Open a terminal and run the following commands

  ```bash
  # Update package index and install docker
  sudo apt-get update && sudo apt-get install docker.io
  ```

- Add your user to the `docker` group to run Docker commands without `sudo`

  ```bash
  sudo usermod -aG docker ${whoami}
  ```

- Log out and log back in to apply the changes

> _Note:_ If you do not have `sudo` privileges, you can install Docker using the [official script](https://docs.docker.com/engine/install/ubuntu/#install-using-the-convenience-script).

```bash
# Download the Docker installation script
curl -fsSL https://get.docker.com -o get-docker.sh

# Append --dry-run to see the commands that will be executed
sudo sh get-docker.sh
```

{% endtab %}

{% endtabs %}

To test whether Docker is installed correctly, run the following command in your terminal:

```bash
# Check Docker version
docker --version

# Test with a simple hello-world container
docker run hello-world
```

There are other ways to run containers, and you can experiment with these as you get more comfortable with Docker.

### A Simple Docker Example

1. Pull a Docker Image

   Here, we will pull a Docker image for the popular bioinformatics tool [`samtools`](https://hub.docker.com/r/biocontainers/samtools).

   ```bash
   # Pull the samtools image from Docker Hub
   docker pull biocontainers/samtools
   ```

2. Run the Docker Container

   ```bash
   docker run -d -name samtools -v /path/to/data:/data biocontainers/samtools tail -f /dev/null
   ```

### Finding Bioinformatics Tool Containers

These registries host a large number of pre-built Docker images for bioinformatics tools:

- <img src="/assets/img/posts/docker-for-bioinformatics/docker-4.svg" width="20" height="20" style="margin-right: 5px;"> [Docker Hub](https://hub.docker.com/)

- <img src="/assets/img/posts/docker-for-bioinformatics/biocontainers-logo.svg" width="20" height="20" style="margin-right: 5px;"> [Biocontainers](https://biocontainers.pro/)

- <img src="/assets/img/posts/docker-for-bioinformatics/quayio-logo.svg" width="20" height="20" style="margin-right: 5px;"> [Quay.io](https://quay.io/organization/biocontainers)

- <img src="/assets/img/posts/docker-for-bioinformatics/pegi3s-logo.svg" width="20" height="20" style="margin-right: 5px;"> [pegi3s Bioinformatics Docker Images Project](http://bdip.i3s.up.pt/)

### Other Resources

If you'd like to read more about reproducible research practices, check out the following resources:

- ["The Turing Way"](https://book.the-turing-way.org/) : A handbook for reporoducible, ethical and collaborative data science.
