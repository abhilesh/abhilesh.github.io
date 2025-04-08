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

Imagine you're baking a cake, but every time you try, your kitchen is missing key ingredients or uses a different oven that bakes at the wrong temperature. Docker is like a self-contained baking kit that comes with all the right ingredients, tools, and even its own portable oven, ensuring your cake turns out exactly the same no matter where you bake it. In bioinformatics, Docker does the same for software by packaging tools, dependencies, and environments so that analyses run reliably across different computing platforms.

### How can Docker help?

The **FAIR** (Findable, Accessible, Interoperable, and Reusable) principles provide guidelines for maximizing the value of research data. Docker aligns bioinformatics workflows with these principles by ensuring software and environments are portable and reproducible:

- **Findability**: Docker images can be easily found in registries like Docker Hub, with clear versioning and documentation for discovery.

- **Accessibility**: Anyone with internet access and Docker installed can retrieve and use containerized tools, promoting open science.

- **Interoperability**: Docker containers provide a consistent runtime environment across different systems, preventing dependency conflicts.

- **Reusability**: Pre-built images allow researchers to reuse workflows without worrying about installation issues, fostering collaboration.

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

### Running Bioinformatics Tools with Docker

We will use the popular tool [`samtools`](http://www.htslib.org/) as an example to demonstrate how to run bioinformatics tools using Docker. `samtools` is a widely used tool for working with Sequence Alignment/Map (SAM) and Binary Alignment/Map (BAM) files.

1. **Pull a Docker Image**

   Here, we will pull a Docker image for [`samtools`](https://hub.docker.com/r/biocontainers/samtools) from Docker Hub.

   ```bash
   # Pull the samtools image from Docker Hub
   docker pull biocontainers/samtools
   ```

   This command will download the `samtools` image and its dependencies to your local machine. We can then use the image to create containers.

2. **Run a single command non-interactively**

   We will use `samtools` to view the first few lines of a BAM file. Replace `/data_dir` with the path to the folder containing your BAM file (`align.bam`)

   ```bash
   # Run samtools view on a BAM file
   docker run --rm biocontainers/samtools -v /data_dir:/data samtools view /data/align.bam | head
   ```

   - `--rm`: Removes the container after execution
   - `-v /data_dir:/data`: Mounts the local directory `/data_dir` to the container directory `/data`

   <br>
   This command is useful for running single commands without needing an interactive shell. The `--rm` flag ensures that the container is removed after the command finishes.

3. **Run a command interactively**

   We can also run an interactive shell within the container to execute multiple commands.

   ```bash
   # Start an interactive shell in the samtools container
   docker run -it --name samtools_container biocontainers/samtools -v /data_dir:/data /bin/bash
   ```

   - `-it`: Starts an interactive terminal session
   - `/bin/bash`: Launches the bash shell in the container

   <br>
   Here, we also used the `--name` flag to give the container a name (`samtools_container`) for easy reference.
   You can now run multiple `samtools` commands within the container. When you exit the shell, the container will be removed.

### Finding Bioinformatics Tool Containers

These registries host a large number of pre-built Docker images for bioinformatics tools:

- <img src="/assets/img/posts/docker-for-bioinformatics/docker-4.svg" width="20" height="20" style="margin-right: 5px;"> [Docker Hub](https://hub.docker.com/)

- <img src="/assets/img/posts/docker-for-bioinformatics/biocontainers-logo.svg" width="20" height="20" style="margin-right: 5px;"> [Biocontainers](https://biocontainers.pro/)

- <img src="/assets/img/posts/docker-for-bioinformatics/quayio-logo.svg" width="20" height="20" style="margin-right: 5px;"> [Quay.io](https://quay.io/organization/biocontainers)

- <img src="/assets/img/posts/docker-for-bioinformatics/pegi3s-logo.svg" width="20" height="20" style="margin-right: 5px;"> [pegi3s Bioinformatics Docker Images Project](http://bdip.i3s.up.pt/)

### Composing Docker Workflows

Docker's real power shines when we use it compose complex workflows with multiple tools. By chaining together containers, we can create reproducible pipelines that can be easily shared and run on different systems.

**Bioinformatics Pipeline using `bash` scripts**

Here's an example of a pipeline that aligns reads to a reference genome using `bwa` and processes the output using `samtools`. These commands can be saved in a `bash` script for easy execution.

```bash
# Pull the bwa and samtools images
docker pull biocontainers/bwa
docker pull biocontainers/samtools

# Run the bwa aligner
docker run --rm -v /data_dir:/data biocontainers/bwa bwa mem /data/ref.fa /data/reads.fq > /data/align.sam

# Run samtools to convert the SAM file to BAM
docker run --rm -v /data_dir:/data biocontainers/samtools samtools view -bS /data/align.sam > /data/align.bam
```

**Bioinformatics Pipeline using `docker compose`**

For more complex workflows, Docker Compose provides a convenient way to define and run multi-container steps. The images and commands for the bioinformatic tools can be defined in a `docker-compose.yml` file, making it easier to manage and reproduce.

<aside>
  <p>Docker compose will typically be installed alongside Docker Desktop for Windows and macOS users. Linux users can install Docker Compose using the following commands</p>
</aside>

### Other Resources

If you'd like to read more about reproducible research practices, check out the following resources:

- ["The Turing Way"](https://book.the-turing-way.org/) : A handbook for reporoducible, ethical and collaborative data science.

Docker is a game-changer for bioinformatics, making workflows more reproducible, scalable, and shareable. Whether you’re running a single tool or a complex pipeline, Docker ensures that your research remains reliable and accessible across different environments.
