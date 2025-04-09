---
layout: distill
title: "Docker for Bioinformatics"
date: 2025-02-27 21:01:00
description: "Portable, Scalable and Reproducible Bioinformatics workflows"
tags: docker bioinformatics
categories: bioinformatics
thumbnail: assets/img/posts/docker-for-bioinformatics/docker-for-bioinformatics-thumbnail.png
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
        {% include figure.liquid loading="eager" path="assets/img/posts/docker-for-bioinformatics/docker-for-bioinformatics-cover.png" class="img-fluid rounded z-depth-1" zoomable=false %}
    </div>
</div>

Bioinformatics analysis often involves complex pipelines with rapidly evolving software tools, each with their own set of dependencies. System compatibility, version mismatches and dependency conflict issues can often be a nightmare, making running and sharing bioinformatic pipelines a challenging task. These challenges not only waste valuable research time but also contribute to irreproducible workflows, where results depend as much on the computing environment as on the analysis itself. Docker offers a powerful solution by packaging software and its dependencies into portable, reproducible containers, ensuring that your bioinformatics pipelines run consistently, whether on your local machine, an HPC cluster, or the cloud.

## What is Docker?

Imagine you're baking a cake, but every time you try, your kitchen is missing key ingredients or uses a different oven that bakes at the wrong temperature. Docker is like a self-contained baking kit that comes with all the right ingredients, tools, and even its own portable oven, ensuring your cake turns out exactly the same no matter where you bake it. In bioinformatics, Docker does the same for software by packaging tools, dependencies, and environments so that analyses run reliably across different computing platforms.

## How can Docker help?

The **FAIR** (Findable, Accessible, Interoperable, and Reusable) principles <d-cite key="wilkinson_fair_2016"></d-cite> provide guidelines for maximizing the value of research data. Docker aligns bioinformatics workflows with these principles by ensuring software and environments are portable and reproducible:

- **Findability**: Docker images can be easily found in registries like Docker Hub, with clear versioning and documentation for discovery.

- **Accessibility**: Anyone with internet access and Docker installed can retrieve and use containerized tools, promoting open science.

- **Interoperability**: Docker containers provide a consistent runtime environment across different systems, preventing dependency conflicts.

- **Reusability**: Pre-built images allow researchers to reuse workflows without worrying about installation issues, fostering collaboration.

## Getting Started

To begin, start by installing Docker on your system. Docker is available for all major operating systems and the installers can be downloaded from the [official website](https://www.docker.com/get-started/). For Windows and macOS users, the recommended approach is to install the Docker Desktop application, while Linux users can install Docker natively for a more lightweight setup.<d-footnote>Docker Desktop creates a Linux virtual machine (VM) on Windows and macOS to run containers, whereas on a Linux machine, Docker runs natively without the need for a VM.</d-footnote>

<hr style="grid-column: text; width: 100%; border: none; border-bottom: 1px solid rgba(0, 0, 0, 0.1); margin-top: 1rem; margin-bottom: 1rem;">

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

**_Note:_** If you do not have `sudo` privileges, you can install Docker using the [official script](https://docs.docker.com/engine/install/ubuntu/#install-using-the-convenience-script).

```bash
# Download the Docker installation script
curl -fsSL https://get.docker.com -o get-docker.sh

# Append --dry-run to see the commands that will be executed
sudo sh get-docker.sh
```

{% endtab %}

{% endtabs %}

<hr style="grid-column: text; width: 100%; border: none; border-bottom: 1px solid rgba(0, 0, 0, 0.1); margin-top: 1rem; margin-bottom: 1rem;">

<br>
To test whether Docker is installed correctly, run the following command in your terminal:

```bash
# Check Docker version
docker --version

# Test with a simple hello-world container
docker run hello-world
```

<aside>
  <p>💡 <strong>Tip:</strong> Running <code>docker --help</code> will display a list of available commands and options.</p>
</aside>

If installed correctly, these commands will print the version of Docker installed on your system and fetch the image and run the `hello-world` container, which prints a message confirming that Docker is working.

## Understanding Key Docker Concepts

- **Images** vs **Containers**:

  A _Docker image_ is a static, read-only blueprint that includes the application and all its dependencies. A _container_ is a live, running instance created from that image. Returning to our baking analogy: the _image_ is your recipe and ingredients kit, while the _container_ is the oven actively baking the cake. You can spin up multiple containers from the same image—just like baking several cakes from one recipe.

- The importance of **Volumes**:

  By default, docker containers are **ephemeral** i.e. once they stop any files written inside them are lost. To effectively manage input/output tasks between the host machine and the container and to persist data, we use _volumes_. We mount a local directory on the host machine to a directory inside the container using the `-v` (or the more flexible `--mount`) flag.

## Running Bioinformatics Tools with Docker

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

   This command is useful for running single commands without needing an interactive shell. The `--rm` flag ensures that the container is removed after the command finishes.

3. **Run a command interactively**

   We can also run an interactive shell within the container to execute multiple commands.

   ```bash
   # Start an interactive shell in the samtools container
   docker run -it --name samtools_container biocontainers/samtools -v /data_dir:/data /bin/bash
   ```

   - `-it`: Starts an interactive terminal session
   - `/bin/bash`: Launches the bash shell in the container

   Here, we also used the `--name` flag to give the container a name (`samtools_container`) for easy reference.

   We can now run multiple `samtools` commands within the container:

   ```bash
   # Check samtools version
   samtools --version

   # Index a reference genome
   samtools faidx /data/ref.fa
   ```

   We can add `--rm` to the interactive `docker run` command to remove the container after exiting the shell.

## Composing Docker Workflows

Docker's real power shines when we use it compose complex workflows with multiple tools. By chaining together containers, we can create reproducible pipelines that can be easily shared and run on different systems.

Let's consider the first step of most bioinformatics workflows: quality control of sequencing reads. This step is often performed using tools like [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)<d-footnote><a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank">FastQC</a> is a quality control tool that analyzes raw sequence data from high throughout sequencing runs</d-footnote> with results conveniently summarized using tools like [`multiqc`](https://seqera.io/multiqc/)<d-footnote><a href="https://seqera.io/multiqc/" target="_blank">MultiQC</a> aggregates results and quality metrics from multiple bioinformatics analysis reports (often including those from <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank">FastQC</a>) into a single, interactive summary report, facilitating comparison across numerous samples or steps.</d-footnote>. These tools can be run in a Docker container, allowing us to easily check the quality of our sequencing data.

To try out the pipeline with real data, we can use test FASTQC files from the [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repository that are ideal for quick pipeline tests. We can run the following commands one-by-one on the command line or save them in a `bash` script to download the test data to the `~/docker-bioinf/data/raw_data` directory. You can replace this with any other directory of your choice.

```bash
# Create directory for test data
mkdir -p ~/docker-bioinf/data/raw_data
cd ~/docker-bioinf/data/raw_data

# Download test FASTQ files
wget https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/test_1.fastq.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/test_2.fastq.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/test2_1.fastq.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/test2_2.fastq.gz
```

<aside>
  <p> 💡 <strong>Note:</strong> If <code>wget</code> is not installed on your system, you can replace it with <code>curl -O</code> in the download commands. For example, <code>wget URL</code> becomes <code>curl -O URL</code>.</p>
</aside>

Check the contents of the `~/docker-bioinf/data/raw_data` directory to confirm that the files have been downloaded successfully.

```bash
# Check the contents of the directory
ls ~/docker-bioinf/data/raw_data
```

### Bioinformatics Pipeline using `bash` scripts

We can chain together multiple docker commands to construct a lightweight, portable pipeline for quality control of sequencing reads. The following `bash` script demonstrates how to run `fastqc` on all FASTQ files in the `~/docker-bioinf/data/raw_data` directory and generate a summary report using `multiqc`. The script will create a new directory called `qc_reports` to store the output reports.

```bash
#!/bin/bash
# Create the output directory
mkdir -p ~/docker-bioinf/data/qc_reports

# Pull the FastQC and MultiQC containers
docker pull quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
docker pull quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0

# Run FastQC on all FASTQ files in the raw_data directory
docker run --rm -v ~/docker-bioinf/data:/data quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0 \
  bash -c 'fastqc /data/raw_data/*.fastq.gz -o /data/qc_reports'

# Run MultiQC to aggregate FastQC reports
docker run --rm -v ~/docker-bioinf/data:/data quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0 \
  multiqc /data/qc_reports -o /data/qc_reports
```

<aside> <p> 💡 <strong>Note:</strong> If your Docker installation does not have root privileges, it may not be able to create new directories inside mounted volumes. To avoid errors, <strong>manually create the <code>qc_reports</code> directory</strong> on the host system before running the pipeline: </p> <pre style="font-size: 0.85em; line-height: 1.4;"><code>mkdir -p ~/docker-bioinf/data/qc_reports</code></pre> </aside>

> Note the use of `bash -c` when running the `fastqc` command which ensures that the command is executed in a shell environment thereby enabling the expansion of wildcards (e.g. _\*fastq.gz_)

The `fastqc` and `multiqc` reports will be saved in the `~/docker-bioinf/data/qc_reports` directory, and you can view them using any web browser. The `fastqc` reports will be in HTML format, while the `multiqc` report will be an interactive HTML file (`multiqc_report.html`) that aggregates the results from all the `fastqc` reports.

### Bioinformatics Pipeline using `docker compose`

For more complex workflows, Docker Compose provides a convenient way to define and run multi-container steps with built-in dependency management. The images and commands for the bioinformatic tools can be defined in a `docker-compose.yml` file, making it easier to manage and reproduce.

<aside>
  <p> Docker compose will typically be installed alongside Docker Desktop for Windows and macOS users. Linux users can install Docker Compose plugin by following this <a href="https://docs.docker.com/compose/install/linux/" target="_blank">link</a>.</p>
</aside>

Let's revisit the earlier bioinformatics task: running FastQC on raw FASTQ files and summarizing the results using MultiQC. Instead of invoking each tool manually with separate `docker run` commands, we can streamline the workflow using a `docker-compose.yml` file:

```yaml
services:
  fastqc:
    image: quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
    volumes:
      - ~/docker-bioinf/data:/data
    entrypoint: bash -c
    command: >
      "mkdir -p /data/qc_reports &&
      fastqc /data/raw_data/*.fastq.gz -o /data/qc_reports"

  multiqc:
    image: quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0
    volumes:
      - ~/docker-bioinf/data:/data
    command: multiqc /data/qc_reports -o /data/qc_reports
    depends_on:
      fastqc:
        condition: service_completed_successfully
```

This `docker-compose.yml` file defines the following key sections:

- `services`: Defines the containers in our pipeline. Each service is a separate Docker container.
- `fastqc`: and `multiqc`: These are our service definitions.
  - `image`: Specifies the Docker image to use (e.g., `quay.io/biocontainers/fastqc`).
  - `volumes`: Mounts a local directory (`~/docker-bioinf/data`) to the container's `/data` directory
  - `entrypoint`: Sets the entrypoint for the container to `bash -c`, allowing us to run multiple commands in a single container.
  - `command`: Sets the command to run within the container.
- `depends_on`: This directive is essential for managing the workflow's order of operations.
  - `fastqc`: Indicates that the `multiqc` service relies on the output of the `fastqc` service.
    - `condition: service_completed_successfully`: By setting this condition, we instruct Docker Compose to only initiate the `multiqc` container once the `fastqc` container has finished its job and exited with a success status. This prevents `multiqc` from running prematurely and encountering an empty or incomplete `qc_reports` directory.

To run the pipeline, execute the following command in the same directory as the `docker-compose.yml` file:

```bash
docker-compose up
```

> 💡 Tip: If you're re-running the pipeline and want a clean start, use <code>docker compose down</code> to remove the containers, or add the <code>--force-recreate</code> flag when running up (e.g. <code> docker compose up --force-recreate</code>).

## Beyond Pre-Built Images: Introducing Dockerfiles

While pre-built Docker images are incredibly helpful, they may not always meet your specific needs. Suppose you want to install a specific version of a tool or dependency, package custom scripts along with the tools or simply that the pre-built image is not available. In such cases, you can create your own Docker images using a `Dockerfile`, which is a text file that contains instructions for building a Docker image.
It specifies the base image, the software to install, and any configuration needed to set up the environment.

We can create our own custom docker image for the same `fastqc` and `multiqc` pipeline using a `Dockerfile`. This provides us complete control to customize the environment, install additional dependencies, and package our scripts along with the tools.

```dockerfile
# Use lightweight linux base
FROM debian:bullseye-slim

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-11-jdk \
    python3 \
    python3-pip \
    bash \
    wget \
    unzip \
    perl \
    libperl-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    mv FastQC /opt/fastqc && \
    chmod +x /opt/fastqc/fastqc && \
    ln -s /opt/fastqc/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.12.1.zip

# Install MultiQC
RUN pip3 install --no-cache-dir multiqc

# Set working directory
WORKDIR /data

# Make shell commands easier to write
ENTRYPOINT ["bash", "-c"]
```

The `Dockerfile` contains the instructions to build an environment with the tools and necessary dependencies for our analysis, the key sections are:

- `FROM <image>`: Specifies the base image to use, this is the image we'll extend by installing our specific tools. In this case, we are using a lightweight Debian image.
<aside><p>Explore more base image options on <a href="https://docs.docker.com/docker-hub/image-library/trusted-content/#docker-official-images" target="_blank">Docker Hub</a>.</p></aside>
- `RUN <command>`: Executes a command in the container during the build process. We use this to install dependencies and tools.
- `WORKDIR <directory>`: Sets the working directory inside the container. This is where the commands will be executed.
- `ENTRYPOINT <command>`: Sets the default command to run when the container starts. In this case, we set it to `bash -c` to allow us to run multiple commands.

A complete list of Dockerfile instructions can be found in the [Dockerfile reference](https://docs.docker.com/engine/reference/builder/).

Next, we build the Docker image using the `docker build` command. The `-t` flag allows us to tag the image with a name (e.g., `my_fastqc_multiqc`).

```bash
# Build the Docker image
docker build -t my_fastqc_multiqc .
```

Once the image is built, we can run it using the same commands as before. The only difference is that we will use our custom image name instead of the pre-built one.

```bash
mkdir -p ~/docker-bioinf/data/qc_reports

# Run FastQC on all FASTQ files in the raw_data directory
docker run --rm -v ~/docker-bioinf/data:/data my_fastqc_multiqc \
  bash -c 'mkdir -p /data/qc_reports && fastqc /data/raw_data/*.fastq.gz -o /data/qc_reports'

# Run MultiQC to aggregate FastQC reports
docker run --rm -v ~/docker-bioinf/data:/data my_fastqc_multiqc \
  multiqc /data/qc_reports -o /data/qc_reports
```

## Best Practices for Docker in Bioinformatics

1. **Always Use Specific Image Versions**:

   Use a versioned image tag (like `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0`) instead of the `latest` tag. This ensures your workflow always uses the exact same version of the tool every time it's run and avoids unexpected changes in behavior, guaranteeing reproducibility.

2. **Leverage Biocontainers**:

   Before searching elsewhere or attempting to build an image, check repositories like [Quay.io/biocontainers](https://quay.io/organization/biocontainers). Utilizing these standardized, pre-built images for common bioinformatics tools saves significant effort and aligns your workflow with community standards.

3. **Handle Data Appropriately**:

   - Separate data from the container: Use volumes to mount data directories from the host system into the container, rather than copying it to the container's filesystem. This keeps the container immutable and promotes reusability with different datasets.
   - Mount reference data as read-only volumes to prevent accidental modifications.
   - Use named volumes for persisting data between container runs and to share data between multiple containers.

## Finding Bioinformatics Tool Containers

These registries host a large number of pre-built Docker images for bioinformatics tools:

- <img src="/assets/img/posts/docker-for-bioinformatics/docker-4.svg" width="20" height="20" style="margin-right: 5px;"> [Docker Hub](https://hub.docker.com/) : The default Docker registry, home to many official and community-maintained bioinformatics tool images.

- <img src="/assets/img/posts/docker-for-bioinformatics/biocontainers-logo.svg" width="20" height="20" style="margin-right: 5px;"> [Biocontainers](https://biocontainers.pro/) : A community-driven project offering thousands of standardized containers for bioinformatics tools, ideal for reproducible pipelines.

- <img src="/assets/img/posts/docker-for-bioinformatics/quayio-logo.svg" width="20" height="20" style="margin-right: 5px;"> [Quay.io](https://quay.io/organization/biocontainers) : Another major container registry, widely used by the BioContainers project to host their images.

- <img src="/assets/img/posts/docker-for-bioinformatics/pegi3s-logo.svg" width="20" height="20" style="margin-right: 5px;"> [pegi3s Bioinformatics Docker Images Project](http://bdip.i3s.up.pt/) : A curated collection of Docker images focused on reproducibility and ease of use.

- <img src="/assets/img/posts/docker-for-bioinformatics/rocker-logo.svg" width="20" height="20" style="margin-right: 5px;"> [Rocker Project](https://rocker-project.org/) : Docker images tailored for R and RStudio users, commonly used in data science and bioinformatics research.

## Resources and Further Reading

**Docker Essentials and Learning:**

- [Getting Started with Docker](https://docs.docker.com/get-started/) : A beginner-friendly guide to Docker, covering installation and basic commands.
- [Docker CLI CheatSheet](https://docs.docker.com/get-started/docker_cheatsheet.pdf) : A quick reference guide for common Docker commands.
- [Docker Documentation](https://docs.docker.com/) : The official Docker documentation is a comprehensive resource for all things Docker.
- [Docker Compose Documentation](https://docs.docker.com/compose/) : A guide to using Docker Compose for multi-container applications.
- [BioContainers Best Practices](https://biocontainers-edu.readthedocs.io/en/latest/best_practices.html): A guide to best practices for using BioContainers in bioinformatics workflows.

**Guided Lessons:**

- [Software Carpentries: Introduction to Docker](https://carpentries-incubator.github.io/docker-introduction/) : A hands-on lesson designed for researchers new to containers.
- [Linux containers in scientific environments (CBG PhD Course)](https://biocorecrg.github.io/PhD_course_containers_2021/): Short hands-on practicum on how to start working using Linux containers.

**Reproducibility in Research Practices:**

- [The FAIR Guiding Principles](https://www.nature.com/articles/sdata201618)<d-cite key="wilkinson_fair_2016"></d-cite>: The original paper outlining the FAIR principles
- [Ten Simple Rules for Reproducible Computational Research](https://doi.org/10.1371/journal.pcbi.1000424)<d-cite key="sandve_ten_2013"></d-cite> : A paper outlining ten simple rules for reproducible research in computational biology.
- [Recommendations for the packaging and containerizing of bioinformatics software](https://f1000research.com/articles/7-742/v2)<d-cite key="gruening_recommendations_2019"></d-cite> : A paper discussing best practices for packaging and containerizing bioinformatics software to ensure reproducibility.
- ["The Turing Way"](https://book.the-turing-way.org/) : A handbook for reporoducible, ethical and collaborative data science.
- [The Open Science Manual](https://arca-dpss.github.io/manual-open-science/): A guide to open science practices, including reproducibility and data sharing.

## Conclusion

Bioinformatic workflows often suffer from "dependency hell", wherein conflicts between software libraries, incompatible versions and platform-specific quirks can make setting up and running analyses a frustrating experience. Containerization technologies like Docker provide a powerful solution by encapsulating the software and it's dependencies along with any necessary configurations into a single, portable package. This ensures that the analysis runs consistently across different environments, making the workflows more reproducible, scalable, and shareable. Whether you're running a single tool or a complex pipeline, Docker ensures that your research remains reliable and accessible across different environments.

Happy Dockering!
