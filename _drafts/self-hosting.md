---
layout: post
title: "Ditch the Cloud: Self-hosting"
date: 2020-08-17 21:01:00
description: "Taking back privacy"
tags: self-hosting
categories: tech
thumbnail: assets/img/posts/ditch-the-cloud/ditch-the-cloud-thumbnail.webp
giscus_comments: true
---

<style>
.tool-icon-large {
  height: 70px;
  width: 100px;
  vertical-align: text-bottom;
  margin-right: 8px;
}
.tool-icon {
  height: 32px;
  width: 32px;
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
        {% include figure.liquid loading="eager" path="assets/img/posts/ditch-the-cloud/ditch-the-cloud-cover.webp" class="img-fluid rounded z-depth-1" zoomable=false %}
    </div>
</div>

We live in a data-rich world and web services form the backbone of our digital lives, helping us store and manage our data as well as communicate and collaborate with others. Large data companies often provide these services for free, but we are still bound by their terms of service and privacy policies. The trade-off for using these services is often our privacy.

> In a nutshell, **Self-hosting** refers to running your own servers and software, rather than relying on third-party services.

### My Self-hosting journey

Through the years, I've experimented with different hardware and software stacks for self-hosting, each with its own set of advantages and challenges.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/RPi_4.png" class="tool-icon-large" alt="Raspberry Pi image" style="width: 80px; height: auto;">
    <a href="https://code.visualstudio.com/" target="_blank" rel="noopener">Raspberry Pi 4B</a>
</span>

I started my self-hosting journey with a Raspberry Pi 4 Model B, a small, affordable and energy-efficient computer. The Raspberry Pi immediately became my sandboxing tool allowing me to experiment with various environments without risking my systems. The Pi ecosystem had matured by this time so that I did not have to be restricted to the official OS, Raspbian, but could install full-blown Ubuntu Server and get experimenting. It served as the perfect learning environment as I still faced a lot of hiccups along the way surrounding compatibility issues with the ARM-based architecture and the limitations of the 32-bit OS that the Pi could support at the time.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/HP_EliteDesk_1.png" class="tool-icon-large" alt="HP EliteDesk 800 image" style="width: 80px; height: auto;">
    <a href="https://code.visualstudio.com/" target="_blank" rel="noopener">HP EliteDesk 800</a>
</span>

While the Pi was perfectly capable and handled most simple tasks well on my local network, as the stack of my self-hosted services grew, I found the Pi quite limiting. The ARM-based architecture meant that sometimes I had to build ARM-compatible versions of software myself and even after the release of a 64-bit OS for the Pi, the Pi's hardware limitations were still a bottleneck for some of the more resource-intensive services I wanted to run. I soon upgraded my setup to a used HP EliteDesk 800 with an Intel i5 processor and 32 GB of RAM. I could install a hypervisor like Proxmox to run multiple virtual machines, testing out different OSes and container orchestration tools. The presence of an integrated GPU meant media processing was also much faster than on the Pi, especially pertinent for remote connections.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/Mac_Mini_M1.png" class="tool-icon-large" alt="Apple Mac Mini M1 image" style="width: 80px; height: auto;">
    <a href="https://code.visualstudio.com/" target="_blank" rel="noopener">Apple Mac Mini M1</a>
</span>

Moving continents, I had to rebuild my self-hosting setup from scratch and I decided to go with an Apple Mac Mini, initially starting with the M1 model and recently upgrading to the Mac Mini M4. The attraction of the Mac Mini was its compact form factor, energy efficiency, native transcoding abilities.

### Important pre-requisites

Self-hosting before the containerization era meant setting up each service manually, installing dependencies directly on the host system, editing configuration files in different places, and hoping one update did not break everything else. The learning curve was steep and the maintenance, a nightmare.

The advent of containerization tools like Docker and orchestration platforms like Kubernetes has revolutionized the self-hosting landscape.

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/docker.svg" class="tool-icon" alt="Docker logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Docker</a>
</span>

Docker forms the foundation of my current self-hosting setup. Each service runs in its own container, with Docker Compose files describing the images, storage volumes, networks and environment variables it requires.

This keeps applications relatively isolated from one another and makes the configuration easier to reproduce if I need to rebuild the server. I keep each service in its own directory, alongside its Compose file and any non-sensitive configuration required to deploy it.

I have covered the underlying concepts and commands in more detail in my post on [Docker for Bioinformatics](https://chatgpt.com/blog/2025/docker-for-bioinformatics/). Although that guide focuses on scientific workflows, the same concepts of images, containers, volumes and Compose files apply here.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/cloudflare.svg" class="tool-icon" alt="Cloudflare logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Cloudflare</a>
</span>

Cloudflare manages DNS for my domains. DNS records connect readable addresses such as `photos.example.com` with the public address through which a service can be reached.

For some services, Cloudflare also sits between visitors and my server as an additional proxy layer. This can provide useful security and traffic-management features, although services involving large uploads, long-running requests or real-time synchronization require more careful configuration.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/nginxproxymanager.svg" class="tool-icon" alt="Nginx Proxy Manager logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Ngnix Proxy Manager</a>
</span>

NGINX Proxy Manager provides a graphical interface for configuring NGINX as a reverse proxy. It receives requests for my public subdomains, obtains and renews HTTPS certificates and forwards each request to the appropriate container.

Without a reverse proxy, I would need to expose a different port for every service and remember addresses such as `example.com:2283` or `example.com:8080`. NGINX Proxy Manager instead allows each application to have a conventional address while keeping its internal port hidden behind the proxy.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/tailscale.svg" class="tool-icon" alt="Tailscale logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Tailscale</a>
</span>

Tailscale creates a private network between my trusted devices. It allows my laptop, phone, tablet and home server to communicate securely even when they are connected to different physical networks.

This means that services intended only for my own use do not need to be exposed to the public internet. As long as a device is connected to my Tailscale network, I can access the service using the server's private Tailscale address.

### Reverse-Proxy

Not everything in my stack needs to be reachable from outside my home network, but some things do — checking Paperless-NGX for a document while I'm out, or sharing a Nextcloud folder with family who aren't going to install Tailscale. I lean on two different approaches depending on who needs access.

### Tailscale, for my own devices

Since Tailscale's already part of my day-to-day, it's the default answer for anything I only need to reach from a device I control — no domain, no certificate, nothing exposed to the public internet. I've written up the actual setup steps for this pattern in my post on [self-hosting a Logseq sync server](/blog/2026/self-hosting-logseq-sync/), which uses exactly this approach for syncing between my own devices.

### DNS + Nginx Proxy Manager, for everyone else

For services I want to share with people who aren't on my tailnet — or reach from a device where installing Tailscale isn't an option — I put Nginx Proxy Manager in front of them, with Cloudflare as the DNS provider pointing subdomains at my server.

The one piece of groundwork worth calling out here: NPM needs its own Docker network, and every service stack I want to expose this way joins that same network rather than publishing ports on the host — NPM then reaches each one internally by container name. Once that's wired up, adding a new service is just: point a subdomain at it, add a proxy host in NPM, done.

I won't repeat the full click-by-click here — I've already walked through creating an NPM proxy host end-to-end, including the nginx directives you need for things like WebSockets and long-running requests, in the [Logseq sync post](/blog/2026/self-hosting-logseq-sync/#-remote-access-with-a-reverse-proxy).

Between the two: Tailscale is zero-setup and stays entirely private, but only reaches my own devices. NPM gets me a real URL anyone can hit, at the cost of a domain and slightly more moving parts. Most of what I run gets both — Tailscale for me, a proxied subdomain for anything I share.

### Self-hosted services I've come to rely on

Some of the self-hosted services that've become a part of my daily workflow include:

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/immich.svg" class="tool-icon" alt="Immich logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Immich</a>
</span>

My Google Photos/iCloud Photos replacement — automatic backup from my phone, plus ML-based search and face grouping.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/karakeep.svg" class="tool-icon" alt="Karakeep logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Karakeep</a>
</span>

A Pocket/Raindrop-style bookmark manager — I dump links, articles, and notes into it and let it auto-tag them.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/nextcloud.svg" class="tool-icon" alt="Nextcloud logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Nextcloud</a>
</span>

My Google Drive/Dropbox replacement for file sync and sharing.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/paperlessngx.svg" class="tool-icon" alt="Paperless-NGX logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Paperless-NGX</a>
</span>

OCRs and tags every scanned document and receipt I throw at it, so I can actually find them again later.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/vaultwarden-light.svg" class="tool-icon" alt="Vaultwarden logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Vaultwarden</a>
</span>

A lightweight, Bitwarden-compatible password manager server — my 1Password/LastPass replacement.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/vikunja.svg" class="tool-icon" alt="Vikunja logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Vikunja</a>
</span>

My Todoist/Trello-style task and project tracker.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/audiobookshelf.svg" class="tool-icon" alt="Audiobookshelf logo">
    <a href="https://www.audiobookshelf.org/" target="_blank" rel="noopener">Audiobookshelf</a>
</span>

Self-hosted audiobook and podcast server — my alternative to Audible/Spotify for that library.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/syncthing-dark.svg" class="tool-icon" alt="Syncthing logo">
    <a href="https://syncthing.net/" target="_blank" rel="noopener">Syncthing</a>
</span>

Peer-to-peer file sync between my own devices, no cloud storage or middleman server involved.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/pairdrop.webp" class="tool-icon" alt="PairDrop logo">
    <a href="https://pairdrop.net/" target="_blank" rel="noopener">PairDrop</a>
</span>

A self-hosted AirDrop-style tool for quickly sending files between devices on the same network.
