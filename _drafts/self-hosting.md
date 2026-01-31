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
        {% include figure.liquid loading="eager" path="assets/img/posts/ditch-the-cloud/ditch-the-cloud-cover.webp" class="img-fluid rounded z-depth-1" zoomable=false %}
    </div>
</div>

We live in a data-rich world and web services form the backbone of our digital lives, helping us store and manage our data as well as communicate and collaborate with others. Large data companies often provide these services for free, but we are still bound by their terms of service and privacy policies. The trade-off for using these services is often our privacy.

> In a nutshell, **Self-hosting** refers to running your own servers and software, rather than relying on third-party services.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/raspberry_pi_icon.png" class="tool-icon" alt="Raspberry Pi image">
    <a href="https://code.visualstudio.com/" target="_blank" rel="noopener">Raspberry Pi 4B</a>
</span>

I started my self-hosting journey with a Raspberry Pi 4 Model B, a small, affordable and energy-efficient computer. The Raspberry Pi immediately became my sandboxing tool allowing me to experiment with various environments without risking my systems. The Pi ecosystem had matured by this time so that I did not have to be restricted to the official OS, Raspbian, but could install full-blown Ubuntu Server and get experimenting.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/hp_elitedesk.png" class="tool-icon" alt="HP EliteDesk 800 image">
    <a href="https://code.visualstudio.com/" target="_blank" rel="noopener">HP EliteDesk 800</a>
</span>

While the Pi was perfectly capable and handled most simple tasks well, as the stack of my self-hosted services grew, I found the Pi quite limiting.

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/apple_mac_mini.png" class="tool-icon" alt="Apple Mac Mini M1 image">
    <a href="https://code.visualstudio.com/" target="_blank" rel="noopener">Apple Mac Mini M1</a>
</span>

My current self-hosting setup is on an Apple Mac Mini M1

Important pre-requisites

<span class="tool">
    <img src="/assets/img/posts/bioinformaticians-toolkit/docker.svg" class="tool-icon" alt="Docker logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Docker</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/cloudflare.svg" class="tool-icon" alt="Cloudflare logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Cloudflare</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/nginxproxymanager.svg" class="tool-icon" alt="Nginx Proxy Manager logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Ngnix Proxy Manager</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/tailscale.svg" class="tool-icon" alt="Tailscale logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Tailscale</a>
</span>

Some of the self-hosted services that've become a part of my daily workflow include:

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/immich.svg" class="tool-icon" alt="Immich logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Immich</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/karakeep.svg" class="tool-icon" alt="Karakeep logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Karakeep</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/nextcloud.svg" class="tool-icon" alt="Nextcloud logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Nextcloud</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/paperlessngx.svg" class="tool-icon" alt="Paperless-NGX logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Paperless-NGX</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/vaultwarden-light.svg" class="tool-icon" alt="Vaultwarden logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Vaultwarden</a>
</span>

<span class="tool">
    <img src="/assets/img/posts/ditch-the-cloud/vikunja.svg" class="tool-icon" alt="Vikunja logo">
    <a href="https://github.com/" target="_blank" rel="noopener">Vikunja</a>
</span>
