---
layout: post
title: "Self-Hosting a Logseq Sync Server"
date: 2026-07-18 10:00:00
description: "Local access with Tailscale, remote access with a reverse proxy"
tags: logseq self-hosting
categories: tech
thumbnail: assets/img/posts/self-hosting-logseq-sync/self-hosting-logseq-sync-thumbnail.webp
giscus_comments: true
disable_animation: true
toc:
  sidebar: left
_styles: >
  .container {
    max-width: 1200px;
  }
---

<style>
.option-icon--tailscale {
  display: inline-block;
  width: 24px;
  height: 24px;
  margin-right: 6px;
  vertical-align: middle;
  background-color: #242424;
  -webkit-mask: url("/assets/img/posts/ditch-the-cloud/tailscale.svg") center / contain no-repeat;
  mask: url("/assets/img/posts/ditch-the-cloud/tailscale.svg") center / contain no-repeat;
}

html[data-theme="dark"] .option-icon--tailscale {
  background-color: #828282;
}

.tip-icon {
  color: #efcc00;
  margin-right: 6px;
}

.info-icon {
  color: #0076df;
  margin-right: 6px;
}

.success-icon {
  color: #28a745;
}
</style>

<div class="row justify-content-center mt-3">
    <div class="col-12 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/self-hosting-logseq-sync-cover.webp" class="img-fluid rounded z-depth-1" zoomable=false %}
    </div>
</div>

[Logseq](https://logseq.com/) has undergone some transformative changes over the past few years, with the [switch to the Database (DB) version](https://logseq.io/p/e3YDyX5AYr) that brings better performance, structured data entry, Real Time Collaboration (RTC) and a host of other features. The move away from storing notes in local markdown files also means that users can no longer rely on solutions like cloud folders (Google Drive, iCloud etc.) or Syncthing for syncing their notes across devices but now have to rely on Logseq's official Sync services to enable cross-device access.

The official Logseq Sync service is a paid service that can currently be enabled by contributing to [Logseq's Open Collective](https://opencollective.com/logseq) at the "Backers" and "Sponsors" tiers. However, in the true spirit of Open Source, the developers have gracefully provided the infrastructure to self-host the Logseq Sync and Publish servers on our own hardware and the community has also built Docker containers ([yshalsager's `logseq-selfhost`](https://github.com/yshalsager/logseq-selfhost)) to make the self-hosted setup simpler.

This post walks through two ways of exposing that server:

- **Private access over [Tailscale](https://tailscale.com/)**: the quickest way to sync between your own devices without exposing the server to the public internet.
- **Public access through a [reverse proxy](https://www.cloudflare.com/learning/cdn/glossary/reverse-proxy/)**: useful when you need to sync from a device that cannot run Tailscale, such as a managed work computer.

Both approaches sit on top of the same [`logseq-selfhost`](https://github.com/yshalsager/logseq-selfhost) containers: choose the [Tailscale](#-private-access-over-tailscale) route if you're syncing between your own devices, or the [reverse proxy](#-public-access-through-a-reverse-proxy) route if you need access from somewhere Tailscale can't reach, like a managed work computer.

### Common Prerequisites

Whichever access route you choose, you will need:

- **Logseq account**: Self-hosting replaces Logseq's hosted Sync and Publish endpoints, but the current clients still require you to sign in with a Logseq account for authentication. You do not need a subscription while using the self-hosted route, only the account.
- **Device capable of running Docker**: A home server, NAS, desktop computer, or VPS.
- **Docker Engine with Compose**: On Linux, follow [Docker's official install guide](https://docs.docker.com/engine/install/); On MacOS and Windows, Docker Desktop includes the Compose plugin
- **Git**: For cloning the [`logseq-selfhost`](https://github.com/yshalsager/logseq-selfhost) repository
- **Logseq DB clients**: Client application for every device you want to sync

### <img src="/assets/img/posts/bioinformaticians-toolkit/logseq.svg" width="24" height="24" style="margin-right: 6px; vertical-align: middle;"> Creating a Logseq Account

A Logseq account can be created using any of Logseq's client applications or by navigating to [test.logseq.com](https://test.logseq.com/) <i class="fa-solid fa-arrow-right"></i> "<i class="fa-solid fa-ellipsis"></i>" <i class="fa-solid fa-arrow-right"></i> **Login** <i class="fa-solid fa-arrow-right"></i> **Sign Up**.

The account created here will be used to authenticate all other clients for the sync operations.

### <span class="option-icon--tailscale" role="img" aria-label="Tailscale logo"></span> Private access over Tailscale

[Tailscale](https://tailscale.com/) creates a private mesh network (a "tailnet") between your devices, so the sync server never needs a public IP, a domain, or a certificate: you just point your Logseq clients at the server's Tailscale IP address.

> <i class="fa-solid fa-circle-info info-icon"></i> Ensure [Tailscale](https://tailscale.com/download) is installed and signed in on the server **and** on every device you want to sync from.

#### Setup

1. **Install Tailscale on the server and on each client device**, then note the Tailscale IP address assigned to the server (visible in the Tailscale admin console, or by running `tailscale ip -4` on the server itself).

<div class="row justify-content-center mt-3">
    <div class="col-sm-8 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/tailscale-admin-console.webp" class="img-fluid rounded z-depth-1" zoomable=true caption="The server's Tailscale IP, as shown in the admin console." %}
    </div>
</div>

1. **Clone the repository and prepare the environment file:**

   ```bash
   git clone https://github.com/yshalsager/logseq-selfhost.git
   cd logseq-selfhost/images/sync
   cp .env.example .env
   ```

2. **Point the sync server at your Tailscale IP.** Open `.env` and change:

   ```
   DB_SYNC_BASE_URL=https://sync.example.com
   ```

   to:

   ```
   DB_SYNC_BASE_URL=http://YOUR_TAILSCALE_IP:8787
   ```

3. **Start the server:**

   ```bash
   docker compose pull
   docker compose up -d
   ```

4. **Verify it's running:**

   ```bash
   curl http://localhost:8787/health
   ```

   should return:

   ```bash
   {"ok":true}
   ```

Once the health endpoint is reachable from another device on your tailnet, continue to [Connecting your Logseq clients](#connecting-your-logseq-clients) and use:

```text
http://YOUR_TAILSCALE_IP:8787
```

> <i class="fa-solid fa-lightbulb tip-icon"></i> You can optionally enable MagicDNS in Tailscale and use a memorable address such as <code>http://my-server:8787</code> instead of the server's Tailscale IP.

One thing worth flagging: the sync server listens on port 8787 by default, and Tailscale does not close that port to the rest of your network, it only gives your devices a private way to reach it. If your server also has a LAN or public IP, take a moment to confirm nothing outside your tailnet can reach port 8787, either with a firewall rule or by checking your router is not forwarding it.

> <i class="fa-solid fa-lightbulb tip-icon"></i> Most firewalls can scope that rule to Tailscale's own IP range (100.64.0.0/10) rather than opening the port to everyone.

### <img src="/assets/img/posts/ditch-the-cloud/nginxproxymanager.svg" width="24" height="24" style="margin-right: 6px; vertical-align: middle;"> Public access through a reverse proxy

If you want to sync from outside your tailnet, or cannot install Tailscale on a device, you can place the same Sync server behind **NGINX Proxy Manager (NPM)** and access it through a normal HTTPS domain.

> <i class="fa-solid fa-circle-info info-icon"></i> I chose NGINX Proxy Manager (NPM) as I already have a server stack running with NPM to connect to. NPM can be replaced with your reverse-proxy manager of choice.

#### Prerequisites

- **NGINX Proxy Manager already installed and running**, with a Docker network that the Logseq container can join.
- **A domain with a DNS provider**, with the ability to create and edit DNS records (this guide uses [Cloudflare](https://dash.cloudflare.com/login))
- **Ports 80 and 443 forwarded to NGINX Proxy Manager** from your router or firewall.

#### Setup

1. **Create a deployment directory** containing `docker-compose.yml` and `.env`, with persistent data stored in a `./data` subdirectory:

   ```bash
   mkdir -p ~/home-server/logseq-sync/data
   cd ~/home-server/logseq-sync
   ```

2. **Create the environment file (`.env`):**

   ```
   DB_SYNC_BASE_URL=https://logseq-sync.example.com
   DB_SYNC_ADMIN_TOKEN=replace-with-random-token
   ```

   Generate the token with:

   ```bash
   openssl rand -hex 32
   ```

3. **Create `docker-compose.yml`:**

   ```yaml
   name: logseq-selfhost-sync

   services:
     logseq-sync:
       image: ghcr.io/yshalsager/logseq-selfhost-sync:latest
       container_name: logseq-selfhost-sync
       restart: unless-stopped
       pull_policy: always

       environment:
         DB_SYNC_PORT: "8787"
         DB_SYNC_BASE_URL: "${DB_SYNC_BASE_URL}"
         DB_SYNC_DATA_DIR: "/app/data"
         DB_SYNC_STORAGE_DRIVER: "sqlite"
         DB_SYNC_ASSETS_DRIVER: "filesystem"
         DB_SYNC_LOG_LEVEL: "info"

         COGNITO_ISSUER: "https://cognito-idp.us-east-1.amazonaws.com/us-east-1_dtagLnju8"
         COGNITO_CLIENT_ID: "69cs1lgme7p8kbgld8n5kseii6"
         COGNITO_JWKS_URL: "https://cognito-idp.us-east-1.amazonaws.com/us-east-1_dtagLnju8/.well-known/jwks.json"

         DB_SYNC_ADMIN_TOKEN: "${DB_SYNC_ADMIN_TOKEN}"

       volumes:
         - ./data:/app/data

       read_only: true

       tmpfs:
         - /tmp

       cap_drop:
         - ALL

       security_opt:
         - no-new-privileges:true

       networks:
         - default
         - nginx-proxy-manager_default

   networks:
     nginx-proxy-manager_default:
       external: true
   ```

   The network name must match the Docker network used by your NPM deployment. You can find it by running:

   ```bash
   docker network ls
   ```

4. **Set up DNS.** Create a record pointing at your infrastructure:

   ```text
   logseq-sync.example.com
   ```

5. **Start the container:**

   ```bash
   docker compose pull
   docker compose up -d
   ```

6. **Create a proxy host in NGINX Proxy Manager:**

   | Setting          | Value                     |
   | ---------------- | ------------------------- |
   | Domain name      | `logseq-sync.example.com` |
   | Scheme           | `http`                    |
   | Forward hostname | `logseq-selfhost-sync`    |
   | Forward port     | `8787`                    |

  <br>

Enable **Websockets Support** and **Block Common Exploits**

   <div class="row justify-content-center mt-3">
       <div class="col-sm-8 mt-3 mt-md-0">
           {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/npm-proxy-host.webp" class="img-fluid rounded z-depth-1 mx-auto d-block" width="70%" zoomable=true caption="The Details tab for the Sync proxy host, with Block Common Exploits and Websockets Support enabled." %}
       </div>
   </div>

Under the "SSL" tab, enable **Force SSL** and **HTTP/2 Support** and request a new certificate:

   <div class="row justify-content-center mt-3">
       <div class="col-sm-8 mt-3 mt-md-0">
           {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/npm-proxy-host-SSL.webp" class="img-fluid rounded z-depth-1 mx-auto d-block" width="70%" zoomable=true caption="The SSL tab for the Sync proxy host, with Force SSL and HTTP/2 Support enabled." %}
       </div>
   </div>

Under **Advanced** (the settings <i class="fa-solid fa-gear"></i> icon), add:

```bash
client_max_body_size 1024m;
proxy_read_timeout 3600s;
proxy_send_timeout 3600s;
```

   <div class="row justify-content-center mt-3">
       <div class="col-sm-8 mt-3 mt-md-0">
           {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/npm-proxy-host-advanced.webp" class="img-fluid rounded z-depth-1 mx-auto d-block" width="70%" zoomable=true caption="The Advanced tab for the Sync proxy host, with the custom NGINX configuration applied." %}
       </div>
   </div>

> <i class="fa-solid fa-circle-info info-icon"></i> The longer <code>proxy_read_timeout</code>/<code>proxy_send_timeout</code> values matter here: Logseq Sync operations on large graphs can take a while, and NPM's defaults will otherwise cut the connection early.

#### Verification

```bash
docker compose ps

curl -i https://logseq-sync.example.com/health
```

should return:

```bash
{"ok":true}
```

Once the public health endpoint is working, continue to [Connecting your Logseq clients](#connecting-your-logseq-clients) and use:

```text
https://logseq-sync.example.com
```

#### A note on security

If you adapt this compose file for your own setup, keep the hardening defaults intact: read-only filesystems, dropped capabilities, and no ports published directly on the host. Together they mean the container has little room to do anything beyond running the sync service, even if something in it were compromised.

### Connecting your Logseq clients

The client setup is the same for both access routes. The only difference is the Sync server URL:

| Access route  | Sync server URL                   |
| :------------ | :-------------------------------- |
| Tailscale     | `http://YOUR_TAILSCALE_IP:8787`   |
| Reverse proxy | `https://logseq-sync.example.com` |

Use the appropriate URL on every device you want to synchronize.

#### Desktop

The latest Desktop client builds can be downloaded from [Logseq's GitHub](https://github.com/logseq/logseq/releases/tag/2.0.1): choose the installer appropriate for your OS.

In Logseq, login with your Logseq account (see [Creating a Logseq Account](#creating-a-logseq-account)), then go to **Settings → Advanced → Sync Server URL** and enter your Sync server URL.

<div class="row justify-content-center mt-3">
    <div class="col-sm-8 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/logseq-custom-sync-server-url.webp" class="img-fluid rounded z-depth-1" zoomable=true caption="Setting the custom Sync Server URL in Logseq's desktop settings." %}
    </div>
</div>

If successful, you should see a cloud icon with a green dot signalling the sync status in the top right corner, next to <i class="fa-solid fa-ellipsis"></i>. Clicking on it will display the current sync status.

#### iPhone and iPad

The iOS/iPadOS clients for Logseq DB can be downloaded by joining the [TestFlight](https://testflight.apple.com/join/eBcJ9Hpc) program for the Logseq Beta app.

In Logseq, login with your Logseq account (see [Creating a Logseq Account](#creating-a-logseq-account)), then add your server url to the **Sync Server URL** setting.

<div class="row justify-content-center mt-3">
    <div class="col-sm-8 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/logseq-custom-sync-server-url_iOS.webp" class="img-fluid rounded z-depth-1 mx-auto d-block" width="50%" zoomable=true caption="Setting the custom Sync Server URL in Logseq for iOS and iPadOS." %}
    </div>
</div>

If successful, you should see a cloud icon with a green dot (<i class="fa-solid fa-circle success-icon"></i>) signalling the sync status in the top right corner, next to <i class="fa-solid fa-ellipsis"></i>. Clicking on it will display the current sync status.

#### Android

The latest Desktop client builds can be downloaded from [Logseq's GitHub](https://github.com/logseq/logseq/releases/tag/2.0.1): choose the `.apk` file from the Assets list.

In Logseq, login with your Logseq account (see [Creating a Logseq Account](#creating-a-logseq-account)), then add your server url to the **Sync Server URL** setting.

 <div class="row justify-content-center mt-3">
    <div class="col-sm-8 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/logseq-custom-sync-server-url_android.webp" class="img-fluid rounded z-depth-1 mx-auto d-block" width="40%" zoomable=true caption="Setting the custom Sync Server URL in Logseq for Android." %}
       </div>
   </div>

If successful, you should see a cloud icon with a green dot (<i class="fa-solid fa-circle success-icon"></i>) signalling the sync status in the top right corner, next to <i class="fa-solid fa-ellipsis"></i>. Clicking on it will display the current sync status.

> <i class="fa-solid fa-circle-info info-icon"></i> If using <code>Tailscale</code>, install <code>Tailscale</code> and sign in to the same tailnet.

> <i class="fa-solid fa-circle-info info-icon"></i> On every device, make sure **Use Logseq Sync Beta** is enabled before creating a new graph or synchronizing an existing one.

### Optional: Web and Publish servers

The Sync server is all you need to synchronize a DB graph. The same `logseq-selfhost` project also provides two optional services:

- `logseq-selfhost-web` — the web application interface.
- `logseq-selfhost-publish` — the Publish server for sharing a graph.

If you need either feature, replace the Sync-only deployment with the full **Sync + Web + Publish** stack in my [`self-hosted-docker-setups` repository](https://github.com/abhilesh/self-hosted-docker-setups/tree/main/logseq-selfhost). The expanded stack uses the same `./data` directory for Sync and attaches all three containers to NPM's network.

Create two additional DNS records:

- `logseq-web.example.com`
- `logseq-publish.example.com`

Then create proxy hosts pointing to:

| Service | Forward hostname          | Internal port |
| ------- | ------------------------- | ------------: |
| Web     | `logseq-selfhost-web`     |        `8080` |
| Publish | `logseq-selfhost-publish` |        `8787` |

For both hosts, enable **Force SSL** and **HTTP/2 Support**. The Web proxy can use the same timeout configuration as Sync, while the Publish proxy should also allow long-running requests.

In Logseq, set:

- **Publish Server URL:** `https://logseq-publish.example.com`

> <i class="fa-solid fa-circle-info info-icon"></i> While the Publish server can also be exposed over Tailscale, the published pages will only be accessible on devices connected to your tailnet.

You can see the Publish server in action in the [previous version of this guide](https://logseq-publish.biophiles.lol/p/fmT_7lx_Kz), which is hosted using the same setup. That version also documents the complete combined Sync, Web, and Publish stack behind NGINX Proxy Manager.

For most personal setups, Tailscale is the simplest place to start. If you later need the reverse proxy route or want to publish Logseq DB pages publicly, replace the Tailscale deployment with the reverse proxy route.

### Useful Links

- [Logseq Open Collective for Official Sync ("Backers" & "Sponsors")](https://opencollective.com/logseq)
- [Logseq DB Documentation](https://github.com/logseq/docs/blob/master/db-version.md)
- [Logseq DB Unofficial FAQ by Luhmann](https://logseq.io/page/e87c7359-51f7-44fe-87b3-4a0cd9f2dee3/695feeec-88be-4c5b-8bf2-572513c2f730)
- [Logseq Discord Server](https://discord.gg/logseq-725182569297215569)
