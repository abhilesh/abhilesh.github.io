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
</style>

<div class="row justify-content-center mt-3">
    <div class="col-12 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/self-hosting-logseq-sync-cover.webp" class="img-fluid rounded z-depth-1" zoomable=false %}
    </div>
</div>

[Logseq](https://logseq.com/) has undergone some transformative changes over the past few years, with the [switch to the Database (DB) version](https://logseq.io/p/e3YDyX5AYr) that brings better performance, structured data entry, Real Time Collaboration (RTC) and a host of other features. The move away from storing notes in local markdown files also means that users can no longer rely on solutions like cloud folders (Google Drive, iCloud etc.) or Syncthing for syncing their notes across devices but now have to rely on Logseq's official Sync services to enable cross-device access.

The official Logseq Sync service is paid service that can currently be enabled by contributing to [Logseq's Open Collective](https://opencollective.com/logseq). In the true spirit of Open Source, the developers have gracefully provided the infrastructure to self-host the Logseq Sync and Publish servers on our own hardware and the community has also built Docker containers ([yshalsager's `logseq-selfhost`](https://github.com/yshalsager/logseq-selfhost)) to make the self-hosted setup simpler.

This post walks through two ways of exposing that server:

- **Local access over [Tailscale](https://tailscale.com/)**, the fastest way to get syncing between your own devices, with nothing exposed to the public internet.
- **Remote access over a [reverse proxy](https://www.cloudflare.com/learning/cdn/glossary/reverse-proxy/)**, useful if you want to sync from anywhere, or cannot install Tailscale (on a work computer for instance).

Both approaches sit on top of the same [`logseq-selfhost`](https://github.com/yshalsager/logseq-selfhost) containers, so you can start with the Tailscale route and layer the reverse proxy on later if you decide you need remote access. The optional Web and Publish services are covered at the end, after the basic Sync server is working.

## Common prerequisites

Whichever access route you choose, you will need:

- **Logseq account**: Self-hosting replaces Logseq's hosted Sync and Publish endpoints, but the current clients still require you to sign in with a Logseq account for authentication. You do not need a subscription while using the self-hosted route, only the account.
- **Device capable of running Docker**: A home server, NAS, desktop computer, or VPS.
- **Docker Engine with Compose**: On Linux, follow [Docker's official install guide](https://docs.docker.com/engine/install/); on macOS and Windows, Docker Desktop includes the Compose plugin.
- **Git**: For cloning the [`logseq-selfhost`](https://github.com/yshalsager/logseq-selfhost) repository.
- **Logseq DB clients**: For every device you want to sync.

## <span class="option-icon--tailscale" role="img" aria-label="Tailscale logo"></span> Local access with Tailscale

[Tailscale](https://tailscale.com/) creates a private mesh network (a "tailnet") between your devices, so the Sync server never needs a public IP, a domain, or a certificate: you just point your Logseq clients at the server's Tailscale IP address.

> <p>💡 Ensure <a href="https://tailscale.com/download">Tailscale</a> is installed and signed in on the server <strong>and</strong> on every device you want to sync from.</p>

### Setup

1. **Install Tailscale on the server and on each client device**, then note the Tailscale IP address assigned to the server (visible in the Tailscale admin console, or by running `tailscale ip -4` on the server itself).

<div class="row justify-content-center mt-3">
    <div class="col-sm-8 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/tailscale-admin-console.webp" class="img-fluid rounded z-depth-1" zoomable=true caption="The server's Tailscale IP, as shown in the admin console." %}
    </div>
</div>

2. **Clone the repository and prepare the environment file:**

   ```bash
   git clone https://github.com/yshalsager/logseq-selfhost.git
   cd logseq-selfhost/images/sync
   cp .env.example .env
   ```

3. **Point the Sync server at your Tailscale IP.** Open `.env` and change:

   ```env
   DB_SYNC_BASE_URL=https://sync.example.com
   ```

   to:

   ```env
   DB_SYNC_BASE_URL=http://<your.tailscale.ip>:8787
   ```

4. **Start the server:**

   ```bash
   docker compose pull
   docker compose up -d
   ```

5. **Verify it is running:**

   ```bash
   curl http://localhost:8787/health
   # {"ok":true}
   ```

### Connecting your clients

- **Desktop:** In Logseq, go to **Settings → Advanced → Sync Server URL** and enter `http://YOUR_TAILSCALE_IP:8787`.

<div class="row justify-content-center mt-3">
    <div class="col-sm-8 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/logseq-custom-sync-server-url.webp" class="img-fluid rounded z-depth-1" zoomable=true caption="Setting the custom Sync Server URL in Logseq's desktop settings." %}
    </div>
</div>

- **iPhone and iPad:** Install [Tailscale for iOS](https://tailscale.com/download/ios), sign in to the same tailnet, and allow iOS to add the VPN configuration. Confirm that `http://YOUR_TAILSCALE_IP:8787/health` opens in Safari, then enter `http://YOUR_TAILSCALE_IP:8787` under **Settings → Advanced → Sync Server URL** in Logseq. Tailscale must be connected whenever Logseq needs to reach the server.

- **Android:** Install [Tailscale for Android](https://tailscale.com/download/android), sign in to the same tailnet, and confirm that the health endpoint opens in your browser. You will first need to log in through [test.logseq.com](https://test.logseq.com), then enter the same server URL under the Sync settings.

- On every device, make sure **Use Logseq Sync Beta** is enabled before creating a new graph or syncing an existing one.

> <p>💡 You can optionally enable MagicDNS in Tailscale and use a memorable address such as <code>http://my-server:8787</code> instead of the server's Tailscale IP.</p>

## <img src="/assets/img/posts/ditch-the-cloud/nginxproxymanager.svg" width="24" height="24" style="margin-right: 6px; vertical-align: middle;"> Remote access with a reverse proxy

If you want to sync from outside your tailnet, or cannot install Tailscale on a device, you can place the same Sync server behind **NGINX Proxy Manager (NPM)** and access it through a normal HTTPS domain.

> <p>💡 I chose NGINX Proxy Manager (NPM) as I already have a server stack running with NPM to connect to. NPM can be replaced with your reverse-proxy manager of choice.</p>

### Prerequisites

- **NGINX Proxy Manager already deployed**, with its own Docker network you can attach other stacks to.
- A domain with a DNS provider you control (this example uses Cloudflare).
- A method for obtaining TLS certificates.
- OpenSSL, for generating a random admin token.

### Setup

1. **Create a deployment directory** containing `docker-compose.yml` and `.env`, with persistent data stored in a `./data` subdirectory:

   ```bash
   mkdir -p ~/home-server/logseq-sync/data
   cd ~/home-server/logseq-sync
   ```

2. **Create the environment file (`.env`):**

   ```env
   DB_SYNC_BASE_URL=https://logseq-sync.example.com
   DB_SYNC_ADMIN_TOKEN=replace-with-random-token
   ```

   Generate the token with:

   ```bash
   openssl rand -hex 32
   ```

3. **Create `docker-compose.yml`:**

<details markdown="1">
<summary><strong>Show the Sync-only Docker Compose file</strong></summary>

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

</details>

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

   Under **Advanced**, add:

   ```nginx
   client_max_body_size 1024m;
   proxy_read_timeout 3600s;
   proxy_send_timeout 3600s;
   ```

   Enable **Websockets Support**, **Block Common Exploits**, **Force SSL**, and **HTTP/2 Support**, then request or select an SSL certificate.

<div class="row justify-content-center mt-3">
    <div class="col-sm-8 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/npm-proxy-host-advanced.webp" class="img-fluid rounded z-depth-1" zoomable=true caption="The Advanced tab for the Sync proxy host, with the custom NGINX configuration and Websockets Support enabled." %}
    </div>
</div>

<aside>
  <p>💡 <strong>Tip:</strong> The longer <code>proxy_read_timeout</code>/<code>proxy_send_timeout</code> values matter here — Logseq Sync operations on large graphs can take a while, and NPM's defaults will otherwise cut the connection early.</p>
</aside>

### Verification

```bash
docker compose ps

curl -i https://logseq-sync.example.com/health
# {"ok":true}
```

### Connecting your clients

In Logseq, go to **Settings → Advanced** and set:

- **Sync Server URL:** `https://logseq-sync.example.com`

<div class="row justify-content-center mt-3">
    <div class="col-sm-8 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/self-hosting-logseq-sync/logseq-advanced-sync-publish-urls.webp" class="img-fluid rounded z-depth-1" zoomable=true caption="The Sync Server URL under Settings → Advanced." %}
    </div>
</div>

### A note on security

The write-up highlights a few security defaults worth keeping when you adapt the Compose file:

- The container runs with a **read-only filesystem** and **`CAP_DROP: ALL`**, with a dedicated `tmpfs` mount for `/tmp`.
- **`no-new-privileges`** is enabled.
- The server is reached through NPM's Docker network, so no container port is published directly on the host.

The Compose file also wires up Cognito-based authentication using Logseq's current authentication service. Review the configuration before exposing the server to the public internet, and keep the contents of `.env` out of version control.

## Optional: Web and Publish servers

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

> <p>💡 While the Publish server can also be exposed over Tailscale, the published pages will only be accessible on devices connected to your tailnet.</p>

## Which one should you use?

- **Tailscale** gets you syncing between your own devices in a few minutes, with nothing exposed publicly — no domain, no certificates, no reverse proxy to maintain. Its limitation is exactly that: only devices on your tailnet can reach it.
- **Reverse proxy** is the better fit if you need to sync from a device you do not control (or cannot install Tailscale on), or if you want to use Logseq's Publish feature to share a graph at a real URL. It asks more of you upfront — a domain, NGINX Proxy Manager, and certificates — but is not limited to your own devices.

Nothing stops you running both against the same underlying data: Tailscale for your own day-to-day syncing, and the proxied Publish endpoint purely for sharing a published graph.

## Sources

This post closely follows two guides — go to them directly for the full, unabridged setup:

- Adri Sood, ["Self-Hosting Logseq Sync Locally"](https://adri99.net/self-hosting-logseq-sync-locally/)
- ["Self-Hosted Logseq Sync Server"](https://logseq-publish.biophiles.lol/p/fmT_7lx_Kz), published with Logseq Publish
- [yshalsager/logseq-selfhost](https://github.com/yshalsager/logseq-selfhost) — the underlying project both guides build on
