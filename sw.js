/* The app is a handful of static files, so the cache only has to answer two questions:
   is the shell still there when the signal is not, and does a new deploy actually land.
   index.html is network-first so a published fix is never hidden behind a stale cache;
   everything else is cache-first, because those files are versioned by CACHE. */
const CACHE = "pdf-reduce-v2";
const SHELL = [
  "./", "./index.html", "./manifest.json",
  "./vendor/pdf-lib.min.js", "./vendor/pako_inflate.min.js",
  "./icon-192.png", "./icon-512.png", "./logo-mark.png"
];

self.addEventListener("install", e => {
  e.waitUntil(caches.open(CACHE).then(c => c.addAll(SHELL)).then(() => self.skipWaiting()));
});

self.addEventListener("activate", e => {
  e.waitUntil(caches.keys()
    .then(keys => Promise.all(keys.filter(k => k !== CACHE).map(k => caches.delete(k))))
    .then(() => self.clients.claim()));
});

self.addEventListener("fetch", e => {
  const req = e.request;
  if (req.method !== "GET" || new URL(req.url).origin !== self.location.origin) return;

  const isPage = req.mode === "navigate" || /\/(index\.html)?$/.test(new URL(req.url).pathname);
  if (isPage){
    e.respondWith(
      fetch(req)
        .then(res => { const copy = res.clone(); caches.open(CACHE).then(c => c.put(req, copy)); return res; })
        .catch(() => caches.match(req).then(hit => hit || caches.match("./index.html")))
    );
    return;
  }
  e.respondWith(caches.match(req).then(hit => hit || fetch(req).then(res => {
    const copy = res.clone();
    caches.open(CACHE).then(c => c.put(req, copy));
    return res;
  })));
});
