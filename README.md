# Maria's PDF Reduce

**Same quality. Smaller files.**

Live: https://jimbodakilla.github.io/maria-pdf-reduce/

Fast, easy and secure PDF compression for your documents. Phone-scanned PDFs are
usually enormous for one avoidable reason: each page is a photograph, but scanners
often store it with Flate compression, which is lossless and designed for flat
graphics and text rather than photographs. A single page can occupy 9 MB that way.

This re-encodes those page images as JPEG and leaves the rest of the file alone.

- **Resolution is never reduced.** Every image keeps its exact pixel dimensions.
- **Page count, page geometry, text layers and colour profiles are preserved.**
- **Nothing is uploaded.** There is no server; all work happens in the browser.
- **Never returns a file larger than the original.** If it cannot help, it says so.
- **Sends straight from the phone.** Where the browser supports it, each finished file
  gets a **Send** button that opens the native share sheet, so the PDF goes directly to
  WhatsApp, Messages or Mail instead of landing in Files to be hunted down later.
- **Installs to the Home Screen and works offline.** No network is needed to compress.

## Measured results

A batch of 20 phone-scanned documents, 602 MB in total:

| Setting | Result | Worst-page PSNR |
|---|---|---|
| Highest quality | 116 MB (about 5× smaller) | 39.9 dB |
| Best for email | 48 MB (about 12× smaller) | 35.7 dB |
| Smallest | ~30 MB (about 18× smaller) | — |

Roughly 40 dB is the usual threshold for "visually lossless". Best for email sits
just below it: indistinguishable when reading or printing, with mild softening
visible only when pixel-peeping.

JPEG is lossy. Keep your originals if a file may ever need to serve as evidentiary proof.

## Running locally

A single HTML file plus icons and `vendor/`, no build step. Serve the folder:

```sh
python3 -m http.server 8000
```

Two dependencies, `pdf-lib` and `pako`, are vendored into `vendor/` rather than loaded
from a CDN: a blocked or slow network used to leave the page dead, and a phone is exactly
where that happens. Only `pako`'s inflate-only build is shipped, since nothing here
deflates. To refresh them:

```sh
npm pack pdf-lib@1.17.1 pako@2.1.0
tar xzf pdf-lib-1.17.1.tgz && cp package/dist/pdf-lib.min.js vendor/
tar xzf pako-2.1.0.tgz && cp package/dist/pako_inflate.min.js vendor/
```

Sharing needs a secure context, so test over `https://` or `127.0.0.1`; `file://` gives
no service worker and no share sheet.

## Licence

MIT
