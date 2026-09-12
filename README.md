# Prensa

Squeeze heavy scanned PDFs down to size, in the browser.

**Live:** https://jimbodakilla.github.io/prensa/

Phone-scanned documents are usually enormous for one avoidable reason: each page is
a photograph, but scanners often store it with Flate compression, which is lossless
and designed for flat graphics and text rather than photographs. A single page can
occupy 9 MB that way.

Prensa re-encodes those page images as JPEG and leaves the rest of the file alone.

- **Resolution is never reduced.** Every image keeps its exact pixel dimensions.
- **Page count, page geometry, text layers and colour profiles are preserved.**
- **Nothing is uploaded.** There is no server; all work happens in the browser.
- **Never returns a file larger than the original.** If it cannot help, it says so.

## Measured results

A batch of 20 phone-scanned documents, 602 MB in total:

| Setting | Result | Worst-page PSNR |
|---|---|---|
| Archive | 116 MB (19%) | 39.9 dB |
| Balanced | 48 MB (8%) | 35.7 dB |
| Smallest | ~30 MB (5%) | — |

Roughly 40 dB is the usual threshold for "visually lossless". Balanced sits just
below it: indistinguishable when reading or printing, with mild softening visible
only when pixel-peeping.

JPEG is lossy. Keep your originals if a file may ever need to serve as evidentiary proof.

## Running locally

It is a single HTML file with no build step. Open `index.html`, or serve the folder:

```sh
python3 -m http.server 8000
```

Only two external dependencies, both from CDN: `pdf-lib` and `pako`.

## Licence

MIT
