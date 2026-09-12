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

A single HTML file plus two icons, no build step. Serve the folder:

```sh
python3 -m http.server 8000
```

Two external dependencies, both from CDN: `pdf-lib` and `pako`. If they are blocked,
the page says so instead of silently doing nothing.

## Licence

MIT
