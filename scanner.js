/* Document scanner geometry + imaging.
   Pure functions over ImageData so they can be tested headlessly.
   Detection runs on a downscaled copy; the warp runs at full resolution. */
(function (root) {
"use strict";

/* ---------- small helpers ---------- */
function toGrayDownscaled(img, targetLong){
  const { width:W, height:H, data } = img;
  const scale = Math.min(1, targetLong / Math.max(W, H));
  const w = Math.max(1, Math.round(W * scale)), h = Math.max(1, Math.round(H * scale));
  const g = new Float32Array(w * h);
  const bx = W / w, by = H / h;
  for (let y = 0; y < h; y++){
    const y0 = Math.floor(y * by), y1 = Math.min(H, Math.floor((y + 1) * by));
    for (let x = 0; x < w; x++){
      const x0 = Math.floor(x * bx), x1 = Math.min(W, Math.floor((x + 1) * bx));
      let s = 0, n = 0;
      for (let yy = y0; yy < y1; yy++){
        let o = (yy * W + x0) * 4;
        for (let xx = x0; xx < x1; xx++, o += 4){
          s += 0.299 * data[o] + 0.587 * data[o+1] + 0.114 * data[o+2]; n++;
        }
      }
      g[y * w + x] = n ? s / n : 0;
    }
  }
  return { g, w, h, scale };
}

function otsu(g){
  const hist = new Float64Array(256);
  for (let i = 0; i < g.length; i++) hist[Math.max(0, Math.min(255, g[i] | 0))]++;
  const total = g.length;
  let sum = 0; for (let i = 0; i < 256; i++) sum += i * hist[i];
  let sumB = 0, wB = 0, best = 0, thr = 127;
  for (let t = 0; t < 256; t++){
    wB += hist[t]; if (!wB) continue;
    const wF = total - wB; if (!wF) break;
    sumB += t * hist[t];
    const mB = sumB / wB, mF = (sum - sumB) / wF;
    const between = wB * wF * (mB - mF) * (mB - mF);
    if (between > best){ best = between; thr = t; }
  }
  return thr;
}

/* largest bright connected component, as a mask */
function largestBrightBlob(g, w, h, thr){
  const lab = new Int32Array(w * h).fill(-1);
  const stack = new Int32Array(w * h);
  let bestId = -1, bestCount = 0, id = 0;
  for (let i = 0; i < w * h; i++){
    if (lab[i] !== -1 || g[i] < thr) continue;
    let sp = 0, count = 0; stack[sp++] = i; lab[i] = id;
    while (sp){
      const p = stack[--sp]; count++;
      const x = p % w, y = (p / w) | 0;
      if (x > 0     && lab[p-1] === -1 && g[p-1] >= thr){ lab[p-1] = id; stack[sp++] = p-1; }
      if (x < w - 1 && lab[p+1] === -1 && g[p+1] >= thr){ lab[p+1] = id; stack[sp++] = p+1; }
      if (y > 0     && lab[p-w] === -1 && g[p-w] >= thr){ lab[p-w] = id; stack[sp++] = p-w; }
      if (y < h - 1 && lab[p+w] === -1 && g[p+w] >= thr){ lab[p+w] = id; stack[sp++] = p+w; }
    }
    if (count > bestCount){ bestCount = count; bestId = id; }
    id++;
  }
  return { lab, bestId, bestCount };
}

/* boundary points of the blob (any pixel with a non-blob 4-neighbour) */
function blobBoundary(lab, bestId, w, h){
  const pts = [];
  for (let y = 0; y < h; y++){
    for (let x = 0; x < w; x++){
      const i = y * w + x;
      if (lab[i] !== bestId) continue;
      if (x === 0 || y === 0 || x === w-1 || y === h-1 ||
          lab[i-1] !== bestId || lab[i+1] !== bestId ||
          lab[i-w] !== bestId || lab[i+w] !== bestId) pts.push([x, y]);
    }
  }
  return pts;
}

function convexHull(pts){
  if (pts.length < 4) return pts.slice();
  const p = pts.slice().sort((a,b) => a[0] - b[0] || a[1] - b[1]);
  const cross = (o,a,b) => (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0]);
  const lower = [];
  for (const q of p){ while (lower.length >= 2 && cross(lower[lower.length-2], lower[lower.length-1], q) <= 0) lower.pop(); lower.push(q); }
  const upper = [];
  for (let i = p.length - 1; i >= 0; i--){ const q = p[i];
    while (upper.length >= 2 && cross(upper[upper.length-2], upper[upper.length-1], q) <= 0) upper.pop(); upper.push(q); }
  lower.pop(); upper.pop();
  return lower.concat(upper);
}

/* thin a hull down to at most `max` points so the quad search stays cheap */
function thinHull(hull, max){
  if (hull.length <= max) return hull;
  const out = []; const step = hull.length / max;
  for (let i = 0; i < max; i++) out.push(hull[Math.floor(i * step)]);
  return out;
}

/* The document is the maximum-area quadrilateral inscribed in the hull.
   For a rectangle photographed at an angle that is exactly its four corners. */
function maxAreaQuad(hull){
  const n = hull.length;
  if (n < 4) return null;
  const area = (a,b,c) => Math.abs((b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1])) / 2;
  let best = null, bestA = -1;
  for (let i = 0; i < n - 3; i++)
    for (let j = i + 1; j < n - 2; j++)
      for (let k = j + 1; k < n - 1; k++)
        for (let l = k + 1; l < n; l++){
          const A = area(hull[i],hull[j],hull[k]) + area(hull[i],hull[k],hull[l]);
          if (A > bestA){ bestA = A; best = [hull[i],hull[j],hull[k],hull[l]]; }
        }
  return best;
}

function orderCorners(q){
  const cx = (q[0][0]+q[1][0]+q[2][0]+q[3][0]) / 4;
  const cy = (q[0][1]+q[1][1]+q[2][1]+q[3][1]) / 4;
  const s = q.slice().sort((a,b) => Math.atan2(a[1]-cy, a[0]-cx) - Math.atan2(b[1]-cy, b[0]-cx));
  let start = 0, bestSum = Infinity;
  for (let i = 0; i < 4; i++){ const v = s[i][0] + s[i][1]; if (v < bestSum){ bestSum = v; start = i; } }
  return [s[start], s[(start+1)%4], s[(start+2)%4], s[(start+3)%4]];   // TL TR BR BL
}

/* ---------- public: find the page ----------
   One brightness threshold is not enough: white paper on a pale desk, or a shadow
   across a corner, both produce a confident-looking but wrong quad. So generate
   several candidates and keep the one whose outline actually sits on an edge. */
function sobel(g, w, h){
  const m = new Float32Array(w * h);
  for (let y = 1; y < h - 1; y++){
    for (let x = 1; x < w - 1; x++){
      const i = y * w + x;
      const gx = -g[i-w-1] - 2*g[i-1] - g[i+w-1] + g[i-w+1] + 2*g[i+1] + g[i+w+1];
      const gy = -g[i-w-1] - 2*g[i-w] - g[i-w+1] + g[i+w-1] + 2*g[i+w] + g[i+w+1];
      m[i] = Math.sqrt(gx*gx + gy*gy);
    }
  }
  return m;
}
/* mean edge strength under the quad's outline: a real page border scores high,
   an arbitrary line across flat desk scores near zero */
function borderScore(quad, mag, w, h){
  let total = 0, n = 0;
  for (let e = 0; e < 4; e++){
    const a = quad[e], b = quad[(e+1) % 4];
    const len = Math.hypot(b[0]-a[0], b[1]-a[1]);
    const steps = Math.max(8, Math.min(160, Math.round(len)));
    for (let s = 0; s <= steps; s++){
      const f = s / steps;
      const x = Math.round(a[0] + (b[0]-a[0]) * f), y = Math.round(a[1] + (b[1]-a[1]) * f);
      if (x < 1 || y < 1 || x >= w-1 || y >= h-1) { n++; continue; }   // off-frame counts as no edge
      let best = 0;                       // look a little either side of the line
      for (let d = -2; d <= 2; d++){
        const xx = Math.min(w-2, Math.max(1, x + d)), yy = Math.min(h-2, Math.max(1, y + d));
        const v = Math.max(mag[y*w + xx], mag[yy*w + x]);
        if (v > best) best = v;
      }
      total += best; n++;
    }
  }
  return n ? total / n : 0;
}
/* A real page boundary has the page on one side and the surface on the other.
   White paper on a white desk does not, which is exactly when brightness
   thresholding fails, so measure the step across the outline directly. */
function edgeStep(quad, g, w, h){
  const cx = (quad[0][0]+quad[1][0]+quad[2][0]+quad[3][0]) / 4;
  const cy = (quad[0][1]+quad[1][1]+quad[2][1]+quad[3][1]) / 4;
  const at = (x, y) => {
    const xi = Math.min(w-1, Math.max(0, Math.round(x))), yi = Math.min(h-1, Math.max(0, Math.round(y)));
    return g[yi*w + xi];
  };
  let sum = 0, n = 0;
  for (let e = 0; e < 4; e++){
    const a = quad[e], b = quad[(e+1) % 4];
    const len = Math.hypot(b[0]-a[0], b[1]-a[1]);
    const steps = Math.max(6, Math.min(90, Math.round(len / 6)));
    for (let s = 1; s < steps; s++){
      const f = s / steps;
      const x = a[0] + (b[0]-a[0]) * f, y = a[1] + (b[1]-a[1]) * f;
      let nx = x - cx, ny = y - cy;                       // outward direction
      const d = Math.hypot(nx, ny) || 1; nx /= d; ny /= d;
      sum += at(x - nx*7, y - ny*7) - at(x + nx*7, y + ny*7);   // inside minus outside
      n++;
    }
  }
  return n ? sum / n : 0;
}
function candidateFrom(g, w, h, thr){
  const { lab, bestId, bestCount } = largestBrightBlob(g, w, h, thr);
  if (bestId < 0 || bestCount < w * h * 0.04) return null;
  const hull = convexHull(blobBoundary(lab, bestId, w, h));
  const quad = maxAreaQuad(thinHull(hull, 36));
  return quad ? orderCorners(quad) : null;
}
function percentile(g, p){
  const hist = new Float64Array(256);
  for (let i = 0; i < g.length; i++) hist[Math.max(0, Math.min(255, g[i] | 0))]++;
  const want = g.length * p;
  let acc = 0;
  for (let v = 0; v < 256; v++){ acc += hist[v]; if (acc >= want) return v; }
  return 255;
}

function detectCorners(img, opts){
  const targetLong = (opts && opts.targetLong) || 640;
  const { g, w, h, scale } = toGrayDownscaled(img, targetLong);
  const mag = sobel(g, w, h);
  let norm = 0;
  for (let i = 0; i < mag.length; i++) norm += mag[i];
  norm = (norm / mag.length) || 1;

  const thrs = [otsu(g), percentile(g, 0.35), percentile(g, 0.5), percentile(g, 0.62), percentile(g, 0.75)];
  const seen = new Set();
  let best = null, bestScore = 0, bestStep = 0;
  for (const thr of thrs){
    const key = thr | 0;
    if (seen.has(key)) continue;
    seen.add(key);
    const quad = candidateFrom(g, w, h, thr);
    if (!quad) continue;
    const area = quadArea(quad), frame = w * h;
    const cov = area / frame;
    if (cov < 0.05 || cov > 0.985) continue;                  // whole frame, or nothing
    const border = borderScore(quad, mag, w, h) / norm;
    const step = Math.abs(edgeStep(quad, g, w, h));
    let score = border * Math.min(2, 0.35 + step / 14);       // outline on an edge AND a real step across it
    if (cov > 0.93) score *= 0.5;                             // probably swallowed the desk
    if (score > bestScore){ bestScore = score; best = quad; bestStep = step; }
  }
  if (!best) return null;
  const inv = 1 / scale;
  const corners = best.map(p => [p[0] * inv, p[1] * inv]);
  const area = quadArea(corners);
  if (area < img.width * img.height * 0.05) return null;
  return { corners, coverage: area / (img.width * img.height), confidence: bestScore, step: bestStep };
}

function quadArea(q){
  let a = 0;
  for (let i = 0; i < 4; i++){ const p = q[i], n = q[(i+1)%4]; a += p[0]*n[1] - n[0]*p[1]; }
  return Math.abs(a) / 2;
}

/* ---------- homography ---------- */
function solve8(A, b){                       // gaussian elimination, 8x8
  const n = 8;
  for (let i = 0; i < n; i++){
    let piv = i;
    for (let r = i + 1; r < n; r++) if (Math.abs(A[r][i]) > Math.abs(A[piv][i])) piv = r;
    if (Math.abs(A[piv][i]) < 1e-12) return null;
    [A[i], A[piv]] = [A[piv], A[i]]; [b[i], b[piv]] = [b[piv], b[i]];
    for (let r = i + 1; r < n; r++){
      const f = A[r][i] / A[i][i];
      for (let c = i; c < n; c++) A[r][c] -= f * A[i][c];
      b[r] -= f * b[i];
    }
  }
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--){
    let s = b[i];
    for (let c = i + 1; c < n; c++) s -= A[i][c] * x[c];
    x[i] = s / A[i][i];
  }
  return x;
}

/* maps destination (flat page) -> source (photo), which is what sampling needs */
function homographyDstToSrc(dst, src){
  const A = [], b = [];
  for (let i = 0; i < 4; i++){
    const [x, y] = dst[i], [u, v] = src[i];
    A.push([x, y, 1, 0, 0, 0, -u*x, -u*y]); b.push(u);
    A.push([0, 0, 0, x, y, 1, -v*x, -v*y]); b.push(v);
  }
  const h = solve8(A, b);
  return h ? [h[0],h[1],h[2],h[3],h[4],h[5],h[6],h[7],1] : null;
}

/* output size that preserves the document's real aspect, from its edge lengths */
function outputSize(corners, maxLong){
  const d = (a,b) => Math.hypot(a[0]-b[0], a[1]-b[1]);
  const wTop = d(corners[0], corners[1]), wBot = d(corners[3], corners[2]);
  const hL   = d(corners[0], corners[3]), hR   = d(corners[1], corners[2]);
  let w = Math.round(Math.max(wTop, wBot)), h = Math.round(Math.max(hL, hR));
  const long = Math.max(w, h);
  if (maxLong && long > maxLong){ const s = maxLong / long; w = Math.round(w*s); h = Math.round(h*s); }
  return { w: Math.max(1, w), h: Math.max(1, h) };
}

/* ---------- public: flatten the page ---------- */
function warp(img, corners, maxLong){
  const { w:OW, h:OH } = outputSize(corners, maxLong);
  const H = homographyDstToSrc([[0,0],[OW,0],[OW,OH],[0,OH]], corners);
  if (!H) return null;
  const src = img.data, SW = img.width, SH = img.height;
  const out = new Uint8ClampedArray(OW * OH * 4);
  for (let y = 0; y < OH; y++){
    const yh = y + 0.5;
    for (let x = 0; x < OW; x++){
      const xh = x + 0.5;
      const den = H[6]*xh + H[7]*yh + H[8];
      const u = (H[0]*xh + H[1]*yh + H[2]) / den;
      const v = (H[3]*xh + H[4]*yh + H[5]) / den;
      const o = (y * OW + x) * 4;
      if (u < 0 || v < 0 || u > SW - 1 || v > SH - 1){ out[o]=out[o+1]=out[o+2]=255; out[o+3]=255; continue; }
      const x0 = u | 0, y0 = v | 0;
      const x1 = Math.min(SW-1, x0+1), y1 = Math.min(SH-1, y0+1);
      const fx = u - x0, fy = v - y0, fx1 = 1-fx, fy1 = 1-fy;
      const i00 = (y0*SW+x0)*4, i10 = (y0*SW+x1)*4, i01 = (y1*SW+x0)*4, i11 = (y1*SW+x1)*4;
      const w00 = fx1*fy1, w10 = fx*fy1, w01 = fx1*fy, w11 = fx*fy;
      for (let c = 0; c < 3; c++)
        out[o+c] = src[i00+c]*w00 + src[i10+c]*w10 + src[i01+c]*w01 + src[i11+c]*w11;
      out[o+3] = 255;
    }
  }
  return { data: out, width: OW, height: OH };
}

/* ---------- public: make it look like a scan ----------
   Divides out the lighting gradient so paper reads as white everywhere, the way a
   flatbed would, then applies the requested look. */
function enhance(img, mode){
  const { width:W, height:H, data } = img;
  if (mode === "photo") return img;

  // background (illumination) estimate: coarse grid of local high percentiles
  const GX = 12, GY = 16;
  const cw = Math.ceil(W / GX), ch = Math.ceil(H / GY);
  const bg = new Float32Array(GX * GY);
  const bucket = new Uint8Array(cw * ch);
  for (let gy = 0; gy < GY; gy++){
    for (let gx = 0; gx < GX; gx++){
      const x0 = gx*cw, x1 = Math.min(W, x0+cw), y0 = gy*ch, y1 = Math.min(H, y0+ch);
      let n = 0;
      for (let y = y0; y < y1; y++){
        let o = (y*W + x0)*4;
        for (let x = x0; x < x1; x++, o += 4)
          bucket[n++] = (0.299*data[o] + 0.587*data[o+1] + 0.114*data[o+2]) | 0;
      }
      if (!n){ bg[gy*GX+gx] = 255; continue; }
      const arr = Array.prototype.slice.call(bucket.subarray(0, n)).sort((a,b)=>a-b);
      bg[gy*GX+gx] = Math.max(40, arr[Math.floor(n * 0.88)]);       // paper, not ink
    }
  }
  const bgAt = (x, y) => {                                          // bilinear over the grid
    const fx = Math.min(GX-1.001, Math.max(0, x / cw - 0.5));
    const fy = Math.min(GY-1.001, Math.max(0, y / ch - 0.5));
    const x0 = fx|0, y0 = fy|0, dx = fx-x0, dy = fy-y0;
    const a = bg[y0*GX+x0], b = bg[y0*GX+x0+1], c = bg[(y0+1)*GX+x0], d = bg[(y0+1)*GX+x0+1];
    return a*(1-dx)*(1-dy) + b*dx*(1-dy) + c*(1-dx)*dy + d*dx*dy;
  };

  const out = new Uint8ClampedArray(data.length);
  const BLACK = 0.42;                       // below this share of local paper -> full ink
  for (let y = 0; y < H; y++){
    for (let x = 0; x < W; x++){
      const o = (y*W + x)*4, b = bgAt(x, y);
      if (mode === "bw"){
        const l = (0.299*data[o] + 0.587*data[o+1] + 0.114*data[o+2]) / b;
        const v = l < 0.62 ? 0 : 255;
        out[o] = out[o+1] = out[o+2] = v; out[o+3] = 255;
      } else if (mode === "gray"){
        const l = (0.299*data[o] + 0.587*data[o+1] + 0.114*data[o+2]) / b;
        const v = (l - BLACK) / (1 - BLACK) * 255;
        out[o] = out[o+1] = out[o+2] = v; out[o+3] = 255;
      } else {                                                      // colour, flattened
        for (let c = 0; c < 3; c++){
          const l = data[o+c] / b;
          out[o+c] = (l - BLACK) / (1 - BLACK) * 255;
        }
        out[o+3] = 255;
      }
    }
  }
  return { data: out, width: W, height: H };
}

root.Scanner = { detectCorners, warp, enhance, outputSize, quadArea, homographyDstToSrc };
})(typeof window !== "undefined" ? window : globalThis);
