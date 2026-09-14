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
   The lighting across a phone photo is smooth, so model it as a smooth surface
   rather than a grid of local guesses. A grid puts a wrong estimate in every cell
   the subject covers, which shows up as blotches in the finished page. */
function solveN(A, b, n){
  for (let i = 0; i < n; i++){
    let piv = i;
    for (let r = i + 1; r < n; r++) if (Math.abs(A[r][i]) > Math.abs(A[piv][i])) piv = r;
    if (Math.abs(A[piv][i]) < 1e-9) return null;
    [A[i], A[piv]] = [A[piv], A[i]]; [b[i], b[piv]] = [b[piv], b[i]];
    for (let r = i + 1; r < n; r++){
      const f = A[r][i] / A[i][i];
      for (let c = i; c < n; c++) A[r][c] -= f * A[i][c];
      b[r] -= f * b[i];
    }
  }
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--){
    let sum = b[i];
    for (let c = i + 1; c < n; c++) sum -= A[i][c] * x[c];
    x[i] = sum / A[i][i];
  }
  return x;
}
const qbasis = (x, y) => [1, x, y, x*x, x*y, y*y];
const lbasis = (x, y) => [1, x, y];
function fitSurface(pts, n){
  const basis = n === 3 ? lbasis : qbasis;
  const A = [], b = new Float64Array(n);
  for (let i = 0; i < n; i++) A.push(new Float64Array(n));
  for (const p of pts){
    const f = basis(p[0], p[1]);
    for (let i = 0; i < n; i++){
      for (let j = 0; j < n; j++) A[i][j] += f[i] * f[j];
      b[i] += f[i] * p[2];
    }
  }
  return solveN(A, b, n);
}
const evalSurface = (c, x, y) => {
  const f = c.length === 3 ? lbasis(x, y) : qbasis(x, y);
  let v = 0;
  for (let i = 0; i < c.length; i++) v += c[i] * f[i];
  return v;
};
/* fit to the bright side only, so ink and dark subjects do not drag the paper level down */
function illumination(lum, W, H, order){
  const n = order === 1 ? 3 : 6;
  const step = Math.max(1, Math.round(Math.min(W, H) / 200));
  const all = [];
  for (let y = 0; y < H; y += step)
    for (let x = 0; x < W; x += step)
      all.push([(x / W) * 2 - 1, (y / H) * 2 - 1, lum[y * W + x]]);
  let coef = fitSurface(all, n);
  for (let it = 0; it < 3 && coef; it++){
    const next = [];
    for (const p of all) if (p[2] >= evalSurface(coef, p[0], p[1])) next.push(p);
    if (next.length < 60) break;
    const c2 = fitSurface(next, n);
    if (!c2) break;
    coef = c2;
  }
  return coef;
}

/* Note: automatic background removal was tried three ways and removed here.
   Largest-bright-region picks the background, because ink fragments the sheet into
   pieces smaller than the unbroken table. Thresholding after illumination correction
   fails because the fitted surface absorbs the sheet-vs-table difference. Growing the
   background from the frame edge climbs the soft paper edge and eats the sheet.
   Cropping to the four corners removes the surface exactly and is already in the UI. */
/* Background removal: attempted four ways, none reliable enough to ship.
   1) largest bright region  -> picks the background: ink fragments the sheet into
      pieces smaller than the unbroken table.
   2) threshold after quadratic illumination correction -> the fitted surface bends
      around a centred subject and absorbs the sheet-vs-table contrast.
   3) grow background inward from the frame edge -> climbs the soft paper edge and
      eats the sheet.
   4) plane-flattened, two-stage Otsu, holes filled before choosing the region
      (the ordering fix for 1) -> separates on a dark table but erased most of the
      subject on a light one, so it was removed rather than shipped.
   Cropping to the four corners removes the surface exactly, and Apple's scanner
   does not remove irregular backgrounds either: it finds a rectangle and crops. */
function enhance(img, mode){
  const { width:W, height:H, data } = img;
  if (mode === "photo") return img;

  const lum = new Float32Array(W * H);
  for (let i = 0, o = 0; i < W * H; i++, o += 4)
    lum[i] = 0.299*data[o] + 0.587*data[o+1] + 0.114*data[o+2];

  const coef = illumination(lum, W, H, 2);
  const bgAt = coef
    ? (x, y) => Math.max(28, evalSurface(coef, (x / W) * 2 - 1, (y / H) * 2 - 1))
    : () => 200;

  // white balance: what the paper reads per channel once the lighting is divided out
  const gains = [1, 1, 1];
  {
    const N = 20000, stride = Math.max(1, Math.floor(W * H / N));
    const samp = [[], [], []];
    for (let i = 0, o = 0; i < W * H; i += stride, o = i * 4){
      const b = bgAt(i % W, (i / W) | 0);
      const r0 = data[o] / b, g0 = data[o+1] / b, b0 = data[o+2] / b;
      const hi = Math.max(r0, g0, b0), lo = Math.min(r0, g0, b0);
      if (hi < 0.55 || hi - lo > 0.16) continue;   // ink, or a coloured mark; not paper
      samp[0].push(r0); samp[1].push(g0); samp[2].push(b0);
    }
    if (samp[0].length < 200){                     // nothing neutral to measure
      for (let i = 0, o = 0; i < W * H; i += stride, o = i * 4){
        const b = bgAt(i % W, (i / W) | 0);
        for (let c = 0; c < 3; c++) samp[c].push(data[o + c] / b);
      }
    }
    for (let c = 0; c < 3; c++){
      samp[c].sort((a, b) => a - b);
      const p = samp[c][Math.floor(samp[c].length * 0.92)] || 1;
      gains[c] = p > 0.2 ? 1 / p : 1;
    }
  }

  const BLACK = mode === "bw" ? 0 : 0.28;      // gentler than before, so colour is not crushed
  const PAPER = 232;                           // above this a pixel is blank paper, not faint ink
  let smooth = lum;
  if (mode === "bw"){                          // 3x3 mean, so grain cannot decide a pixel
    smooth = new Float32Array(W * H);
    for (let y = 0; y < H; y++){
      const y0 = y > 0 ? y - 1 : 0, y1 = y < H - 1 ? y + 1 : H - 1;
      for (let x = 0; x < W; x++){
        const x0 = x > 0 ? x - 1 : 0, x1 = x < W - 1 ? x + 1 : W - 1;
        let sum = 0, k = 0;
        for (let yy = y0; yy <= y1; yy++)
          for (let xx = x0; xx <= x1; xx++){ sum += lum[yy*W + xx]; k++; }
        smooth[y*W + x] = sum / k;
      }
    }
  }
  const out = new Uint8ClampedArray(data.length);
  for (let y = 0; y < H; y++){
    for (let x = 0; x < W; x++){
      const o = (y * W + x) * 4, b = bgAt(x, y);
      if (mode === "bw"){
        const r = (smooth[y*W + x] / b) * ((gains[0] + gains[1] + gains[2]) / 3);
        const v = r < 0.74 ? 0 : 255;
        out[o] = out[o+1] = out[o+2] = v;
      } else if (mode === "gray"){
        const r = (lum[y*W + x] / b) * ((gains[0] + gains[1] + gains[2]) / 3);
        let v = (r - BLACK) / (1 - BLACK) * 255;
        if (v >= PAPER) v = 255;                   // blank paper
        out[o] = out[o+1] = out[o+2] = v;
      } else {
        let r0 = ((data[o]   / b) * gains[0] - BLACK) / (1 - BLACK) * 255;
        let g0 = ((data[o+1] / b) * gains[1] - BLACK) / (1 - BLACK) * 255;
        let b0 = ((data[o+2] / b) * gains[2] - BLACK) / (1 - BLACK) * 255;
        const hi2 = Math.max(r0, g0, b0), lo2 = Math.min(r0, g0, b0);
        if (lo2 >= PAPER && hi2 - lo2 <= 20){ r0 = g0 = b0 = 255; }   // blank paper
        out[o] = r0; out[o+1] = g0; out[o+2] = b0;
      }
      out[o+3] = 255;
    }
  }
  if (mode === "bw"){
    const src = out.slice();
    for (let y = 1; y < H - 1; y++){
      for (let x = 1; x < W - 1; x++){
        const o = (y*W + x) * 4;
        let dark = 0;
        for (let dy = -1; dy <= 1; dy++)
          for (let dx = -1; dx <= 1; dx++)
            if (src[((y+dy)*W + (x+dx)) * 4] < 128) dark++;
        const self = src[o] < 128;
        if (self && dark <= 1){ out[o] = out[o+1] = out[o+2] = 255; }       // lone black speck
        else if (!self && dark >= 8){ out[o] = out[o+1] = out[o+2] = 0; }   // lone white pinhole
      }
    }
  }
  return { data: out, width: W, height: H };
}

/* ---------- 1-bit PNG ----------
   Canvas only writes 8-bit RGBA, so a pure black-and-white page costs four channels
   per pixel: about 250 KB for a text page. The same bitmap as a 1-bit greyscale PNG
   is about 50 KB. For a 40-page filing that is the difference between 10 MB and 2 MB. */
const CRCT = (() => {
  const t = new Uint32Array(256);
  for (let n = 0; n < 256; n++){
    let c = n;
    for (let k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
    t[n] = c >>> 0;
  }
  return t;
})();
function crc32(bytes){
  let c = 0xFFFFFFFF;
  for (let i = 0; i < bytes.length; i++) c = CRCT[(c ^ bytes[i]) & 0xFF] ^ (c >>> 8);
  return (c ^ 0xFFFFFFFF) >>> 0;
}
function chunk(type, data){
  const out = new Uint8Array(12 + data.length);
  const dv = new DataView(out.buffer);
  dv.setUint32(0, data.length);
  for (let i = 0; i < 4; i++) out[4 + i] = type.charCodeAt(i);
  out.set(data, 8);
  const body = out.subarray(4, 8 + data.length);
  dv.setUint32(8 + data.length, crc32(body));
  return out;
}
/* img must already be black or white; anything >=128 is treated as white */
function packBilevel(img, withFilterByte){
  const { width:W, height:H, data } = img;
  const rowBytes = (W + 7) >> 3;
  const stride = rowBytes + (withFilterByte ? 1 : 0);
  const raw = new Uint8Array(stride * H);
  let p = 0;
  for (let y = 0; y < H; y++){
    if (withFilterByte) raw[p++] = 0;                 // PNG filter: none
    let byte = 0, bit = 7;
    for (let x = 0; x < W; x++){
      if (data[(y * W + x) * 4] >= 128) byte |= (1 << bit);   // 1 = white
      if (bit === 0){ raw[p++] = byte; byte = 0; bit = 7; } else bit--;
    }
    if (bit !== 7) raw[p++] = byte;                   // flush a partial last byte
  }
  return { bytes: raw, width: W, height: H, rowBytes };
}
function png1bit(img){
  if (typeof pako === "undefined") return null;      // fall back to canvas PNG
  const { width:W, height:H } = img;
  const raw = packBilevel(img, true).bytes;
  const idat = pako.deflate(raw, { level: 9 });
  const ihdr = new Uint8Array(13);
  const dv = new DataView(ihdr.buffer);
  dv.setUint32(0, W); dv.setUint32(4, H);
  ihdr[8] = 1;    // bit depth
  ihdr[9] = 0;    // colour type: greyscale
  ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;
  const parts = [new Uint8Array([137,80,78,71,13,10,26,10]),
                 chunk("IHDR", ihdr), chunk("IDAT", idat), chunk("IEND", new Uint8Array(0))];
  let total = 0; for (const a of parts) total += a.length;
  const out = new Uint8Array(total);
  let q = 0; for (const a of parts){ out.set(a, q); q += a.length; }
  return new Blob([out], { type: "image/png" });
}

/* ---------- mixed raster ----------
   A scanned page is crisp black text on a mostly blank sheet, with a little colour.
   Encoding all of that as one JPEG spends almost everything on the text edges.
   Split it: the text becomes a 1-bit stencil at full size, the colour becomes a
   small JPEG underneath. Text ends up sharper than the JPEG version and the page
   costs about a sixth as much. */
function splitLayers(img, inkAt, chromaAt){
  const { width:W, height:H, data } = img;
  const mask = { data: new Uint8ClampedArray(W * H * 4), width: W, height: H };
  const bg   = { data: new Uint8ClampedArray(W * H * 4), width: W, height: H };
  let ink = 0;
  for (let i = 0, o = 0; i < W * H; i++, o += 4){
    const r = data[o], g = data[o+1], b = data[o+2];
    const l = 0.299*r + 0.587*g + 0.114*b;
    const chroma = Math.max(r, g, b) - Math.min(r, g, b);
    const isInk = l < inkAt && chroma < chromaAt;
    mask.data[o] = mask.data[o+1] = mask.data[o+2] = isInk ? 0 : 255;
    mask.data[o+3] = 255;
    if (isInk){ bg.data[o] = bg.data[o+1] = bg.data[o+2] = 255; ink++; }
    else { bg.data[o] = r; bg.data[o+1] = g; bg.data[o+2] = b; }
    bg.data[o+3] = 255;
  }
  return { mask, bg, inkFraction: ink / (W * H) };
}

root.Scanner = { detectCorners, warp, enhance, outputSize, quadArea, homographyDstToSrc, png1bit, packBilevel, splitLayers };
})(typeof window !== "undefined" ? window : globalThis);
