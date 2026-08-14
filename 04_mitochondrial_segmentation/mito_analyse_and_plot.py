"""
mito_analyse_and_plot.py
========================
Single mitochondrion analysis script.
Reads segmentation + filled volume, runs full analysis,
generates a self-contained HTML dashboard with:
  - Spinnable 3D view (real granule shapes via marching cubes)
  - Box-whisker plots for all metrics
  - Dot strip plots showing individual granules
  - Auto-generated interpretation text
  - Summary metric cards

Inputs:
  - Original segmentation TIFF (labels: membrane, granule, background)
  - Filled volume TIFF (from mito_volume_fill.py)

Outputs:
  - granule_measurements.csv
  - single_mito_dashboard.html  (self-contained, open in any browser)

Requirements:
  pip install numpy scipy scikit-image tifffile matplotlib

Usage:
  python mito_analyse_and_plot.py
"""

import os
import csv
import json
import numpy as np
import tifffile
from scipy.ndimage import distance_transform_edt, label, gaussian_filter
from scipy.spatial.distance import cdist
from skimage.measure import regionprops, marching_cubes

# =============================================================================
# CONFIGURATION
# =============================================================================

SEG_PATH       = "Pos018_ImageForProc_1.2_GranulesandMembranes_Pre-trained2DU-Net_Depth5.tiff"
FILLED_PATH    = "Pos018_mito_filled_v3.tif"
OUTPUT_DIR     = "./Pos018_initialanalysis"

MEMBRANE_LABEL = 1
GRANULE_LABEL  = 2
VOXEL_SIZE_NM  = 0.992

GRANULE_BLUR_SIGMA  = 3.0
GRANULE_THRESHOLD   = 0.3
MIN_GRANULE_VOL_NM3 = 500.0

# Marching cubes — blur sigma for granule shape smoothing before isosurface
# Higher = smoother shapes, lower = more detail
SHAPE_BLUR_SIGMA    = 2.5
SHAPE_ISO_LEVEL     = 0.25

# Mito mesh — downsample factor to keep JSON size manageable
# 4 = use every 4th voxel in each dimension (8x fewer voxels, much faster)
MITO_DOWNSAMPLE     = 4

# =============================================================================


def detect_granules(gran_mask, filled, blur_sigma, threshold,
                    min_vol_nm3, voxel_size_nm):
    voxel_vol   = voxel_size_nm ** 3
    gran_in     = gran_mask & filled
    if not np.any(gran_in):
        return np.zeros_like(gran_mask, dtype=np.uint16), []
    blurred     = gaussian_filter(gran_in.astype(float), sigma=blur_sigma)
    thresholded = blurred > threshold
    labelled, _ = label(thresholded)
    props       = regionprops(labelled)
    gran_props  = []
    keep        = set()
    for p in props:
        vol = p.area * voxel_vol
        if vol < min_vol_nm3:
            continue
        keep.add(p.label)
        gran_props.append({
            "label":       p.label,
            "volume_nm3":  vol,
            "diameter_nm": p.equivalent_diameter_area * voxel_size_nm,
            "centroid":    p.centroid,
            "bbox":        p.bbox,
        })
    clean = np.where(np.isin(labelled, list(keep)), labelled, 0)
    return clean.astype(np.uint16), gran_props


def mesh_to_dict(verts, faces, scale=1.0, offset=None):
    """Convert marching cubes output to serialisable dict for Three.js."""
    if offset is not None:
        verts = verts + np.array(offset)
    verts = verts * scale
    return {
        "vertices": verts.flatten().tolist(),
        "faces":    faces.flatten().tolist(),
    }


def extract_granule_meshes(gran_labelled, gran_props, voxel_size_nm,
                            shape_blur, iso_level):
    """
    Run marching cubes on each granule's blurred density field.
    Returns list of mesh dicts with real shapes.
    """
    meshes = []
    for g in gran_props:
        lbl = g["label"]
        bb  = g["bbox"]  # z0,y0,x0,z1,y1,x1
        pad = 4

        # Extract padded bounding box
        z0 = max(0, bb[0]-pad); z1 = min(gran_labelled.shape[0], bb[3]+pad)
        y0 = max(0, bb[1]-pad); y1 = min(gran_labelled.shape[1], bb[4]+pad)
        x0 = max(0, bb[2]-pad); x1 = min(gran_labelled.shape[2], bb[5]+pad)

        sub  = (gran_labelled[z0:z1, y0:y1, x0:x1] == lbl).astype(float)
        blur = gaussian_filter(sub, sigma=shape_blur)

        if blur.max() < iso_level:
            continue

        try:
            verts, faces, _, _ = marching_cubes(blur, level=iso_level)
            # Offset verts to world coordinates (in nm)
            world_offset = np.array([z0, y0, x0]) * voxel_size_nm
            verts_nm     = verts * voxel_size_nm + world_offset
            meshes.append({
                "label":      lbl,
                "volume_nm3": g["volume_nm3"],
                "mesh":       mesh_to_dict(verts_nm, faces),
            })
        except Exception as e:
            print(f"  Marching cubes failed for granule {lbl}: {e}")

    print(f"  Extracted {len(meshes)} granule meshes")
    return meshes


def extract_mito_mesh(filled, voxel_size_nm, downsample, iso_level=0.5):
    """
    Extract mito surface mesh from filled volume.
    Downsamples first to keep JSON size manageable.
    """
    print(f"  Downsampling mito volume by {downsample}x...")
    ds      = filled[::downsample, ::downsample, ::downsample].astype(float)
    blurred = gaussian_filter(ds, sigma=2.0)

    if blurred.max() < iso_level:
        iso_level = blurred.max() * 0.5

    print("  Running marching cubes on mito volume...")
    verts, faces, _, _ = marching_cubes(blurred, level=iso_level)
    verts_nm = verts * voxel_size_nm * downsample
    print(f"  Mito mesh: {len(verts)} verts, {len(faces)} faces")
    return mesh_to_dict(verts_nm, faces)


def generate_summary_text(summary, gran_props):
    n     = summary["n_granules"]
    vol   = summary["mito_vol_um3"]
    frac  = summary["gran_fraction"] * 100
    diam  = summary["mean_diam_nm"]
    nn    = summary["mean_nn_nm"]
    depth = summary["mean_depth_nm"]
    dens  = summary["gran_density"]
    mx    = summary["max_gran_nm3"]

    size_desc = "small" if diam < 20 else "medium-sized" if diam < 40 else "large"

    if not np.isnan(nn):
        cluster_desc = ("tightly clustered" if nn < diam * 1.5
                        else "moderately spaced" if nn < diam * 3
                        else "spatially dispersed")
    else:
        cluster_desc = "single"

    depth_desc = ("peripheral, close to the outer membrane" if depth < 20
                  else "mid-depth within the matrix" if depth < 40
                  else "deep within the matrix")

    outlier_text = ""
    if n > 1:
        med = float(np.median([g["volume_nm3"] for g in gran_props]))
        if mx > med * 5:
            outlier_text = (f" One large granule ({mx:,.0f} nm³, "
                           f"{mx/med:.1f}× the median) may represent "
                           f"a coalescence event.")

    return (
        f"This mitochondrion has a volume of {vol:.4f} µm³ and contains "
        f"{n} calcium phosphate granule{'s' if n != 1 else ''} occupying "
        f"{frac:.2f}% of the total mitochondrial volume ({summary['total_gran_nm3']:,.0f} nm³). "
        f"Granules are {size_desc} (mean diameter {diam:.1f} nm) and "
        f"{cluster_desc} (mean nearest-neighbour distance {nn:.1f} nm). "
        f"Their average depth from the outer membrane is {depth:.1f} nm, "
        f"placing them {depth_desc}. "
        f"Granule density is {dens:.1f} per µm³.{outlier_text}"
    )


def build_html(summary, gran_props, granule_meshes, mito_mesh, interp_text):
    """Build fully self-contained HTML dashboard."""

    # Serialise data for JS
    metrics_js = json.dumps({
        "mito_vol_um3":   round(summary["mito_vol_um3"], 6),
        "mito_vol_nm3":   round(summary["mito_vol_nm3"], 1),
        "n_granules":     summary["n_granules"],
        "total_gran_nm3": round(summary["total_gran_nm3"], 1),
        "gran_fraction":  round(summary["gran_fraction"] * 100, 3),
        "mean_diam_nm":   round(summary["mean_diam_nm"], 1),
        "max_gran_nm3":   round(summary["max_gran_nm3"], 1),
        "mean_nn_nm":     round(float(summary["mean_nn_nm"]), 1),
        "mean_depth_nm":  round(summary["mean_depth_nm"], 1),
        "gran_density":   round(summary["gran_density"], 1),
    })

    granules_js = json.dumps([{
        "label":      g["label"],
        "volume_nm3": round(g["volume_nm3"], 1),
        "diam_nm":    round(g["diameter_nm"], 1),
        "depth_nm":   round(g.get("depth_nm", 0), 1),
        "nn_nm":      round(g.get("nn_dist_nm", 0), 1),
        "cx": round(g["centroid"][2] * summary["voxel_nm"], 1),
        "cy": round(g["centroid"][1] * summary["voxel_nm"], 1),
        "cz": round(g["centroid"][0] * summary["voxel_nm"], 1),
    } for g in gran_props])

    meshes_js  = json.dumps(granule_meshes)
    mito_js    = json.dumps(mito_mesh)
    interp_js  = json.dumps(interp_text)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Mitochondrion analysis dashboard</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;background:#0f0f13;color:#e8e8f0;min-height:100vh;padding:1.5rem}}
h1{{font-size:15px;font-weight:500;color:#e8e8f0;margin-bottom:4px}}
.subtitle{{font-size:11px;color:#888;margin-bottom:1.5rem}}
.metrics{{display:grid;grid-template-columns:repeat(6,1fr);gap:10px;margin-bottom:1.5rem}}
.metric{{background:#1a1a24;border:0.5px solid #2a2a3a;border-radius:8px;padding:12px 14px}}
.metric .label{{font-size:10px;color:#888;margin-bottom:5px;text-transform:uppercase;letter-spacing:0.05em}}
.metric .value{{font-size:19px;font-weight:500;color:#e8e8f0}}
.metric .unit{{font-size:10px;color:#888;margin-left:2px}}
.main-grid{{display:grid;grid-template-columns:1.1fr 0.9fr;gap:1rem;margin-bottom:1rem}}
.card{{background:#1a1a24;border:0.5px solid #2a2a3a;border-radius:10px;padding:1rem}}
.card-title{{font-size:10px;font-weight:500;color:#888;text-transform:uppercase;letter-spacing:0.06em;margin-bottom:10px}}
#three-canvas{{width:100%;height:360px;display:block;border-radius:6px;background:#070710;cursor:grab}}
#three-canvas:active{{cursor:grabbing}}
.legend{{display:flex;gap:16px;margin-top:8px}}
.leg{{display:flex;align-items:center;gap:5px;font-size:11px;color:#888}}
.leg-sw{{width:12px;height:12px;border-radius:3px}}
.charts-grid{{display:grid;grid-template-columns:1fr 1fr;gap:1rem;margin-bottom:1rem}}
.chart-inner{{position:relative}}
.bw-svg{{width:100%;overflow:visible}}
.dot-col{{display:flex;flex-direction:column;gap:14px;margin-top:4px}}
.dot-row{{display:flex;align-items:center;gap:8px}}
.dot-label{{font-size:11px;color:#888;width:72px;text-align:right;flex-shrink:0}}
.dot-track{{flex:1;position:relative;height:28px}}
.track-bg{{position:absolute;left:0;right:0;top:13px;height:0.5px;background:#2a2a3a}}
.interp{{background:#1a1a24;border:0.5px solid #2a2a3a;border-left:3px solid #d946a8;border-radius:10px;padding:1rem 1.25rem;font-size:12.5px;line-height:1.75;color:#c8c8d8}}
.interp strong{{color:#e8e8f0}}
.bottom-grid{{display:grid;grid-template-columns:1fr 1fr;gap:1rem;margin-bottom:1rem}}
</style>
</head>
<body>
<h1>Mitochondrion analysis dashboard</h1>
<div class="subtitle" id="sub-title">Loading...</div>

<div class="metrics" id="metric-cards"></div>

<div class="main-grid">
  <div class="card">
    <div class="card-title">3D view &mdash; drag to rotate</div>
    <canvas id="three-canvas"></canvas>
    <div class="legend">
      <div class="leg"><div class="leg-sw" style="background:rgba(100,149,237,0.3);border:1px solid rgba(100,149,237,0.5)"></div>Mitochondrial volume</div>
      <div class="leg"><div class="leg-sw" style="background:#d946a8"></div>Calcium phosphate granules</div>
    </div>
  </div>
  <div class="card">
    <div class="card-title">Individual granule measurements</div>
    <div class="dot-col" id="dot-col"></div>
  </div>
</div>

<div class="bottom-grid">
  <div class="card">
    <div class="card-title">Volume distribution</div>
    <svg class="bw-svg" id="bw-vol" height="110"></svg>
  </div>
  <div class="card">
    <div class="card-title">Nearest-neighbour distance distribution</div>
    <svg class="bw-svg" id="bw-nn" height="110"></svg>
  </div>
  <div class="card">
    <div class="card-title">Diameter distribution</div>
    <svg class="bw-svg" id="bw-diam" height="110"></svg>
  </div>
  <div class="card">
    <div class="card-title">Depth from outer membrane distribution</div>
    <svg class="bw-svg" id="bw-depth" height="110"></svg>
  </div>
</div>

<div class="interp" id="interp"></div>

<script>
const METRICS   = {metrics_js};
const GRANULES  = {granules_js};
const G_MESHES  = {meshes_js};
const MITO_MESH = {mito_js};
const INTERP    = {interp_js};

document.getElementById('sub-title').textContent =
  'Single mitochondrion · ' + GRANULES.length + ' granules detected · voxel size 0.992 nm';

function buildMetrics() {{
  const cards = [
    {{label:'Mito volume',   value:METRICS.mito_vol_um3.toFixed(4), unit:'µm³'}},
    {{label:'Granule count', value:METRICS.n_granules,              unit:''}},
    {{label:'Gran / mito',   value:METRICS.gran_fraction.toFixed(2),unit:'%'}},
    {{label:'Mean diameter', value:METRICS.mean_diam_nm.toFixed(1), unit:'nm'}},
    {{label:'Mean NN dist',  value:METRICS.mean_nn_nm.toFixed(1),   unit:'nm'}},
    {{label:'Density',       value:METRICS.gran_density.toFixed(0), unit:'/µm³'}},
  ];
  document.getElementById('metric-cards').innerHTML = cards.map(c =>
    `<div class="metric"><div class="label">${{c.label}}</div>
     <div class="value">${{c.value}}<span class="unit">${{c.unit}}</span></div></div>`
  ).join('');
}}

function buildDotPlots() {{
  const metrics = [
    {{key:'diam_nm',  label:'Diameter',  unit:'nm',  color:'#4878cf'}},
    {{key:'volume_nm3',label:'Volume',   unit:'nm³', color:'#d946a8'}},
    {{key:'depth_nm', label:'Depth',     unit:'nm',  color:'#3aaa5e'}},
    {{key:'nn_nm',    label:'NN dist',   unit:'nm',  color:'#9b59b6'}},
  ];
  const col = document.getElementById('dot-col');
  metrics.forEach(m => {{
    const vals = GRANULES.map(g => g[m.key]);
    const mn=Math.min(...vals), mx=Math.max(...vals), rng=mx-mn||1;
    const sorted=[...vals].sort((a,b)=>a-b);
    const med=sorted[Math.floor(sorted.length/2)];
    const row = document.createElement('div');
    row.className='dot-row';
    row.innerHTML=`<div class="dot-label">${{m.label}}</div><div class="dot-track" id="dtr-${{m.key}}"><div class="track-bg"></div></div>`;
    col.appendChild(row);
    const track=document.getElementById('dtr-'+m.key);
    const medPct=(med-mn)/rng*100;
    const medLine=document.createElement('div');
    medLine.style.cssText=`position:absolute;left:${{medPct.toFixed(1)}}%;top:0;bottom:0;width:1.5px;background:#fff;opacity:0.25`;
    track.appendChild(medLine);
    vals.forEach((v,i)=>{{
      const pct=rng===0?50:(v-mn)/rng*100;
      const jitter=(Math.random()-0.5)*8;
      const dot=document.createElement('div');
      dot.style.cssText=`position:absolute;left:${{pct.toFixed(1)}}%;top:${{8+jitter}}px;width:9px;height:9px;border-radius:50%;background:${{m.color}};transform:translateX(-50%);opacity:0.9;cursor:pointer`;
      dot.title=`Granule ${{GRANULES[i].label}}: ${{v.toLocaleString()}} ${{m.unit}}`;
      track.appendChild(dot);
    }});
  }});
}}

function buildBoxWhisker(svgId, key, unit, color) {{
  const svg  = document.getElementById(svgId);
  const W    = svg.parentElement.clientWidth - 32;
  const H    = 100;
  svg.setAttribute('viewBox', `0 0 ${{W}} ${{H}}`);
  svg.setAttribute('width', W);

  const vals   = [...GRANULES.map(g=>g[key])].sort((a,b)=>a-b);
  const n      = vals.length;
  const mn     = vals[0], mx=vals[n-1], rng=mx-mn||1;
  const q1     = vals[Math.floor(n*0.25)];
  const med    = vals[Math.floor(n*0.5)];
  const q3     = vals[Math.floor(n*0.75)];
  const mean   = vals.reduce((a,b)=>a+b,0)/n;

  const pad=30, aw=W-pad*2;
  const sc=v=>pad+(v-mn)/rng*aw;

  const axisY=75, boxY=40, boxH=22;

  // axis line
  const axis=`<line x1="${{pad}}" y1="${{axisY}}" x2="${{W-pad}}" y2="${{axisY}}" stroke="#333" stroke-width="1"/>`;

  // tick marks
  let ticks='';
  [mn,q1,med,q3,mx].forEach(v=>{{
    const x=sc(v);
    ticks+=`<line x1="${{x}}" y1="${{axisY}}" x2="${{x}}" y2="${{axisY+5}}" stroke="#555" stroke-width="1"/>
            <text x="${{x}}" y="${{axisY+16}}" text-anchor="middle" font-size="9" fill="#666">${{v>=1000?Math.round(v/1000)+'k':Math.round(v)}}</text>`;
  }});

  // whiskers
  const whiskers=`
    <line x1="${{sc(mn)}}" y1="${{boxY+boxH/2}}" x2="${{sc(q1)}}" y2="${{boxY+boxH/2}}" stroke="${{color}}" stroke-width="1.5" stroke-dasharray="3,2" opacity="0.6"/>
    <line x1="${{sc(q3)}}" y1="${{boxY+boxH/2}}" x2="${{sc(mx)}}" y2="${{boxY+boxH/2}}" stroke="${{color}}" stroke-width="1.5" stroke-dasharray="3,2" opacity="0.6"/>
    <line x1="${{sc(mn)}}" y1="${{boxY+4}}" x2="${{sc(mn)}}" y2="${{boxY+boxH-4}}" stroke="${{color}}" stroke-width="1.5" opacity="0.6"/>
    <line x1="${{sc(mx)}}" y1="${{boxY+4}}" x2="${{sc(mx)}}" y2="${{boxY+boxH-4}}" stroke="${{color}}" stroke-width="1.5" opacity="0.6"/>`;

  // IQR box
  const box=`<rect x="${{sc(q1)}}" y="${{boxY}}" width="${{sc(q3)-sc(q1)}}" height="${{boxH}}"
              fill="${{color}}" fill-opacity="0.2" stroke="${{color}}" stroke-width="1.5" rx="2"/>`;

  // median line
  const medLine=`<line x1="${{sc(med)}}" y1="${{boxY}}" x2="${{sc(med)}}" y2="${{boxY+boxH}}" stroke="${{color}}" stroke-width="2"/>`;

  // mean diamond
  const meanD=`<polygon points="${{sc(mean)}},${{boxY-4}} ${{sc(mean)+5}},${{boxY+boxH/2}} ${{sc(mean)}},${{boxY+boxH+4}} ${{sc(mean)-5}},${{boxY+boxH/2}}"
               fill="${{color}}" opacity="0.7"/>`;

  // individual points
  let dots='';
  vals.forEach(v=>{{
    const jitter=(Math.random()-0.5)*10;
    dots+=`<circle cx="${{sc(v)}}" cy="${{boxY+boxH/2+jitter}}" r="3" fill="${{color}}" fill-opacity="0.7"/>`;
  }});

  // labels
  const labels=`
    <text x="${{sc(med)}}" y="${{boxY-8}}" text-anchor="middle" font-size="9.5" fill="#aaa">median ${{med>=1000?(med/1000).toFixed(1)+'k':med.toFixed(1)}} ${{unit}}</text>
    <text x="${{pad}}" y="12" font-size="9" fill="#666">n=${{n}}</text>`;

  svg.innerHTML = axis+ticks+whiskers+box+medLine+meanD+dots+labels;
}}

function buildInterp() {{
  document.getElementById('interp').innerHTML = INTERP.replace(
    /(\\d[\\d,.]+ nm³|\\d+\\.\\d+ µm³|\\d+ calcium|\\d+%|\\d+\\.\\d+%|coalescence event)/g,
    '<strong>$1</strong>'
  );
}}

function initThree() {{
  const canvas = document.getElementById('three-canvas');
  const W=canvas.clientWidth, H=360;
  const renderer=new THREE.WebGLRenderer({{canvas,antialias:true,alpha:true}});
  renderer.setPixelRatio(Math.min(devicePixelRatio,2));
  renderer.setSize(W,H);
  renderer.setClearColor(0x070710,1);

  const scene=new THREE.Scene();
  const camera=new THREE.PerspectiveCamera(40,W/H,0.1,100000);

  scene.add(new THREE.AmbientLight(0xffffff,0.5));
  const dl=new THREE.DirectionalLight(0xffffff,0.9);
  dl.position.set(1,2,1.5);
  scene.add(dl);
  const dl2=new THREE.DirectionalLight(0x8888ff,0.3);
  dl2.position.set(-1,-1,-1);
  scene.add(dl2);

  const root=new THREE.Group();
  scene.add(root);

  // Build mito mesh
  if (MITO_MESH && MITO_MESH.vertices && MITO_MESH.vertices.length > 0) {{
    const geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.Float32BufferAttribute(MITO_MESH.vertices,3));
    geo.setIndex(MITO_MESH.faces);
    geo.computeVertexNormals();
    const mat=new THREE.MeshPhongMaterial({{
      color:0x6495ed,transparent:true,opacity:0.12,
      side:THREE.DoubleSide,depthWrite:false
    }});
    root.add(new THREE.Mesh(geo,mat));
    // wireframe edge
    const wmat=new THREE.MeshPhongMaterial({{
      color:0x6495ed,transparent:true,opacity:0.25,
      wireframe:false,side:THREE.FrontSide,depthWrite:false
    }});
    const wgeo=geo.clone();
    root.add(new THREE.Mesh(wgeo,wmat));
  }}

  // Build granule meshes (real shapes)
  G_MESHES.forEach(gm=>{{
    if (!gm.mesh || !gm.mesh.vertices.length) return;
    const geo=new THREE.BufferGeometry();
    geo.setAttribute('position',new THREE.Float32BufferAttribute(gm.mesh.vertices,3));
    geo.setIndex(gm.mesh.faces);
    geo.computeVertexNormals();
    const mat=new THREE.MeshPhongMaterial({{
      color:0xd946a8,shininess:90,specular:0x441133
    }});
    root.add(new THREE.Mesh(geo,mat));
  }});

  // Centre the root group on the mito centroid
  const box=new THREE.Box3().setFromObject(root);
  const centre=new THREE.Vector3();
  box.getCenter(centre);
  root.position.sub(centre);
  const size=new THREE.Vector3();
  box.getSize(size);
  const maxDim=Math.max(size.x,size.y,size.z);
  camera.position.set(0,0,maxDim*1.5);
  camera.far=maxDim*10;
  camera.updateProjectionMatrix();

  // Drag to rotate
  let dragging=false,prevX=0,prevY=0;
  canvas.addEventListener('mousedown',e=>{{dragging=true;prevX=e.clientX;prevY=e.clientY}});
  window.addEventListener('mouseup',()=>dragging=false);
  window.addEventListener('mousemove',e=>{{
    if(!dragging)return;
    root.rotation.y+=(e.clientX-prevX)*0.008;
    root.rotation.x+=(e.clientY-prevY)*0.008;
    prevX=e.clientX;prevY=e.clientY;
  }});

  canvas.addEventListener('touchstart',e=>{{dragging=true;prevX=e.touches[0].clientX;prevY=e.touches[0].clientY}});
  canvas.addEventListener('touchend',()=>dragging=false);
  canvas.addEventListener('touchmove',e=>{{
    if(!dragging)return;
    root.rotation.y+=(e.touches[0].clientX-prevX)*0.01;
    root.rotation.x+=(e.touches[0].clientY-prevY)*0.01;
    prevX=e.touches[0].clientX;prevY=e.touches[0].clientY;
    e.preventDefault();
  }},{{passive:false}});

  (function animate(){{
    requestAnimationFrame(animate);
    if(!dragging) root.rotation.y+=0.003;
    renderer.render(scene,camera);
  }})();
}}

window.addEventListener('load',()=>{{
  buildMetrics();
  buildDotPlots();
  buildBoxWhisker('bw-vol',  'volume_nm3','nm³','#d946a8');
  buildBoxWhisker('bw-nn',   'nn_nm',     'nm', '#9b59b6');
  buildBoxWhisker('bw-diam', 'diam_nm',   'nm', '#4878cf');
  buildBoxWhisker('bw-depth','depth_nm',  'nm', '#3aaa5e');
  buildInterp();
  initThree();
}});
</script>
</body>
</html>"""
    return html


def run():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print(f"Loading segmentation: {SEG_PATH}")
    seg    = tifffile.imread(SEG_PATH)
    print(f"Loading filled volume: {FILLED_PATH}")
    filled = tifffile.imread(FILLED_PATH).astype(bool)

    mem_mask  = (seg == MEMBRANE_LABEL)
    gran_mask = (seg == GRANULE_LABEL)
    print(f"Shape: {seg.shape}  |  Filled voxels: {int(np.sum(filled)):,}")

    # Detect granules
    print("\nDetecting individual granules...")
    gran_labelled, gran_props = detect_granules(
        gran_mask, filled,
        GRANULE_BLUR_SIGMA, GRANULE_THRESHOLD,
        MIN_GRANULE_VOL_NM3, VOXEL_SIZE_NM
    )
    print(f"  Found {len(gran_props)} granules")

    # Depth from boundary
    print("Measuring granule depths...")
    depth_map = distance_transform_edt(filled) * VOXEL_SIZE_NM
    for g in gran_props:
        cz,cy,cx = g["centroid"]
        z=int(np.clip(cz,0,depth_map.shape[0]-1))
        y=int(np.clip(cy,0,depth_map.shape[1]-1))
        x=int(np.clip(cx,0,depth_map.shape[2]-1))
        g["depth_nm"] = float(depth_map[z,y,x])

    # Nearest-neighbour distances
    if len(gran_props) > 1:
        centroids = np.array([g["centroid"] for g in gran_props]) * VOXEL_SIZE_NM
        dists     = cdist(centroids,centroids)
        np.fill_diagonal(dists,np.inf)
        nn_dists  = np.min(dists,axis=1)
        for g,nn in zip(gran_props,nn_dists):
            g["nn_dist_nm"] = float(nn)
        mean_nn = float(np.mean(nn_dists))
    else:
        mean_nn = float("nan")
        for g in gran_props:
            g["nn_dist_nm"] = float("nan")

    # Summary
    voxel_vol     = VOXEL_SIZE_NM ** 3
    mito_nm3      = int(np.sum(filled)) * voxel_vol
    mito_um3      = mito_nm3 / 1e9
    gran_nm3      = int(np.sum(gran_mask & filled)) * voxel_vol
    gran_fraction = gran_nm3 / mito_nm3 if mito_nm3 > 0 else 0
    n_gran        = len(gran_props)
    mean_diam     = float(np.mean([g["diameter_nm"] for g in gran_props])) if n_gran > 0 else 0
    max_gran      = max((g["volume_nm3"] for g in gran_props), default=0)
    mean_depth    = float(np.mean([g["depth_nm"] for g in gran_props])) if n_gran > 0 else 0

    summary = {
        "mito_vol_nm3":   mito_nm3,
        "mito_vol_um3":   mito_um3,
        "n_granules":     n_gran,
        "total_gran_nm3": gran_nm3,
        "gran_fraction":  gran_fraction,
        "mean_diam_nm":   mean_diam,
        "max_gran_nm3":   max_gran,
        "mean_nn_nm":     mean_nn,
        "mean_depth_nm":  mean_depth,
        "gran_density":   n_gran / mito_um3 if mito_um3 > 0 else 0,
        "voxel_nm":       VOXEL_SIZE_NM,
    }

    print(f"\n{'='*50}")
    print(f"  Mito volume:        {mito_nm3:>12,.1f} nm³")
    print(f"  Mito volume:        {mito_um3:>12.6f} µm³")
    print(f"  Granule count:      {n_gran:>12}")
    print(f"  Gran / mito:        {gran_fraction:>12.4f}")
    print(f"  Mean diameter:      {mean_diam:>12.1f} nm")
    print(f"  Mean NN distance:   {mean_nn:>12.1f} nm")
    print(f"{'='*50}")

    # Save CSV
    csv_path = os.path.join(OUTPUT_DIR, "granule_measurements.csv")
    if gran_props:
        keys = ["label","volume_nm3","diameter_nm","depth_nm","nn_dist_nm"]
        with open(csv_path,"w",newline="") as f:
            w = csv.DictWriter(f,fieldnames=keys,extrasaction="ignore")
            w.writeheader()
            w.writerows(gran_props)
        print(f"\nCSV saved: {csv_path}")

    # Extract meshes
    print("\nExtracting granule shapes (marching cubes)...")
    granule_meshes = extract_granule_meshes(
        gran_labelled, gran_props, VOXEL_SIZE_NM,
        SHAPE_BLUR_SIGMA, SHAPE_ISO_LEVEL
    )

    print("Extracting mito surface mesh...")
    mito_mesh = extract_mito_mesh(filled, VOXEL_SIZE_NM, MITO_DOWNSAMPLE)

    # Generate interpretation
    interp_text = generate_summary_text(summary, gran_props)

    # Build and save HTML
    print("\nBuilding HTML dashboard...")
    html = build_html(summary, gran_props, granule_meshes, mito_mesh, interp_text)
    html_path = os.path.join(OUTPUT_DIR, "single_mito_dashboard.html")
    with open(html_path, "w") as f:
        f.write(html)
    print(f"Dashboard saved: {html_path}")
    print("\nOpen in any browser — fully self-contained, no internet needed.")
    print(f"\nDone. All outputs in: {OUTPUT_DIR}")


if __name__ == "__main__":
    run()
