"""Independent OCC/Netgen C3D10 and CalculiX reference for the bracket."""

import FreeCAD as App
import Part
import Fem  # Registers the native Netgen document object.

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import re
import resource
import subprocess
import time

V = App.Vector


def shape():
    outline = [(0, 0), (60, 0), (60, 8), (14, 8)]
    edges = [Part.makeLine(V(*a, 0), V(*b, 0)) for a, b in zip(outline, outline[1:])]
    mid = 14 - 6 / math.sqrt(2)
    edges.append(Part.Arc(V(14,8,0), V(mid,mid,0), V(8,14,0)).toShape())
    end = [(8,14), (8,50), (0,50), (0,0)]
    edges.extend(Part.makeLine(V(*a,0), V(*b,0)) for a,b in zip(end,end[1:]))
    body = Part.Face(Part.Wire(edges)).extrude(V(0,0,40))
    body = body.cut(Part.makeCylinder(4,10,V(40,-1,20),V(0,1,0)))
    body = body.cut(Part.makeCylinder(4,10,V(-1,32,20),V(1,0,0)))
    assert body.isValid() and len(body.Solids)==1
    assert abs(body.Volume - (34080-616*math.pi)) < 1e-7
    return body


def load_and_probe(mesh, nodes):
    """Integrate uniform traction on planar quadratic triangles."""
    weights = defaultdict(float)
    interpolation = None
    area = 0.0
    for face in mesh.Faces:
        ids = mesh.getElementNodes(face)
        if not all(abs(nodes[i].x-60) < 1e-8 for i in ids):
            continue
        assert len(ids)==6
        a,b,c = [nodes[i] for i in ids[:3]]
        for (i,j),node in zip(((a,b),(b,c),(c,a)),ids[3:]):
            assert (nodes[node]-(i+j)*0.5).Length < 1e-7
        tri = (b-a).cross(c-a).Length/2
        area += tri
        for i in ids[3:]:
            weights[i] += tri/3
        # On the planar loaded end, invert the affine y-z triangle.
        by,bz,cy,cz = b.y-a.y,b.z-a.z,c.y-a.y,c.z-a.z
        dy,dz=4-a.y,20-a.z
        denom=by*cz-bz*cy
        l2,l3=(dy*cz-dz*cy)/denom,(by*dz-bz*dy)/denom
        l1=1-l2-l3
        if min(l1,l2,l3) >= -1e-9:
            interpolation = dict(zip(ids,(l1*(2*l1-1),l2*(2*l2-1),l3*(2*l3-1),
                                          4*l1*l2,4*l2*l3,4*l3*l1)))
    assert abs(area-320)<1e-7 and interpolation is not None
    return weights, interpolation, area


def block_ids(ids):
    ids = list(sorted(ids))
    return "\n".join(",".join(map(str,ids[i:i+12])) for i in range(0,len(ids),12))


def node_report(dat, name, set_name):
    rows = {}
    active = False
    for line in dat.splitlines():
        if f"for set {set_name} " in line:
            active = name in line
            continue
        if active:
            match = re.fullmatch(r"\s*(\d+)\s+([-+.\dEe]+)\s+([-+.\dEe]+)\s+([-+.\dEe]+)\s*",line)
            if match:
                rows[int(match[1])] = list(map(float,match.groups()[1:]))
            elif line.strip() and rows:
                break
    assert rows, (name,set_name)
    return rows


def run(args):
    root = args.directory.resolve()
    root.mkdir(parents=True, exist_ok=True)
    begin = time.monotonic()
    doc=App.newDocument("Reference")
    feature=doc.addObject("PartDesign::Feature","ReferenceShape")
    feature.Shape=shape()
    mesh_obj=doc.addObject("Fem::FemMeshShapeNetgenObject","Tet10")
    mesh_obj.Shape=feature
    mesh_obj.MaxSize=args.size
    mesh_obj.SecondOrder=True
    doc.recompute()
    mesh=mesh_obj.FemMesh
    assert mesh.VolumeCount and len(mesh.getElementNodes(mesh.Volumes[0]))==10
    nodes=mesh.Nodes
    clamp={i for i,p in nodes.items() if abs(p.y-50)<1e-8}
    tip={i for i,p in nodes.items() if abs(p.x-60)<1e-8}
    weights, probe, area = load_and_probe(mesh,nodes)
    deck=root/"bracket.inp"
    mesh.writeABAQUS(str(deck),1,False)
    with deck.open("a") as stream:
        stream.write("\n*NSET,NSET=CLAMP\n"+block_ids(clamp)+"\n*NSET,NSET=TIP\n"+block_ids(tip)+"\n")
        stream.write("*MATERIAL,NAME=STEEL\n*ELASTIC\n210000,0.3\n*SOLID SECTION,ELSET=Eall,MATERIAL=STEEL\n")
        stream.write("*STEP\n*STATIC\n*BOUNDARY\nCLAMP,1,3\n*CLOAD\n")
        for i,w in sorted(weights.items()):
            stream.write(f"{i},2,{-1000*w/area:.17g}\n")
        stream.write("*NODE PRINT,NSET=TIP\nU\n*NODE PRINT,NSET=CLAMP\nRF\n*END STEP\n")
    mesh_seconds=time.monotonic()-begin
    with (root/"ccx.log").open("w") as log:
        subprocess.run([args.ccx,"-i","bracket"],cwd=root,stdout=log,stderr=subprocess.STDOUT,check=True)
    output=(root/"bracket.dat").read_text()
    disp=node_report(output,"displacements","TIP")
    react=node_report(output,"forces","CLAMP")
    reaction=[sum(v[k] for v in react.values()) for k in range(3)]
    assert math.dist(reaction,[0,1000,0])<0.01, reaction
    report=dict(size_mm=args.size,occ_volume_mm3=feature.Shape.Volume,nodes=len(nodes),
                tet10=mesh.VolumeCount,loaded_area_mm2=area,reaction_n=reaction,
                mean_tip_y_mm=sum(w*disp[i][1] for i,w in weights.items())/area,
                point_tip_y_mm=sum(w*disp[i][1] for i,w in probe.items()),
                geometry_mesh_deck_s=mesh_seconds,solve_s=time.monotonic()-begin-mesh_seconds,
                freecad_version=App.Version()[:3],
                solver_peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                mesh_process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (root/"result.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--size",type=float,required=True)
    parser.add_argument("--ccx",required=True)
    parser.add_argument("directory",type=Path)
    run(parser.parse_args())
