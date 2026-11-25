"""Measure the exposed walls of a fixed two-circle union from the saved model."""
import json
import math
from pathlib import Path
import sys
from j1_geometry import evaluate, gauss, topology, volume


def measure(model):
    blocks = model["blocks"]
    target_volume = 8*(400-16*math.pi/3-2*math.sqrt(3))
    volumes = {str(q):volume(blocks,q) for q in [8,12,16]}
    nodes, weights=gauss(16)
    wall_areas = [0.,0.]
    radial_error = 0.
    minimum_distance_other = math.inf
    intersections = set()
    for block in blocks:
        center_point,_=evaluate(block,[.5,1.,0.])
        sign = 1 if center_point[0]>0 else -1
        for u in [i/32 for i in range(33)]:
            (x,y,z),_=evaluate(block,[u,1.,0.])
            radial_error=max(radial_error,abs(math.hypot(x-sign,y)-2))
            minimum_distance_other=min(minimum_distance_other,math.hypot(x+sign,y))
            if abs(x)<1e-12:
                intersections.add((round(x,12),round(y,12)))
        for u,wu in zip(nodes,weights):
            _,(du,_,dw)=evaluate(block,[u,1.,.5])
            cross=[du[(i+1)%3]*dw[(i+2)%3]-du[(i+2)%3]*dw[(i+1)%3] for i in range(3)]
            wall_areas[int(sign>0)]+=wu*math.hypot(*cross)
    topo=topology(blocks)
    expected_junctions={(0,round(y,12)) for y in [-math.sqrt(3),math.sqrt(3)]}
    # The dual graph of physical whole-face adjacencies must be connected.
    from j1_geometry import face_samples
    adjacency={i:set() for i in range(len(blocks))}
    seen={}
    for i,block in enumerate(blocks):
        for axis in range(3):
            for side in (0,1):
                key=tuple(sorted(face_samples(block,axis,side)))
                if key in seen:
                    adjacency[i].add(seen[key])
                    adjacency[seen[key]].add(i)
                else:
                    seen[key]=i
    visited=set()
    components=0
    for i in adjacency:
        if i in visited:
            continue
        pending=[i]
        components+=1
        while pending:
            k=pending.pop()
            if k not in visited:
                visited.add(k)
                pending.extend(adjacency[k]-visited)
    target_wall=64*math.pi/3
    assert abs(volumes["16"]/target_volume-1)<1e-10
    assert radial_error<=1e-12 and minimum_distance_other>=2-1e-12
    assert intersections==expected_junctions and components==1
    assert topo == {"interfaces":8,"exterior_faces":32,"nonmanifold_faces":0,"boundary_euler_sampled":0}
    assert all(abs(a/target_wall-1)<1e-11 for a in wall_areas)
    return {"target_volume_mm3":target_volume,"integrated_volume_mm3":volumes,
            "target_exposed_area_per_circle_mm2":target_wall,"wall_areas_mm2":wall_areas,
            "sample_maximum_radial_error_mm":radial_error,"sample_minimum_distance_to_other_center_mm":minimum_distance_other,
            "cross_section_junctions_mm":sorted(intersections),"sampled_components":components,**topo}

if __name__=="__main__":
    print(json.dumps(measure(json.loads(Path(sys.argv[1]).read_text())),indent=2))
