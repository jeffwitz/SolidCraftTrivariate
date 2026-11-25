//! Independent exact coverage and incidence audit for polygonal prisms.
use crate::{
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    incidence::{enumerate_faces,Incidence},
    polygon_prism::{PolygonPrism,Vertex},
};
use std::collections::{BTreeMap,BTreeSet};

#[derive(Debug,Clone,PartialEq)]
pub struct PolygonCertificate{
    pub certificate:BracketCertificate,
    pub volume_sixths_mm3:i64,
    pub projected_area_twice_mm2:i64,
}

#[derive(Clone,Copy)]
struct Tri { ids:[usize;3] }

fn cross(a:[i32;2],b:[i32;2],c:[i32;2])->i64{
    i64::from(b[0]-a[0])*i64::from(c[1]-a[1])-i64::from(b[1]-a[1])*i64::from(c[0]-a[0])
}
fn on_segment(a:[i32;2],b:[i32;2],p:[i32;2])->bool{
    cross(a,b,p)==0&&(0..2).all(|d|(a[d].min(b[d])..=a[d].max(b[d])).contains(&p[d]))
}
fn segments_contact(a:[i32;2],b:[i32;2],c:[i32;2],d:[i32;2])->bool{
    let (x,y,z,w)=(cross(a,b,c),cross(a,b,d),cross(c,d,a),cross(c,d,b));
    (x.signum()!=y.signum()&&z.signum()!=w.signum())||on_segment(a,b,c)||on_segment(a,b,d)||on_segment(c,d,a)||on_segment(c,d,b)
}
fn inside_triangle(p:[i32;2],t:[[i32;2];3])->bool{cross(t[0],t[1],p)>0&&cross(t[1],t[2],p)>0&&cross(t[2],t[0],p)>0}
fn validate_contour(recipe:&PolygonPrism)->Result<(Vec<usize>,i64),&'static str>{
    let n=recipe.vertices.len();
    if !(3..=24).contains(&n)||!(1..=64).contains(&recipe.height_mm)||recipe.top_shift_mm.iter().any(|x|!(-32..=32).contains(x)){return Err("polygon_prism_range")}
    let mut ids=BTreeMap::new();let mut positions=BTreeSet::new();
    for (i,v) in recipe.vertices.iter().enumerate(){
        if v.id.is_empty()||v.id.len()>24||!v.id.bytes().all(|b|b.is_ascii_alphanumeric()||b==b'_')||ids.insert(v.id.as_str(),i).is_some(){return Err("polygon_vertex_id")}
        if v.xy_mm.iter().any(|x|!(-64..=64).contains(x))||!positions.insert(v.xy_mm){return Err("polygon_vertex_coordinates")}
    }
    let mut area2=0_i64;
    for i in 0..n{
        let a=recipe.vertices[i].xy_mm;let b=recipe.vertices[(i+1)%n].xy_mm;let c=recipe.vertices[(i+2)%n].xy_mm;
        if cross(a,b,c)==0{return Err("polygon_collinear_adjacent_edges")}
        area2+=i64::from(a[0])*i64::from(b[1])-i64::from(a[1])*i64::from(b[0]);
        for j in i+1..n{
            if j==i+1||(i==0&&j==n-1){continue}
            let c=recipe.vertices[j].xy_mm;let d=recipe.vertices[(j+1)%n].xy_mm;
            if segments_contact(a,b,c,d){return Err("polygon_nonadjacent_edge_contact")}
        }
    }
    if area2==0{return Err("polygon_zero_area")}
    let mut order=(0..n).collect::<Vec<_>>();if area2<0{order.reverse();area2 = -area2;}
    let first=order.iter().enumerate().min_by_key(|(_,i)|&recipe.vertices[**i].id).unwrap().0;
    order.rotate_left(first);
    Ok((order,area2))
}
fn parse_triangle_id(id:&str,ids:&BTreeMap<&str,usize>)->Result<([usize;3],usize),&'static str>{
    let fields=id.split('/').collect::<Vec<_>>();if fields.len()!=5||fields[0]!="tri"{return Err("polygon_block_id")}
    let a=*ids.get(fields[1]).ok_or("polygon_block_vertex")?;let b=*ids.get(fields[2]).ok_or("polygon_block_vertex")?;let c=*ids.get(fields[3]).ok_or("polygon_block_vertex")?;let at=*ids.get(fields[4]).ok_or("polygon_block_vertex")?;
    if a==b||b==c||c==a{return Err("polygon_block_id")}
    Ok(([a,b,c],at))
}
fn fields_for_index<'a>(ids:&'a BTreeMap<&'a str,usize>,index:usize)->&'a str{ids.iter().find(|(_,i)|**i==index).unwrap().0}
fn expected_controls(tri:[usize;3],at:usize,vertices:&[Vertex],height:i32,shift:[i32;2])->Option<Vec<[i64;4]>>{
    let place=tri.iter().position(|x|*x==at)?;let v=tri[place];let prev=tri[(place+2)%3];let next=tri[(place+1)%3];
    let a=vertices[v].xy_mm;let p=vertices[prev].xy_mm;let n=vertices[next].xy_mm;
    let centroid=[2_i64*(i64::from(a[0])+i64::from(p[0])+i64::from(n[0])),2_i64*(i64::from(a[1])+i64::from(p[1])+i64::from(n[1]))];
    let points=[
        [6*i64::from(a[0]),6*i64::from(a[1])],
        [3*i64::from(a[0]+p[0]),3*i64::from(a[1]+p[1])],
        [3*i64::from(a[0]+n[0]),3*i64::from(a[1]+n[1])],
        centroid,
    ];
    Some(points.into_iter().flat_map(|xy|{
        let bottom=[xy[0],xy[1],0,6];
        let top=[xy[0]+6*i64::from(shift[0]),xy[1]+6*i64::from(shift[1]),6*i64::from(height),6];
        [bottom,top]
    }).collect())
}
fn block_corner_xy(block:&BracketBlock,index:usize)->[i64;2]{[block.controls[2*index][0],block.controls[2*index][1]]}
fn check_quad_jacobian(block:&BracketBlock)->Result<i128,&'static str>{
    let p00=block_corner_xy(block,0);let p01=block_corner_xy(block,1);let p10=block_corner_xy(block,2);let p11=block_corner_xy(block,3);
    let a=[p10[0]-p00[0],p10[1]-p00[1]];let b=[p01[0]-p00[0],p01[1]-p00[1]];
    let c=[p11[0]-p10[0]-p01[0]+p00[0],p11[1]-p10[1]-p01[1]+p00[1]];
    let det=|x:[i64;2],y:[i64;2]|i128::from(x[0])*i128::from(y[1])-i128::from(x[1])*i128::from(y[0]);
    let values=[det(a,b),det(a,[b[0]+c[0],b[1]+c[1]]),det([a[0]+c[0],a[1]+c[1]],b),det([a[0]+c[0],a[1]+c[1]],[b[0]+c[0],b[1]+c[1]])];
    let min=*values.iter().min().unwrap();if min<=0{return Err("polygon_nonpositive_jacobian")};Ok(min)
}
fn edge_key(a:usize,b:usize)->(usize,usize){if a<b{(a,b)}else{(b,a)}}
fn validate_triangle_intersections(tris:&[Tri],vertices:&[Vertex])->Result<(),&'static str>{
    for i in 0..tris.len(){for j in i+1..tris.len(){
        let a=tris[i].ids;let b=tris[j].ids;
        for x in 0..3{for y in 0..3{
            let (ai,an)=(a[x],a[(x+1)%3]);let (bi,bn)=(b[y],b[(y+1)%3]);
            if edge_key(ai,an)==edge_key(bi,bn){continue}
            let (pa,pan,pb,pbn)=(vertices[ai].xy_mm,vertices[an].xy_mm,vertices[bi].xy_mm,vertices[bn].xy_mm);
            let proper=cross(pa,pan,pb).signum()*cross(pa,pan,pbn).signum()<0
                &&cross(pb,pbn,pa).signum()*cross(pb,pbn,pan).signum()<0;
            if proper{return Err("polygon_triangle_overlap")}
            for (point, endpoints) in [(pa,(pb,pbn)),(pan,(pb,pbn)),(pb,(pa,pan)),(pbn,(pa,pan))]{
                if on_segment(endpoints.0,endpoints.1,point)&&point!=endpoints.0&&point!=endpoints.1{
                    return Err("polygon_triangle_t_junction")
                }
            }
        }}
        for &v in &a{if inside_triangle(vertices[v].xy_mm,b.map(|x|vertices[x].xy_mm)){return Err("polygon_triangle_overlap")}}
        for &v in &b{if inside_triangle(vertices[v].xy_mm,a.map(|x|vertices[x].xy_mm)){return Err("polygon_triangle_overlap")}}
    }}
    Ok(())
}

/// Independently reconstruct triangle coverage and audit all quadrilateral prism blocks.
pub fn audit_polygon_prism(recipe:&PolygonPrism,blocks:&[BracketBlock])->Result<PolygonCertificate,&'static str>{
    let (contour,area2)=validate_contour(recipe)?;
    let idmap=recipe.vertices.iter().enumerate().map(|(i,v)|(v.id.as_str(),i)).collect::<BTreeMap<_,_>>();
    if blocks.len()!=3*(recipe.vertices.len()-2){return Err("polygon_block_count")}
    let mut unique=BTreeSet::new();let mut groups=BTreeMap::<[usize;3],BTreeMap<usize,&BracketBlock>>::new();let mut min_det=i128::MAX;
    for block in blocks{
        if !unique.insert(&block.id){return Err("polygon_duplicate_block")}
        if block.degrees!=[1,1,1]||block.controls.len()!=8{return Err("polygon_altered_block")}
        let (triangle,at)=parse_triangle_id(&block.id,&idmap)?;
        if block.id!=format!("tri/{}/{}/{}/{}",fields_for_index(&idmap,triangle[0]),fields_for_index(&idmap,triangle[1]),fields_for_index(&idmap,triangle[2]),fields_for_index(&idmap,at))
            ||triangle.iter().map(|i|fields_for_index(&idmap,*i)).min().unwrap()!=fields_for_index(&idmap,triangle[0])
        {return Err("polygon_block_id")}
        let a=recipe.vertices[triangle[0]].xy_mm;let b=recipe.vertices[triangle[1]].xy_mm;let c=recipe.vertices[triangle[2]].xy_mm;
        let area_tri=cross(a,b,c);if area_tri<=0{return Err("polygon_triangle_orientation")}
        let expected=expected_controls(triangle,at,&recipe.vertices,recipe.height_mm,recipe.top_shift_mm).ok_or("polygon_block_id")?;
        if block.controls!=expected{return Err("polygon_altered_block")}
        min_det=min_det.min(check_quad_jacobian(block)?);
        if groups.entry(triangle).or_default().insert(at,block).is_some(){return Err("polygon_duplicate_corner_patch")}
    }
    let mut triangles=Vec::new();let mut triangle_area_sum=0_i64;let mut edge_uses=BTreeMap::<(usize,usize),Vec<(usize,usize)>>::new();
    for (ids,patches) in groups{
        if patches.len()!=3||ids.iter().any(|id|!patches.contains_key(id)){return Err("polygon_missing_corner_patch")}
        let area=cross(recipe.vertices[ids[0]].xy_mm,recipe.vertices[ids[1]].xy_mm,recipe.vertices[ids[2]].xy_mm);
        if area<=0{return Err("polygon_triangle_orientation")}
        triangle_area_sum=triangle_area_sum.checked_add(area).ok_or("polygon_area_overflow")?;
        for k in 0..3{let a=ids[k];let b=ids[(k+1)%3];edge_uses.entry(edge_key(a,b)).or_default().push((a,b));}
        triangles.push(Tri{ids});
    }
    if triangles.len()!=recipe.vertices.len()-2||triangle_area_sum!=area2{return Err("polygon_section_coverage")}
    validate_triangle_intersections(&triangles,&recipe.vertices)?;
    let boundary=(0..contour.len()).map(|i|edge_key(contour[i],contour[(i+1)%contour.len()])).collect::<BTreeSet<_>>();
    for (edge,uses) in &edge_uses{
        if boundary.contains(edge){if uses.len()!=1||uses[0] != contour.iter().copied().zip(contour.iter().copied().cycle().skip(1)).find(|(a,b)|edge_key(*a,*b)==*edge).unwrap(){return Err("polygon_boundary_incidence")}}
        else if uses.len()!=2||uses[0].0!=uses[1].1||uses[0].1!=uses[1].0{return Err("polygon_internal_incidence")}
    }
    if edge_uses.keys().filter(|e|boundary.contains(e)).count()!=boundary.len(){return Err("polygon_section_coverage")}

    let incidence:Incidence=enumerate_faces(blocks)?;
    let vertex_id=recipe.vertices.iter().map(|v|(v.xy_mm,v.id.as_str())).collect::<BTreeMap<_,_>>();
    let mut design_faces=Vec::with_capacity(incidence.exterior.len());
    for side in &incidence.exterior{
        let parts=side.block.split('/').collect::<Vec<_>>();if parts.len()!=5||parts[0]!="tri"{return Err("polygon_boundary_provenance")}
        let tri=[*idmap.get(parts[1]).ok_or("polygon_boundary_provenance")?,*idmap.get(parts[2]).ok_or("polygon_boundary_provenance")?,*idmap.get(parts[3]).ok_or("polygon_boundary_provenance")?];
        let at=*idmap.get(parts[4]).ok_or("polygon_boundary_provenance")?;let place=tri.iter().position(|x|*x==at).ok_or("polygon_boundary_provenance")?;
        let prev=tri[(place+2)%3];let next=tri[(place+1)%3];
        let tag=match (side.axis,side.upper){
            (2,false)=>"start".to_owned(),(2,true)=>"end".to_owned(),
            (0,false)=>edge_tag(prev,at,&boundary,&recipe.vertices,&vertex_id)?,
            (1,false)=>edge_tag(at,next,&boundary,&recipe.vertices,&vertex_id)?,
            _=>return Err("polygon_internal_face_exterior"),
        };
        design_faces.push((side.clone(),tag));
    }
    design_faces.sort_by(|a,b|a.0.cmp(&b.0));
    let volume=area2.checked_mul(3).and_then(|x|x.checked_mul(i64::from(recipe.height_mm))).ok_or("polygon_volume_overflow")?;
    let det_num=min_det.checked_mul(i128::from(recipe.height_mm)).ok_or("polygon_jacobian_overflow")?;
    let det_bound=next_down(det_num as f64/36.);
    let cert=BracketCertificate{incidence,determinant_lower_bound_mm3:det_bound,radial_error_upper_mm:0.,design_faces};
    Ok(PolygonCertificate{certificate:cert,volume_sixths_mm3:volume,projected_area_twice_mm2:area2})
}

fn next_down(value:f64)->f64{
    if value.is_nan()||value==f64::NEG_INFINITY{return value}
    if value==0.0{return -f64::from_bits(1)}
    let bits=value.to_bits();
    f64::from_bits(if value>0.0{bits-1}else{bits+1})
}

fn edge_tag(a:usize,b:usize,boundary:&BTreeSet<(usize,usize)>,vertices:&[Vertex],ids:&BTreeMap<[i32;2],&str>)->Result<String,&'static str>{
    if !boundary.contains(&edge_key(a,b)){return Err("polygon_internal_face_exterior")}
    let x=ids.get(&vertices[a].xy_mm).ok_or("polygon_boundary_provenance")?;let y=ids.get(&vertices[b].xy_mm).ok_or("polygon_boundary_provenance")?;
    let (lo,hi)=if x<y{(x,y)}else{(y,x)};Ok(format!("edge/{lo}/{hi}"))
}
