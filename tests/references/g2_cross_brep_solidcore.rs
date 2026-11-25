use solidcore::curve::curve::Curve;
use solidcore::math::{Point3, Vec3};
use solidcore::query::brep::BrepClassifier;
use solidcore::intersect::point_class::PointClassType;
use solidcore::error::EntityRef;
use solidcore::topology::face::p2::{surface as face_surface, uv_domain};
use solidcore::nurbs::b_spline_surface_n_lib::p2::EvaluateFlags;
use solidcore::eval::plan::basis::Side;
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::math::CPoint;
use solidcore::sweep::api_brep::create_linear_sweep;
use solidcore::mass::integrator::properties::brep_volume_default;
use solidcore::boolean::api_brep::{boolean_intersection, boolean_union, boolean_difference};
use solidcore::Kernel;
use solidcore_ops::Body;
use std::{env,fs,process};
use solidcore_step::write::product::{solidcore_bodies,Header};

include!("input.rs");

fn profile_tool(k:&mut Kernel)->Result<solidcore::Brep,String>{
    let mut profile=Vec::<Curve>::new();
    for span in PROFILE_ARCS {
        let cp=span.map(|p| CPoint::new(p[0]*p[3],p[1]*p[3],p[2]*p[3],p[3]));
        profile.push(Curve::from(BSplineCurve::from_data(3,3,7,3,&cp,&[0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0]).map_err(|e|format!("rational Bezier input: {e:?}"))?));
    }
    let refs=profile.iter().map(Some).collect::<Vec<_>>();
    create_linear_sweep(k,&refs,Vec3::Z,FREE_LENGTH,true).map_err(|e|format!("SolidCore profile sweep: {e:?}"))
}

fn rect_tool(k:&mut Kernel)->Result<solidcore::Brep,String>{
    let rect=RECT_POINTS.windows(2).map(|p| BSplineCurve::create_line_segment(3,Point3::new(p[0][0],p[0][1],p[0][2]),Point3::new(p[1][0],p[1][1],p[1][2])).map(Curve::from).map_err(|e|format!("degree-one source line: {e:?}"))).collect::<Result<Vec<_>,_>>()?;
    let refs=rect.iter().map(Some).collect::<Vec<_>>();
    create_linear_sweep(k,&refs,Vec3::X,CROSS_LENGTH,true).map_err(|e|format!("SolidCore rectangle sweep: {e:?}"))
}

fn stock(k:&mut Kernel)->Result<solidcore::Brep,String>{
    let x0=BOUNDS[0][0];let x1=BOUNDS[0][1];let y0=BOUNDS[1][0];let y1=BOUNDS[1][1];let z0=BOUNDS[2][0];let z1=BOUNDS[2][1];
    let pts=[[x0,y0,z0],[x1,y0,z0],[x1,y1,z0],[x0,y1,z0],[x0,y0,z0]];
    let edges=pts.windows(2).map(|p| BSplineCurve::create_line_segment(3,Point3::new(p[0][0],p[0][1],p[0][2]),Point3::new(p[1][0],p[1][1],p[1][2])).map(Curve::from).map_err(|e|format!("stock edge: {e:?}"))).collect::<Result<Vec<_>,_>>()?;
    let refs=edges.iter().map(Some).collect::<Vec<_>>();
    create_linear_sweep(k,&refs,Vec3::Z,z1-z0,true).map_err(|e|format!("stock sweep: {e:?}"))
}

fn volume(k:&Kernel,b:&mut solidcore::Brep)->Result<f64,String>{brep_volume_default(k,b,1e-8).map_err(|e|format!("B-rep volume: {e:?}"))}

fn clamp(v:f64,a:f64,b:f64)->f64{v.max(a).min(b)}

fn closest_face_jet(k:&Kernel,b:&solidcore::Brep,target:[f64;3])->Option<(String,f64,[f64;3],f64,f64,f64,[f64;3],[f64;3],[f64;3],[f64;3])> {
    let target=Point3::new(target[0],target[1],target[2]);
    let mut best:Option<(String,f64,[f64;3],f64,f64,f64,[f64;3],[f64;3],[f64;3],[f64;3])>=None;
    for fid in b.face_list(){
        let Ok(Some(srf))=face_surface(b,*fid) else{continue};
        let Ok(domain)=uv_domain(b,*fid) else{continue};let ui=domain.u_interval();let vi=domain.v_interval();
        if !(ui.min().is_finite()&&ui.max().is_finite()&&vi.min().is_finite()&&vi.max().is_finite()&&ui.min()<ui.max()&&vi.min()<vi.max()){continue}
        let mut seed=(ui.min(),vi.min(),f64::INFINITY);
        for iu in 0..5{for iv in 0..5{
            let u=ui.min()+(ui.max()-ui.min())*iu as f64/4.0;let v=vi.min()+(vi.max()-vi.min())*iv as f64/4.0;
            if let Ok(p)=srf.evaluate_point(k,solidcore::math::Point2::new(u,v)){let d=(p-target).length_squared();if d<seed.2{seed=(u,v,d)}}
        }}
        let(mut u,mut v,_)=seed;
        for _ in 0..24{
            let Ok(e)=srf.evaluate(k,solidcore::math::Point2::new(u,v),2,2,EvaluateFlags::default()) else{break};
            let p=e[0][0];let du=e[1][0];let dv=e[0][1];let r=target-p;
            let(a,c,d)=(du.dot(du),du.dot(dv),dv.dot(dv));let(det)=(a*d-c*c);
            if det.abs()<1e-20{break}
            let ru=du.dot(r);let rv=dv.dot(r);
            let step_u=(d*ru-c*rv)/det;let step_v=(a*rv-c*ru)/det;
            let nu=clamp(u+step_u,ui.min(),ui.max());let nv=clamp(v+step_v,vi.min(),vi.max());
            if (nu-u).abs()+(nv-v).abs()<1e-13{u=nu;v=nv;break}u=nu;v=nv;
        }
        let uv=solidcore::math::Point2::new(u,v);
        if let Ok(e)=srf.evaluate(k,uv,2,2,EvaluateFlags::default()){
            let residual=e[0][0]-target;let dist=residual.length();let n=e[1][0].cross(e[0][1]);
            let normal=if n.length()>0.0{n/n.length()}else{Vec3::ZERO};
            if let Ok(g)=srf.evaluate_geometric(k,uv,solidcore::nurbs::knots::Side::Right,solidcore::nurbs::knots::Side::Right){
                if best.as_ref().is_none_or(|x|dist<x.1){best=Some((format!("{fid:?}"),dist,[normal.x,normal.y,normal.z],g.gaussian_curvature,g.principle_curvature1,g.principle_curvature2,[g.principle_curvature_vector1.x,g.principle_curvature_vector1.y,g.principle_curvature_vector1.z],[g.principle_curvature_vector2.x,g.principle_curvature_vector2.y,g.principle_curvature_vector2.z],[e[1][0].x,e[1][0].y,e[1][0].z],[e[0][1].x,e[0][1].y,e[0][1].z]))}
            }
        }
    }
    best
}

fn run() -> Result<(), String> {
    let mut k=Kernel::default();
    let mut p=profile_tool(&mut k)?;let vp=volume(&k,&mut p)?;
    let mut r=rect_tool(&mut k)?;let vr=volume(&k,&mut r)?;
    let pa=profile_tool(&mut k)?;let rb=rect_tool(&mut k)?;
    let mut common=boolean_intersection(&mut k,pa,rb).map_err(|e|format!("SolidCore tool intersection: {e:?}"))?;
    let vi=volume(&k,&mut common)?;
    let pa=profile_tool(&mut k)?;let rb=rect_tool(&mut k)?;
    let union=boolean_union(&mut k,pa,rb).map_err(|e|format!("SolidCore tool union: {e:?}"))?;
    let mut union=union;
    let vu=volume(&k,&mut union)?;
    let stock=stock(&mut k)?;
    let mut remaining=boolean_difference(&mut k,stock,union).map_err(|e|format!("SolidCore stock minus tool union: {e:?}"))?;
    let vrm=volume(&k,&mut remaining)?;
    let nfaces=remaining.face_list().len();let nedges=remaining.edge_list().len();let nverts=remaining.vertex_list().len();
    let closed_manifold=solidcore::query::api_queries::brep_is_closed_manifold(&remaining).map_err(|e|format!("closed-manifold check: {e:?}"))?;
    let manifold_solid=solidcore::query::api_queries::brep_is_manifold_solid(&remaining).map_err(|e|format!("manifold-solid check: {e:?}"))?;
    let mut classifier=BrepClassifier::of_brep(&remaining);
    let classes=CLASS_POINTS.iter().map(|p|{
        match classifier.point_3d_classify(&k,Point3::new(p[0],p[1],p[2]),2e-7,true){
            Ok(c)=>match (c.point_class(),c.object()){
                (PointClassType::Region,Some(EntityRef::Region(r)))=>if remaining.region(r).is_ok_and(|x|x.is_void==0){"Inside"}else{"Outside"},
                (PointClassType::Region,_)=>"RegionUnknown",
                (PointClassType::Face|PointClassType::Edge|PointClassType::Vertex,_)=>"OnBoundary",
                (PointClassType::Nothing,_)=>"Outside",
                _=>"Unclassified",
            },
            Err(_)=>"ERROR",
        }
    }).collect::<Vec<_>>();
    let jets=JET_TARGETS.iter().filter_map(|t|closest_face_jet(&k,&remaining,[t[0],t[1],t[2]]).map(|(face,d,n,gaussian,k1,k2,v1,v2,du,dv)|format!("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",t[3],t[4],face,d,n[0],n[1],n[2],gaussian,k1,k2,v1[0],v1[1],v1[2],v2[0],v2[1],v2[2],du[0],du[1],du[2],dv[0],dv[1],dv[2],t[2]))).collect::<Vec<_>>();
    let seams=SEAM_PROBES.iter().filter_map(|q|{
        let class=classifier.point_3d_classify(&k,Point3::new(q[3],q[4],q[5]),2e-7,true).ok().map(|c|match(c.point_class(),c.object()){
            (PointClassType::Region,Some(EntityRef::Region(r)))=>if remaining.region(r).is_ok_and(|x|x.is_void==0){"Inside"}else{"Outside"},
            (PointClassType::Region,_)=>"RegionUnknown",
            (PointClassType::Face|PointClassType::Edge|PointClassType::Vertex,_)=>"OnBoundary",
            (PointClassType::Nothing,_)=>"Outside",_=>"Unclassified"})?;
        let distance=closest_face_jet(&k,&remaining,[q[3],q[4],q[5]]).map(|x|x.1)?;
        Some(format!("{},{},{},{},{}",q[0],q[1],q[2],class,distance))
    }).collect::<Vec<_>>();
    let tol_probes=SEAM_TOLERANCE_PROBES.iter().filter_map(|q|{
        let class=classifier.point_3d_classify(&k,Point3::new(q[4],q[5],q[6]),2e-7,true).ok().map(|c|match(c.point_class(),c.object()){
            (PointClassType::Region,Some(EntityRef::Region(r)))=>if remaining.region(r).is_ok_and(|x|x.is_void==0){"Inside"}else{"Outside"},
            (PointClassType::Region,_)=>"RegionUnknown",
            (PointClassType::Face|PointClassType::Edge|PointClassType::Vertex,_)=>"OnBoundary",
            (PointClassType::Nothing,_)=>"Outside",_=>"Unclassified"})?;
        Some(format!("{},{},{},{},{}",q[0],q[1],q[2],q[3],class))
    }).collect::<Vec<_>>();
    let args=env::args().collect::<Vec<_>>();if args.len()!=2{return Err("usage: solidcore STEP".into())}
    let body=Body::from_kernel(remaining);
    let step=solidcore_bodies(&[("g2-cross-plate",&body)],&Header{file_name:"g2-cross-solidcore.step".into(),..Header::default()}).map_err(|e|format!("SolidCore STEP export: {e:?}"))?;
    fs::write(&args[1],step).map_err(|e|format!("write SolidCore STEP: {e}"))?;
    println!("profile_tool_volume_mm3={vp:.17}\nrectangle_tool_volume_mm3={vr:.17}\ntool_intersection_volume_mm3={vi:.17}\ntool_union_volume_mm3={vu:.17}\nremaining_solid_volume_mm3={vrm:.17}\nfaces={nfaces}\nedges={nedges}\nvertices={nverts}\nclosed_manifold={closed_manifold}\nmanifold_solid={manifold_solid}\nclassifications={}",classes.join(","));
    println!("boundary_jets={}",jets.join(";"));
    println!("seam_probes={}",seams.join(";"));
    println!("seam_tolerance_probes={}",tol_probes.join(";"));
    Ok(())
}

fn main(){if let Err(e)=run(){eprintln!("{e}");process::exit(1)}}
