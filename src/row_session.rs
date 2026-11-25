//! Transactional whole-face row updates with reusable local Q2 stiffness matrices.
use crate::{
    analysis_space::{from_certified_faces, AnalysisBlock, AnalysisSpace},
    bracket::BracketBlock,
    bracket_audit::BracketCertificate,
    elasticity::{assemble, gauss, Matrix, Stiffness},
    row_bores::build_row_bores,
    row_bores_audit::audit_row,
};
use std::{collections::BTreeMap, sync::Arc, time::Instant};

#[derive(Clone, Debug, PartialEq)]
pub struct RowRequest {
    pub half_widths: [i32; 2],
    pub centers: Vec<[i32; 2]>,
    pub radii: Vec<i32>,
    pub young_mpa: f64,
    pub poisson: f64,
    pub gauss_order: usize,
}
#[derive(Debug, Default, Clone, Copy)]
pub struct RowUpdate {
    pub reused: usize,
    pub integrated: usize,
    pub build_s: f64,
    pub audit_numbering_s: f64,
    pub element_s: f64,
    pub merge_s: f64,
}
struct Element {
    controls: Vec<[i128; 4]>,
    stiffness: Stiffness,
}
pub struct RowState {
    pub request: RowRequest,
    pub blocks: Vec<BracketBlock>,
    pub certificate: BracketCertificate,
    pub space: AnalysisSpace,
    pub stiffness: Stiffness,
    pub update: RowUpdate,
    contributions: BTreeMap<String, Arc<Element>>,
}
#[derive(Default)]
pub struct RowAnalysisSession {
    state: Option<RowState>,
}

impl RowAnalysisSession {
    pub fn state(&self) -> Option<&RowState> {
        self.state.as_ref()
    }

    /// Rebuild topology and numbering, re-integrating only changed local data.
    /// All candidates are owned locally until the final successful assignment.
    pub fn update(&mut self, request: RowRequest) -> Result<&RowState, &'static str> {
        if !request.young_mpa.is_finite()
            || request.young_mpa <= 0.
            || !(-1. < request.poisson && request.poisson < 0.5)
        {
            return Err("row_session_invalid_material");
        }
        gauss(request.gauss_order)?;
        let mut timing = RowUpdate::default();
        let t = Instant::now();
        let blocks = build_row_bores(
            request.half_widths[0],
            request.half_widths[1],
            &request.centers,
            &request.radii,
        )
        .map_err(|_| "row_session_geometry_refused")?
        .blocks;
        timing.build_s = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let certificate = audit_row(
            request.half_widths,
            &request.centers,
            &request.radii,
            &blocks,
        )?;
        let space =
            from_certified_faces(&blocks, &certificate.incidence, &certificate.design_faces)?;
        timing.audit_numbering_s = t.elapsed().as_secs_f64();
        let old = self.state.as_ref().filter(|s| {
            s.request.young_mpa == request.young_mpa
                && s.request.poisson == request.poisson
                && s.request.gauss_order == request.gauss_order
        });
        let mut contributions = BTreeMap::new();
        let t = Instant::now();
        for block in &space.blocks {
            let cached = old
                .and_then(|s| s.contributions.get(&block.id))
                .filter(|e| e.controls == block.controls);
            let element = if let Some(e) = cached {
                timing.reused += 1;
                Arc::clone(e)
            } else {
                let local = AnalysisSpace {
                    ndof: 27,
                    blocks: vec![AnalysisBlock {
                        id: block.id.clone(),
                        controls: block.controls.clone(),
                        dofs: (0..27).collect(),
                    }],
                    boundary_dofs: BTreeMap::new(),
                };
                let stiffness = assemble(
                    &local,
                    request.young_mpa,
                    request.poisson,
                    request.gauss_order,
                )?;
                timing.integrated += 1;
                Arc::new(Element {
                    controls: block.controls.clone(),
                    stiffness,
                })
            };
            contributions.insert(block.id.clone(), element);
        }
        timing.element_s = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let mut rows = vec![BTreeMap::<usize, Matrix>::new(); space.ndof];
        let mut integrated_volume = 0.;
        for block in &space.blocks {
            let element = &contributions[&block.id].stiffness;
            integrated_volume += element.integrated_volume;
            for (local, &global) in block.dofs.iter().enumerate() {
                for k in element.offsets[local]..element.offsets[local + 1] {
                    let col = block.dofs[element.columns[k]];
                    let entry = rows[global].entry(col).or_default();
                    for (target, source) in entry.iter_mut().zip(&element.values[k]) {
                        for (value, addend) in target.iter_mut().zip(source) {
                            *value += addend;
                        }
                    }
                }
            }
        }
        let mut stiffness = Stiffness {
            offsets: vec![0],
            columns: Vec::new(),
            values: Vec::new(),
            integrated_volume,
        };
        for row in rows {
            for (col, v) in row {
                stiffness.columns.push(col);
                stiffness.values.push(v);
            }
            stiffness.offsets.push(stiffness.columns.len());
        }
        stiffness.validate()?;
        if !integrated_volume.is_finite() || integrated_volume <= 0. {
            return Err("row_session_volume");
        }
        timing.merge_s = t.elapsed().as_secs_f64();
        self.state = Some(RowState {
            request,
            blocks,
            certificate,
            space,
            stiffness,
            update: timing,
            contributions,
        });
        Ok(self.state.as_ref().unwrap())
    }
}
