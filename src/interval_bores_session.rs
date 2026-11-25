//! Transactional depth edits on a fixed polygon profile and fixed bore axes.
use crate::{
    bracket::BracketBlock,
    incidence::BlockSide,
    interval_bores::{
        build_layer, layer_intervals, prepare_profile, sort_partition, validate, IntervalBore,
        IntervalBorePrism,
    },
    interval_bores_audit::{audit_interval_bore_prism, IntervalBoreCertificate},
};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
    time::Instant,
};

type LayerKey = ([i32; 2], u8);
type Layers = BTreeMap<LayerKey, Arc<Vec<BracketBlock>>>;
type Selections = BTreeMap<String, Vec<BlockSide>>;

/// A meeting with a fixed axial plane, along the linear edit trajectory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AxialEvent {
    pub moving_cap: String,
    pub fixed_plane: String,
    pub height_mm: i32,
    /// Reduced exact fraction in [0,1], including departures from coincidences.
    pub time: [i32; 2],
}

#[derive(Debug)]
pub struct IntervalEditReport {
    pub events: Vec<AxialEvent>,
    pub old_class: &'static str,
    pub new_class: &'static str,
    pub created_features: Vec<String>,
    pub lost_features: Vec<String>,
    pub remapped_features: Vec<String>,
    pub reused_layers: usize,
    pub built_layers: usize,
    pub reused_blocks: usize,
    pub built_blocks: usize,
    /// Signed coefficient of pi in V(new) - V(old), in cubic millimetres.
    pub volume_change_pi_mm3: i64,
    pub build_s: f64,
    pub audit_s: f64,
}

/// A published state and its independently established certificate.
pub struct IntervalState {
    pub recipe: IntervalBorePrism,
    pub blocks: Vec<BracketBlock>,
    pub certificate: IntervalBoreCertificate,
    selections: Selections,
    layers: Layers,
}

pub struct IntervalBoreSession {
    profile: Vec<BracketBlock>,
    current: IntervalState,
}

fn class([lo, hi]: [i32; 2], height: i32) -> &'static str {
    match (lo == 0, hi == height) {
        (true, true) => "through",
        (true, false) => "open_start",
        (false, true) => "open_end",
        (false, false) => "enclosed",
    }
}

fn fraction(mut a: i32, mut b: i32) -> [i32; 2] {
    if b < 0 {
        a = -a;
        b = -b;
    }
    let (mut x, mut y) = (a, b);
    while y != 0 {
        (x, y) = (y, x % y);
    }
    [a / x, b / x]
}

fn crossings(recipe: &IntervalBorePrism, bore: &IntervalBore, next: [i32; 2]) -> Vec<AxialEvent> {
    let mut planes = vec![
        ("start".to_owned(), 0),
        ("end".to_owned(), recipe.height_mm),
    ];
    for other in &recipe.bores {
        if other.bore.id != bore.bore.id {
            for (k, which) in ["lower", "upper"].into_iter().enumerate() {
                planes.push((
                    format!("bore/{}/{}", other.bore.id, which),
                    other.z_interval_mm[k],
                ));
            }
        }
    }
    let mut events = Vec::new();
    for (k, which) in ["lower", "upper"].into_iter().enumerate() {
        let old = bore.z_interval_mm[k];
        let new = next[k];
        if old == new {
            continue;
        }
        for (name, z) in &planes {
            if (old.min(new)..=old.max(new)).contains(z) {
                events.push(AxialEvent {
                    moving_cap: format!("bore/{}/{}", bore.bore.id, which),
                    fixed_plane: name.clone(),
                    height_mm: *z,
                    time: fraction(*z - old, new - old),
                });
            }
        }
    }
    events.sort_by(|a, b| {
        (a.time[0] * b.time[1])
            .cmp(&(b.time[0] * a.time[1]))
            .then_with(|| a.moving_cap.cmp(&b.moving_cap))
            .then_with(|| a.fixed_plane.cmp(&b.fixed_plane))
    });
    events
}

struct Assembly {
    state: IntervalState,
    reused_layers: usize,
    reused_blocks: usize,
    build_s: f64,
    audit_s: f64,
}

fn assemble(
    profile: &[BracketBlock],
    recipe: IntervalBorePrism,
    old: &Layers,
) -> Result<Assembly, &'static str> {
    let started = Instant::now();
    let bores = validate(&recipe)?;
    let mut layers = Layers::new();
    let mut blocks = Vec::new();
    let (mut reused_layers, mut reused_blocks) = (0, 0);
    for span in layer_intervals(recipe.height_mm, &bores) {
        let mask = bores.iter().enumerate().fold(0u8, |mask, (i, b)| {
            let [lo, hi] = b.z_interval_mm;
            mask | (u8::from(lo <= span[0] && span[1] <= hi) << i)
        });
        let key = (span, mask);
        let layer = if let Some(previous) = old.get(&key) {
            reused_layers += 1;
            reused_blocks += previous.len();
            Arc::clone(previous)
        } else {
            Arc::new(build_layer(profile, &recipe, &bores, span[0], span[1])?)
        };
        blocks.extend(layer.iter().cloned());
        layers.insert(key, layer);
    }
    sort_partition(&mut blocks)?;
    let build_s = started.elapsed().as_secs_f64();
    let started = Instant::now();
    let certificate = audit_interval_bore_prism(&recipe, &blocks)?;
    let mut selections = Selections::new();
    for (side, feature) in &certificate.certificate.design_faces {
        selections
            .entry(feature.clone())
            .or_default()
            .push(side.clone());
    }
    let audit_s = started.elapsed().as_secs_f64();
    Ok(Assembly {
        state: IntervalState {
            recipe,
            blocks,
            certificate,
            selections,
            layers,
        },
        reused_layers,
        reused_blocks,
        build_s,
        audit_s,
    })
}

impl IntervalBoreSession {
    pub fn new(recipe: IntervalBorePrism) -> Result<Self, &'static str> {
        let bores = validate(&recipe)?;
        let profile = prepare_profile(&recipe, &bores)?;
        let current = assemble(&profile, recipe, &Layers::new())?.state;
        Ok(Self { profile, current })
    }

    pub fn current(&self) -> &IntervalState {
        &self.current
    }

    pub fn into_current(self) -> IntervalState {
        self.current
    }

    /// Change one excavation. Required features must exist before and after the edit.
    pub fn apply(
        &mut self,
        id: &str,
        interval: [i32; 2],
        required: &[String],
    ) -> Result<IntervalEditReport, String> {
        let mut recipe = self.current.recipe.clone();
        let index = recipe
            .bores
            .iter()
            .position(|b| b.bore.id == id)
            .ok_or("unknown interval bore")?;
        let old = recipe.bores[index].z_interval_mm;
        recipe.bores[index].z_interval_mm = interval;
        if interval[0] < 0 || interval[0] >= interval[1] || interval[1] > recipe.height_mm {
            return Err("interval_bore_z_interval".into());
        }
        for name in required {
            if !self.current.selections.contains_key(name) {
                return Err(format!("required feature absent before edit: {name}"));
            }
        }
        let events = crossings(
            &self.current.recipe,
            &self.current.recipe.bores[index],
            interval,
        );
        let height = recipe.height_mm;
        let Assembly {
            state: next,
            reused_layers,
            reused_blocks,
            build_s,
            audit_s,
        } = assemble(&self.profile, recipe, &self.current.layers)?;
        for name in required {
            if !next.selections.contains_key(name) {
                return Err(format!("required feature lost during edit: {name}"));
            }
        }
        let previous: BTreeSet<_> = self.current.selections.keys().collect();
        let current: BTreeSet<_> = next.selections.keys().collect();
        let report = IntervalEditReport {
            events,
            old_class: class(old, height),
            new_class: class(interval, height),
            created_features: current
                .difference(&previous)
                .map(|id| (**id).clone())
                .collect(),
            lost_features: previous
                .difference(&current)
                .map(|id| (**id).clone())
                .collect(),
            remapped_features: previous
                .intersection(&current)
                .filter(|id| self.current.selections[**id] != next.selections[**id])
                .map(|id| (**id).clone())
                .collect(),
            reused_layers,
            built_layers: next.layers.len() - reused_layers,
            reused_blocks,
            built_blocks: next.blocks.len() - reused_blocks,
            volume_change_pi_mm3: self.current.certificate.removed_pi_volume_coefficient_mm3
                - next.certificate.removed_pi_volume_coefficient_mm3,
            build_s,
            audit_s,
        };
        self.current = next;
        Ok(report)
    }
}
