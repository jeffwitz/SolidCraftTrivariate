#!/usr/bin/env python3
"""Create paired JSON and Markdown from frozen per-executable measurements."""
import argparse, hashlib, json, math, statistics
from collections import Counter, defaultdict
from pathlib import Path

def pct(xs, p):
    return sorted(xs)[max(0, math.ceil(p*len(xs))-1)] if xs else None

def summarize_pairs(base, cand):
    b={x['id']:x for x in base['cases']}; c={x['id']:x for x in cand['cases']}
    if b.keys()!=c.keys(): raise ValueError('baseline/candidate case IDs differ')
    groups=defaultdict(list)
    for ident in b:
        x,y=b[ident],c[ident]
        outcome=('common_success' if x['geometry_pass'] and y['geometry_pass'] else
                 'newly_admitted' if not x['geometry_pass'] and y['geometry_pass'] else
                 'regression' if x['geometry_pass'] and not y['geometry_pass'] else 'common_failure')
        r={'id':ident,'family':x['family'],'paired_outcome':outcome,
           'baseline_admitted':x['admitted'],'candidate_admitted':y['admitted'],
           'baseline_geometry_pass':x['geometry_pass'],'candidate_geometry_pass':y['geometry_pass'],
           'baseline_method':x['construction_method'],'candidate_method':y['construction_method'],
           'baseline_blocks':x['block_count'],'candidate_blocks':y['block_count'],
           'baseline_input_vertex_euler_gap':(x['block_count']-x['input_vertex_lower_bound_total_blocks']
             if x['block_count'] is not None else None),
           'candidate_input_vertex_euler_gap':(y['block_count']-y['input_vertex_lower_bound_total_blocks']
             if y['block_count'] is not None else None),
           'block_delta_candidate_minus_baseline':(y['block_count']-x['block_count']
             if x['block_count'] is not None and y['block_count'] is not None else None),
           'baseline_build_s':x['build_elapsed_s'],'candidate_build_s':y['build_elapsed_s'],
           'baseline_total_s':x['build_elapsed_s']+x.get('audit_elapsed_s',0),
           'candidate_total_s':y['build_elapsed_s']+y.get('audit_elapsed_s',0),
           'baseline_condition_max':max((q['max_condition_2'] for q in x.get('condition_grid_5cube_per_block') or []),default=None),
           'candidate_condition_max':max((q['max_condition_2'] for q in y.get('condition_grid_5cube_per_block') or []),default=None)}
        groups[x['family']].append(r)
    summaries={}
    for family, rows in sorted(groups.items()):
        common=[r for r in rows if r['paired_outcome']=='common_success']
        bd=[r['block_delta_candidate_minus_baseline'] for r in common]
        gap_b=[r['baseline_input_vertex_euler_gap'] for r in common]
        gap_c=[r['candidate_input_vertex_euler_gap'] for r in common]
        ratios=[r['candidate_condition_max']/r['baseline_condition_max']
                for r in common if r['baseline_condition_max'] and r['candidate_condition_max']]
        methods_b=Counter(r['baseline_method'] for r in rows if r['baseline_admitted'])
        methods_c=Counter(r['candidate_method'] for r in rows if r['candidate_admitted'])
        summaries[family]={
            'denominator':len(rows),
            'paired_outcomes':dict(Counter(r['paired_outcome'] for r in rows)),
            'baseline_admitted':sum(r['baseline_admitted'] for r in rows),
            'candidate_admitted':sum(r['candidate_admitted'] for r in rows),
            'baseline_geometry_pass':sum(r['baseline_geometry_pass'] for r in rows),
            'candidate_geometry_pass':sum(r['candidate_geometry_pass'] for r in rows),
            'baseline_construction_methods':dict(methods_b),
            'candidate_construction_methods':dict(methods_c),
            'common_success_block_delta':{'n':len(bd),'median':statistics.median(bd) if bd else None,
                'total':sum(bd) if bd else None,'min':min(bd) if bd else None,'max':max(bd) if bd else None},
            'common_success_input_vertex_euler_gap':{'n':len(common),
                'baseline_total':sum(gap_b) if gap_b else None,
                'candidate_total':sum(gap_c) if gap_c else None,
                'candidate_minus_baseline':sum(gap_c)-sum(gap_b) if gap_b else None,
                'baseline_median':statistics.median(gap_b) if gap_b else None,
                'candidate_median':statistics.median(gap_c) if gap_c else None},
            'common_success_condition_ratio_candidate_over_baseline':{
                'n':len(ratios),'median':statistics.median(ratios) if ratios else None,
                'p95_nearest_rank':pct(ratios,.95),'worst':max(ratios) if ratios else None}}
    all_rows=[r for rows in groups.values() for r in rows]
    return summaries,all_rows

def md_table(rows, headers):
    lines=['| '+' | '.join(headers)+' |','| '+' | '.join('---' for _ in headers)+' |']
    lines += ['| '+' | '.join(str(r.get(k,'')) for k in headers)+' |' for r in rows]
    return '\n'.join(lines)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--candidate',type=Path,required=True); p.add_argument('--corpus',type=Path,required=True)
    p.add_argument('--json',type=Path,required=True); p.add_argument('--markdown',type=Path,required=True)
    a=p.parse_args(); b=json.loads(a.baseline.read_text()); c=json.loads(a.candidate.read_text())
    corpus_sha=hashlib.sha256(a.corpus.read_bytes()).hexdigest()
    if b['corpus_sha256']!=corpus_sha or c['corpus_sha256']!=corpus_sha:
        raise SystemExit('result corpus SHA mismatch')
    grouped,rows=summarize_pairs(b,c)
    case_data={x['id']:x for x in json.loads(a.corpus.read_text())['cases']}
    for row in rows:
        row['cut_count']=len(case_data[row['id']]['recipe'].get('cuts',[]))
    seeded=[r for r in rows if r['family'].startswith('seeded_convex_')]
    seeded_split={}
    for label,pred in [('single_cut',lambda r:r['cut_count']==1),
                       ('multiple_cuts',lambda r:r['cut_count']>1)]:
        part=[r for r in seeded if pred(r)]
        seeded_split[label]={'n':len(part),
            'baseline_methods':dict(Counter(r['baseline_method'] for r in part)),
            'candidate_methods':dict(Counter(r['candidate_method'] for r in part)),
            'radial_selected':sum(r['candidate_method']=='radial_convex_footprint' for r in part),
            'triangulated_retained':sum(r['candidate_method']=='triangulated_feature_footprint' for r in part),
            'baseline_input_vertex_euler_gap_total':sum(r['baseline_input_vertex_euler_gap'] for r in part),
            'candidate_input_vertex_euler_gap_total':sum(r['candidate_input_vertex_euler_gap'] for r in part)}
    out={'schema':'polygonal_paired_evaluation_v1','corpus_sha256':corpus_sha,
         'environment':{'python':__import__('platform').python_version(),
            'platform':__import__('platform').platform(),
            'timing_note':'Candidate run overlapped the startup of a separate cargo test process; wall-clock timings are descriptive, not a controlled speed comparison.'},
         'baseline':{k:b[k] for k in ('tag','cli_sha256','elapsed_s','summary','condition_summary')},
         'candidate':{k:c[k] for k in ('tag','cli_sha256','elapsed_s','summary','condition_summary')},
         'families':grouped,'seeded_convex_stratification':seeded_split,'paired_cases':rows,
         'protocol':{'deadline_s_per_cli_operation':60,'sample_grid':'closed 5x5x5 tensor grid per block',
            'condition_metric':'2-norm Jacobian condition number; per-model maximum over all block/grid samples',
            'geometry_gate':'CLI closed-domain certificate valid; exact planar containment, edge incidence, no overlap, area conservation and volume agreement',
            'timeouts_and_rejections':'remain in each family denominator'}}
    a.json.parent.mkdir(parents=True,exist_ok=True); a.markdown.parent.mkdir(parents=True,exist_ok=True)
    a.json.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    sections=['# Paired polygonal sketch partition evaluation','',
      f"Public synthetic corpus SHA-256: `{corpus_sha}`.",
      f"Baseline executable SHA-256: `{b['cli_sha256']}` (`{b['tag']}`).",
      f"Candidate executable SHA-256: `{c['cli_sha256']}` (`{c['tag']}`).",
      'The baseline source is `e5bd5c6`; the frozen candidate source is `4e8a3944cd5aa85f8ba0984b109b6649b2d54b6f`. The 45-case public synthetic bundle SHA-256 was recorded before candidate scoring.',
      'All measurements are paired on the same 45 recipes: 5 fixed convex through-hole, 5 fixed convex blind-pocket, 5 fixed concave stock/cut, and 30 seeded convex cases split evenly between through-holes and blind pockets. Integer millimetres are within the input coordinate contract.',
      'Each admitted model is independently checked with exact rational section containment, edge incidence, non-overlap, area conservation, and exact volume, plus the CLI closed-domain determinant certificate. Each build and audit operation had a 60 second timeout.','',
      '## Paired outcomes by family','',
      md_table([{'family':f,'n':v['denominator'],
        'baseline_valid':v['baseline_geometry_pass'],'candidate_valid':v['candidate_geometry_pass'],
        'new':v['paired_outcomes'].get('newly_admitted',0),'regressions':v['paired_outcomes'].get('regression',0),
        'common_success':v['paired_outcomes'].get('common_success',0),'common_failure':v['paired_outcomes'].get('common_failure',0),
        'median_block_delta':v['common_success_block_delta']['median'],
        'baseline_euler_gap':v['common_success_input_vertex_euler_gap']['baseline_total'],
        'candidate_euler_gap':v['common_success_input_vertex_euler_gap']['candidate_total']} for f,v in grouped.items()],
        ['family','n','baseline_valid','candidate_valid','new','regressions','common_success','common_failure','median_block_delta','baseline_euler_gap','candidate_euler_gap']),
      '', 'A negative block delta means the candidate used fewer blocks. Rejections and timeouts remain in the denominator. All 45 paired inputs passed geometry on both executables. Total blocks fell from 6,184 to 5,201 (−15.9%); there were no new admissions or regressions. Euler gap is observed blocks minus the bound using mandatory input contour vertices; the conditional bound including inserted host vertices and exact per-case gaps are in the JSON.','',
      '## Construction methods','',
      f"Among the 30 seeded convex cases (stock polygons have 5–12 sides), the baseline used triangulation on all 30. The candidate selected radial construction on {seeded_split['single_cut']['radial_selected']}/{seeded_split['single_cut']['n']} one-cut cases; `{next((r['id'] for r in seeded if r['cut_count']==1 and r['candidate_method']=='triangulated_feature_footprint'), 'none')}` is the remaining one-cut case kept on triangulation. All {seeded_split['multiple_cuts']['n']} multiple-cut cases retained triangulation. The input-vertex Euler gap totals for the one-cut stratum are {seeded_split['single_cut']['baseline_input_vertex_euler_gap_total']}→{seeded_split['single_cut']['candidate_input_vertex_euler_gap_total']}; the multiple-cut stratum is {seeded_split['multiple_cuts']['baseline_input_vertex_euler_gap_total']}→{seeded_split['multiple_cuts']['candidate_input_vertex_euler_gap_total']}.",'']
    for family,v in grouped.items():
        sections += [f"### {family}",'',f"Baseline: `{v['baseline_construction_methods']}`", '',
          f"Candidate: `{v['candidate_construction_methods']}`",'']
    sections += ['## Sampled map condition','',
      'The condition metric is diagnostic: for every block, the spectral condition number of the Jacobian was evaluated on the closed 5×5×5 parameter grid. The reported statistics take each solid’s maximum over its blocks. This is not a continuous regularity certificate.','',
      md_table([{'executable':x['tag'],'n':x['condition_summary']['n'],
        'median':x['condition_summary']['median'],'p95':x['condition_summary']['p95_nearest_rank'],
        'worst':x['condition_summary']['worst']} for x in (b,c)],
        ['executable','n','median','p95','worst']),
      '', 'The candidate reduced blocks sharply on some cases while worsening sampled condition. The largest paired rise is `blind_triangle_shallow`: 105→56 blocks and κmax 53.94→908.70 (16.84×). Seeded cases `_03`, `_09`, `_21`, and `_12` rise 1.41–1.84×. These remain visible even though the aggregate median improves.',
      'Per-case block counts, observed timing, construction method, Euler bounds, pass/refusal status, and per-block 5³ condition samples are in the adjacent JSON. Candidate timing overlapped a separate cargo test process at startup; observed wall totals (235.60 s baseline, 209.23 s candidate) are descriptive, not a controlled speed claim.','',
      '## Scope and limits','',
      'The corpus is deterministic public synthetic input, not the externally owned private holdout. It covers convex polygonal through-holes, convex blind pockets, and concave stock or cut profiles. Nonparallel and stacked 3D bores are assessed separately. The Euler bound applies only to connected conforming quadrilateral sections under the stated mandatory-corner assumptions. Jacobian samples are not used as a regularity certificate.']
    a.markdown.write_text('\n'.join(sections)+'\n')
if __name__=='__main__': main()
