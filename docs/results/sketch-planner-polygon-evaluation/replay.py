#!/usr/bin/env python3
"""Replay frozen baseline model JSON through a candidate CLI's audit command."""
import argparse, hashlib, json, subprocess, time
from pathlib import Path

CORPUS=Path(__file__).with_name('cases.json')

def main():
    p=argparse.ArgumentParser(); p.add_argument('--cli',type=Path,required=True)
    p.add_argument('--baseline-models',type=Path,required=True)
    p.add_argument('--baseline-result',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True); p.add_argument('--timeout',type=float,default=60)
    a=p.parse_args(); cli=a.cli.resolve(); sha=hashlib.sha256(cli.read_bytes()).hexdigest()
    corpus=json.loads(CORPUS.read_text()); base=json.loads(a.baseline_result.read_text())
    if base['corpus_sha256']!=hashlib.sha256(CORPUS.read_bytes()).hexdigest():
        raise SystemExit('baseline result corpus SHA mismatch')
    rows=[]; start=time.monotonic()
    for case in corpus['cases']:
        model=a.baseline_models/f"{case['id']}.model.json"
        t=time.monotonic()
        try:
            proc=subprocess.run([str(cli),'audit',str(model)],text=True,capture_output=True,timeout=a.timeout)
            elapsed=time.monotonic()-t
            try: report=json.loads(proc.stdout)
            except Exception: report=None
            rows.append({'id':case['id'],'family':case['family'],'returncode':proc.returncode,
                'elapsed_s':elapsed,'status':report.get('status') if report else None,
                'method':report.get('method') if report else None,
                'determinant_lower_bound_mm3':report.get('determinant_lower_bound_mm3') if report else None,
                'volume_exact_mm3':report.get('volume_exact_mm3') if report else None,
                'stdout':None if report else proc.stdout.strip(),
                'stderr':proc.stderr.strip()})
        except subprocess.TimeoutExpired as exc:
            rows.append({'id':case['id'],'family':case['family'],'returncode':None,
                         'elapsed_s':time.monotonic()-t,'status':'timeout','stderr':str(exc)})
        print(case['id'],rows[-1].get('status'),rows[-1].get('returncode'),flush=True)
    if hashlib.sha256(cli.read_bytes()).hexdigest()!=sha: raise RuntimeError('CLI changed during replay')
    summary={}
    for fam in sorted({x['family'] for x in rows}):
        group=[x for x in rows if x['family']==fam]
        summary[fam]={'denominator':len(group),'valid':sum(x.get('returncode')==0 and x.get('status')=='valid' for x in group),
          'nonvalid':sum(x.get('returncode')!=0 or x.get('status')!='valid' for x in group)}
    output={'schema':'baseline_model_backward_replay_v1','cli_sha256':sha,
      'baseline_cli_sha256':base['cli_sha256'],'corpus_sha256':hashlib.sha256(CORPUS.read_bytes()).hexdigest(),
      'elapsed_s':time.monotonic()-start,'denominator':len(rows),'summary':summary,'cases':rows,
      'note':'This replay audits the archived baseline models. It is reported separately and does not alter fresh-recipe geometry scores.'}
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
if __name__=='__main__': main()
