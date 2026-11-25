from pathlib import Path
import argparse,json,subprocess,time,hashlib
root=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description="Compare the same M11 jobs on two release binaries")
for name in ("baseline","candidate","work","output"):parser.add_argument(name,type=Path)
args=parser.parse_args()
dst=args.work;dst.mkdir(exist_ok=True,parents=True)
base=args.baseline.resolve();current=args.candidate.resolve()
rows=[]
for name in ['oblique_circle_stress','extreme_shear_circle_pressure','extreme_shear_circle_vector','extreme_shear_circle_stress']:
 for version,binary in [('25676c2',base),('working',current)]:
  for budget in [65536,262144]:
   source=root/'build/adaptive-surface-load';model=source/'budget-replay'/name/'model.json'
   job=json.loads((source/name/'job.json').read_text());job['surface_quadrature']['max_panels']=budget
   jp=dst/(name+'.json');jp.write_text(json.dumps(job));out=dst/'result.json';out.write_bytes(b'keep')
   st=time.monotonic();p=subprocess.run([str(binary),'solve-job',str(model),str(jp),str(out)],capture_output=True,text=True)
   r=dict(case=name,version=version,sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
          budget=budget,code=p.returncode,diagnostic=p.stdout.strip(),elapsed_s=time.monotonic()-st)
   if p.returncode:r['preserved']=out.read_bytes()==b'keep'
   else:r['quadrature']=json.loads(out.read_text())['loads'][0]['integration']
   rows.append(r);print(name,version,budget,r.get('quadrature',{}).get('panels'),flush=True)
args.output.write_text(json.dumps(rows,indent=2)+'\n')
