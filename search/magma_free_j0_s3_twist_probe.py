#!/usr/bin/env python3
"""Exact finite-field rank probes for the S3-equivariant j=0 twist family.

For
  q=t^2-t+1, r=t(t-1), a=q(mu), b=q(mu)^3/r(mu),
  E_mu: y^2=x^3+b^2r^2-a^3q^3,
we probe the representative packet twist by -t:
  Y^2=X^3-t^3(b^2r^2-a^3q^3).

The other two normalized packet channels are Q(t)-isomorphic by S3 symmetry
with constants (-1,1,-1).  Magma's finite-field function-field analytic
information is discovery evidence; characteristic-zero promotion still
requires explicit sections and an independent certificate.
"""
from __future__ import annotations
import argparse,json,time,traceback,urllib.parse,urllib.request,xml.etree.ElementTree as ET
from pathlib import Path
SERVER='http://magma.maths.usyd.edu.au/xml/calculator.xml';REFERER='http://magma.maths.usyd.edu.au/calc/'

def code(p,mu):
 return f'''SetColumns(0);
F:=GF({p}); K<t>:=FunctionField(F);
mu:=F!({mu});
q:=t^2-t+1; r:=t*(t-1);
a:=mu^2-mu+1; rmu:=mu*(mu-1);
assert rmu ne 0;
b:=a^3/rmu;
B:=K!(b^2)*r^2-K!(a^3)*q^3;
assert B ne 0;
T:=EllipticCurve([K|0,0,0,0,-t^3*B]);
assert Discriminant(T) ne 0;
AI:=AnalyticInformation(T);
print "R30_PRIME",{p};
print "R30_MU",{mu};
print "R30_ANALYTIC_INFORMATION",AI;
print "R30_ARITHMETIC_RANK",AI[1];
print "R30_GEOMETRIC_RANK",AI[2];
print "R30_J0_S3_TWIST_PASS";
'''
def submit(name,program,outdir):
 st=time.time();data=urllib.parse.urlencode({'input':program}).encode();req=urllib.request.Request(SERVER,data=data,headers={'Content-Type':'application/x-www-form-urlencoded','Accept':'application/xml,text/xml,text/html','Referer':REFERER,'User-Agent':'elliptic-rank-30-j0-probe/1.0'},method='POST');raw=b''
 try:
  with urllib.request.urlopen(req,timeout=180) as resp:raw=resp.read();http=resp.status
  (outdir/f'{name}.xml').write_bytes(raw);root=ET.fromstring(raw);warn=' '.join(''.join(w.itertext()).strip() for w in root.findall('.//warning')).strip();lines=[''.join(x.itertext()) for x in root.findall('.//results/line')];text='\n'.join(lines);(outdir/f'{name}.txt').write_text(text+'\n')
  rank=grank=None
  for line in text.splitlines():
   z=line.split()
   if z[:1]==['R30_ARITHMETIC_RANK'] and len(z)>1:rank=int(z[1])
   if z[:1]==['R30_GEOMETRIC_RANK'] and len(z)>1:grank=int(z[1])
  passed='R30_J0_S3_TWIST_PASS' in text and not warn
  return {'status':'pass' if passed else 'fail','rank':rank,'geometric_rank':grank,'warning':warn,'http_status':http,'elapsed_seconds':time.time()-st,'output_tail':text[-3000:]}
 except Exception:return {'status':'exception','elapsed_seconds':time.time()-st,'traceback':traceback.format_exc()}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--prime',type=int,required=True);ap.add_argument('--mus',default='2,3,4,5');ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--summary',type=Path,required=True);a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True);rows=[]
 for mu in map(int,a.mus.split(',')):
  if mu%a.prime in (0,1):continue
  name=f'j0_s3_p{a.prime}_mu{mu}';program=code(a.prime,mu);(a.output_dir/f'{name}.m').write_text(program);rec=submit(name,program,a.output_dir);rec.update({'prime':a.prime,'mu':mu});rows.append(rec)
 out={'status':'pass','prime':a.prime,'runs':rows,'best':max(rows,key=lambda r:r.get('rank',-1) if r.get('rank') is not None else -1),'truth_note':'finite-characteristic exact rank data only; no Q(t) lower-bound claim'};a.summary.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'prime':a.prime,'runs':[(r['mu'],r['status'],r.get('rank'),r.get('geometric_rank')) for r in rows]}));return 0
if __name__=='__main__':raise SystemExit(main())
