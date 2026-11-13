const {test}=require('node:test');
const assert=require('node:assert/strict');
const M=require('../src/boundary-value-model.cjs');
const close=(a,b,tol=1e-9)=>assert.ok(Math.abs(a-b)<tol,`${a} != ${b}`);
test('the three endpoint targets have the claimed solution counts',()=>{
  const s=M.fresh();assert.equal(M.classification(s).kind,'unique');close(M.classification(s).C,1);
  for(const C of [-2,-1,0,.6,1,2]){
    const q={...s,C};close(M.endpointAt(q,0).y,0);close(M.endpointAt(q,0).dy,2*C);
    close(M.endpointAt({...q,caseId:'infinite'},Math.PI/2).y,0);
    close(M.endpointAt({...q,caseId:'none'},Math.PI/2).y-1,-1);
  }
  assert.equal(M.classification({...s,caseId:'infinite'}).kind,'infinite');
  assert.equal(M.classification({...s,caseId:'none'}).kind,'none');
});
test('all parameter signs solve the ODE and the left boundary, with nonzero amplitude',()=>{
  const h=1e-5;
  for(const bc of ['values','mixed'])for(const length of ['pi','two','two-pi'])for(const lambda of [-4,-1,-1e-10,0,1e-10,.25,1,2,4,9,30]){
    const s={...M.fresh(),bc,length,lambda},ell=M.L(s),a=M.candidate(s,0),b=M.candidate(s,ell);
    close(bc==='values'?a.y:a.dy,0);
    for(const x of [.05,.4,.8].map(t=>t*ell)){
      const left=M.candidate(s,x-h),p=M.candidate(s,x),right=M.candidate(s,x+h);
      close((right.y-left.y)/(2*h),p.dy,2e-7);
      close(-(right.dy-left.dy)/(2*h),lambda*p.y,2e-6);
      close(p.out,lambda*p.y);
    }
    const data=M.sample({...s,view:'eigenvalues'},1201);const peak=Math.max(...data.map(p=>Math.abs(p.y)));
    assert.ok(peak>.9999&&peak<=1+1e-12);
    if(lambda<=0)assert.ok(Math.abs(b.y)>.99);
  }
});
test('Dirichlet and mixed spectra are selected by their respective boundary conditions',()=>{
  for(const bc of ['values','mixed'])for(const length of ['pi','two','two-pi']){
    const s={...M.fresh(),bc,length,view:'eigenvalues'},ell=M.L(s);
    for(let n=1;n<=3;n++){
      const lambda=M.eigenvalue(s,n),q={...s,lambda};
      assert.equal(M.mode(q),n);close(M.residual(q),0,1e-12);
      const p=M.candidate(q,ell);close(p.y,0);
      const expected=((bc==='values'?n:(2*n-1)/2)*Math.PI/ell)**2;close(lambda,expected);
      assert.equal(M.mode({...q,lambda:lambda+.001}),null);
    }
  }
  const s={...M.fresh(),lambda:1};close(M.eigenvalue({...s,length:'two-pi'},1),M.eigenvalue(s,1)/4);
  const tiny=M.candidate({...s,lambda:0},1);close(tiny.y,1/Math.PI);
});
test('a boundary-admissible polynomial fails the eigenfunction identity on the interval',()=>{
  const h=1e-5;
  for(const bc of ['values','mixed'])for(const lambda of [0,1,4,9]){
    const s={...M.fresh(),bc,lambda,poly:true},ell=M.L(s),check=M.operatorCheck(s);
    close(check.left,0);close(check.right,0);assert.ok(check.error>0);assert.equal(check.eigenfunction,false);
    for(const x of [0,.3,.5,1].map(t=>t*ell)){
      const a=M.polynomial(s,x-h),p=M.polynomial(s,x),b=M.polynomial(s,x+h);
      close(-(b.dy-a.dy)/(2*h),p.out,1e-8);
    }
    assert.equal(M.operatorCheck({...s,poly:false}).relation,true);
  }
  // At this one point the polynomial matches, yet the identity still fails elsewhere.
  const s={...M.fresh(),lambda:8/Math.PI**2,poly:true},p=M.polynomial(s,Math.PI/2);
  close(p.out,p.scaled);assert.equal(M.operatorCheck(s).eigenfunction,false);
  assert.equal(M.operatorCheck({...M.fresh(),lambda:2}).eigenfunction,false);
  assert.equal(M.operatorCheck({...M.fresh(),lambda:4}).eigenfunction,true);
});
test('the scan selects the next exact root and saved state preserves exact mode values',()=>{
  for(const bc of ['values','mixed']){
    const s={...M.fresh(),bc,lambda:-4};let q=s;
    for(let n=1;n<=3;n++){
      const p=M.next(q);assert.equal(p.n,n);q={...q,lambda:p.lambda};assert.equal(M.mode(q),n);
    }
    assert.deepEqual(M.restore({version:1,...q}),q);
  }
  assert.deepEqual(M.restore(null),M.fresh());
  const s=M.restore({version:1,A:0,lambda:Infinity,C:'2',progress:4,bc:'other',length:'other'});
  assert.equal(s.A,.25);assert.equal(s.lambda,2);assert.equal(s.progress,1);assert.equal(s.C,.6);
});
