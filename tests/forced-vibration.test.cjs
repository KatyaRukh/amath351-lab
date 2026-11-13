const {test}=require('node:test');
const assert=require('node:assert/strict');
const M=require('../src/forced-vibration-model.cjs');
const close=(a,b,tol=1e-9)=>assert.ok(Math.abs(a-b)<tol,`${a} != ${b}`);

test('the six lecture presets match independently derived solutions',()=>{
  for(const t of [0,.2,1,3,7,12]){
    close(M.at(M.preset('bounded'),t).x,(Math.cos(t)-Math.cos(2*t))/3);
    close(M.at(M.preset('beats'),t).x,(Math.cos(1.8*t)-Math.cos(2*t))/(4-1.8**2));
    close(M.at(M.preset('resonance'),t).x,t*Math.sin(2*t)/4);
    close(M.at(M.preset('damped'),t).x,Math.sin(2*t)/4-Math.exp(-t)*Math.sin(Math.sqrt(3)*t)/(2*Math.sqrt(3)));
    close(M.at(M.preset('practice-free'),t).x,t*Math.sin(3*t));
    const b=Math.sqrt(8);
    close(M.at(M.preset('practice-damped'),t).x,Math.sin(3*t)-3/b*Math.exp(-t)*Math.sin(b*t));
  }
});
test('the full solution satisfies the ODE, initial values, and decomposition',()=>{
  const h=1e-5;
  for(const m of [.5,1,3])for(const c of [0,.05,1,4,6])for(const omega of [.1,1,2,3,4]){
    const s={...M.fresh(),m,c,omega,x0:-.7,v0:1.3};
    close(M.at(s,0).x,s.x0);close(M.at(s,0).v,s.v0);
    for(const t of [.2,2,7]){
      const a=M.at(s,t-h),p=M.at(s,t),b=M.at(s,t+h);
      close(p.x,p.xc+p.xp,1e-10);close(p.v,p.vc+p.vp,1e-10);
      close((b.x-a.x)/(2*h),p.v,1e-6);
      close(m*(b.v-a.v)/(2*h)+c*p.v+s.k*p.x,p.force,3e-6);
      close(m*(b.vp-a.vp)/(2*h)+c*p.vp+s.k*p.xp,p.force,3e-6);
    }
  }
});
test('near-resonance evaluation is stable and the envelope contains zero-data responses',()=>{
  const s=M.preset('resonance');
  for(const t of [0,.01,1,24,65,120])for(const eps of [-1e-9,0,1e-9]){
    const p=M.at({...s,omega:2+eps},t);
    assert.ok(Number.isFinite(p.x));close(p.x,M.at(s,t).x,4e-6);
  }
  for(const id of ['bounded','beats','resonance','practice-free']){
    const s=M.preset(id);for(const p of M.sample(s))assert.ok(Math.abs(p.x)<=M.envelopeAt(s,p.t)+1e-10);
  }
  assert.equal(M.envelopeAt({...s,x0:1},2),null);
  close(M.envelopeAt(s,120),30);
});
test('positive damping removes unbounded resonance and the transient decays',()=>{
  for(const id of ['damped','practice-damped']){
    const s=M.preset(id),g=M.coefficients(s);
    assert.ok(Number.isFinite(g.R));close(M.at(s,100).xc,0,1e-12);
    for(const t of [80,81,90])assert.ok(Math.abs(M.at(s,t).xp)<=g.R+1e-10);
  }
  assert.equal(M.response(M.preset('resonance'),2),Infinity);
  close(M.response(M.preset('damped'),2),.25);
  close(M.response(M.preset('practice-damped'),3),1);
  for(const c of [.05,.5,1,2,4,6]){
    const s={...M.fresh(),c},p=M.peak(s);
    for(const w of [0,.2,1,2,3,4])assert.ok(M.response(s,w)<=p.R+1e-10);
    if(c<Math.sqrt(8))assert.ok(p.omega<2&&p.omega>0);else assert.equal(p.omega,0);
  }
});
test('saved values are bounded and quantized to controls; each preset can be restored',()=>{
  assert.deepEqual(M.restore(null),M.fresh());
  assert.deepEqual(M.restore({version:2,c:2}),M.fresh());
  for(const id of Object.keys(M.presets))assert.deepEqual(M.restore({version:1,...M.preset(id)}),M.preset(id));
  const s=M.restore({version:1,m:10,k:100,c:1e-30,omega:500,F:-2,x0:Infinity,v0:'4',time:999,duration:120,compare:true});
  assert.equal(s.m,1);assert.equal(s.k,4);assert.equal(s.c,0);assert.equal(s.omega,4);assert.equal(s.F,0);assert.equal(s.time,120);assert.equal(s.compare,true);
  for(const p of M.sample(s))close(p.x,0);
});
