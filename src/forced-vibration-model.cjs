// Analytic responses of m*x'' + c*x' + k*x = F*cos(omega*t).
(function(root){
  const Free=typeof module==='object'&&module.exports?require('./vibration-model.cjs'):root.VibrationModel;
  const defaults={m:1,k:4,c:0,F:1,omega:1,x0:0,v0:0,time:0,duration:18,speed:1,parts:false,envelope:true,compare:false};
  const presets=Object.freeze({
    bounded:{title:'A bounded forced response',duration:18,omega:1},
    beats:{title:'Nearby frequencies produce beats',duration:65,speed:4,omega:1.8},
    resonance:{title:'Resonance: a growing envelope',duration:24,speed:2,omega:2},
    damped:{title:'Transient and steady-state response',duration:12,omega:2,c:2,parts:true},
    'practice-free':{title:'Practice: no damping',duration:12,k:9,F:6,omega:3,compare:true},
    'practice-damped':{title:'Practice: with damping',duration:12,k:9,F:6,omega:3,c:2,compare:true}
  });
  const fresh=()=>({...defaults});
  const preset=id=>{const {title,...values}=presets[id]||presets.bounded;return {...fresh(),...values};};
  const sinc=u=>Math.abs(u)<1e-4?1-u*u/6+u**4/120:Math.sin(u)/u;
  const resonant=s=>s.c===0&&Math.abs(s.omega-Math.sqrt(s.k/s.m))<1e-12;
  function coefficients(s){
    const omega0=Math.sqrt(s.k/s.m),D=s.k-s.m*s.omega*s.omega,H=s.c*s.omega,den=D*D+H*H;
    if(resonant(s))return {omega0,resonant:true,A:0,B:0,R:s.F===0?0:Infinity};
    const A=s.F*D/den,B=s.F*H/den;
    return {omega0,resonant:false,A,B,R:s.F/Math.hypot(D,H)};
  }
  function at(s,t){
    const g=coefficients(s),w=s.omega,w0=g.omega0;
    let x,v,xp,vp,xc,vc;
    if(s.c===0){
      // Product form avoids cancellation as omega approaches omega0.
      const sum=w0+w,diff=w0-w,half=sum*t/2,q=diff*t/2;
      const driven=s.F/s.m*t/sum*Math.sin(half)*sinc(q);
      const drivenV=s.F/s.m/sum*(Math.sin(half)*Math.cos(q)+t*sum/2*Math.cos(half)*sinc(q));
      const initial=Free.at(s,t);x=initial.x+driven;v=initial.v+drivenV;
      xp=g.resonant?driven:g.A*Math.cos(w*t);
      vp=g.resonant?drivenV:-w*g.A*Math.sin(w*t);
      xc=x-xp;vc=v-vp;
    }else{
      xp=g.A*Math.cos(w*t)+g.B*Math.sin(w*t);
      vp=w*(-g.A*Math.sin(w*t)+g.B*Math.cos(w*t));
      const hom=Free.at({...s,x0:s.x0-g.A,v0:s.v0-w*g.B},t);
      xc=hom.x;vc=hom.v;x=xc+xp;v=vc+vp;
    }
    if(t===0){x=s.x0;v=s.v0;}
    const force=s.F*Math.cos(w*t);
    return {t,x,v,xc,vc,xp,vp,force,a:(force-s.c*v-s.k*x)/s.m};
  }
  const envelopeAvailable=s=>s.c===0&&s.x0===0&&s.v0===0&&s.F>0;
  function envelopeAt(s,t){
    if(!envelopeAvailable(s))return null;
    const w0=Math.sqrt(s.k/s.m);
    return s.F/s.m*t/(w0+s.omega)*Math.abs(sinc((w0-s.omega)*t/2));
  }
  function sample(s,count=Math.max(721,Math.ceil(s.duration*40)+1)){
    return Array.from({length:count},(_,i)=>at(s,s.duration*i/(count-1)));
  }
  function response(s,w,c=s.c){
    const den=Math.hypot(s.k-s.m*w*w,c*w);
    return s.F===0?0:den<1e-12?Infinity:s.F/den;
  }
  function peak(s,c=s.c){
    if(c===0)return {omega:Math.sqrt(s.k/s.m),R:s.F===0?0:Infinity};
    const omega=Math.sqrt(Math.max(0,s.k/s.m-c*c/(2*s.m*s.m)));
    return {omega,R:response(s,omega,c)};
  }
  function matchingPreset(s){
    return Object.keys(presets).find(id=>{const p=preset(id);return ['m','k','c','F','omega','x0','v0'].every(k=>Math.abs(s[k]-p[k])<1e-10);})||null;
  }
  function restore(raw){
    const s=fresh();if(raw?.version!==1)return s;
    // The UI keeps mass fixed; stiffness is selected by the lecture/practice preset.
    if(raw.k===4||raw.k===9)s.k=raw.k;
    const limits={c:[0,6,.05],F:[0,6,.1],omega:[.1,4,.01],x0:[-2,2,.1],v0:[-4,4,.1]};
    for(const [key,[lo,hi,step]] of Object.entries(limits))if(typeof raw[key]==='number'&&Number.isFinite(raw[key]))s[key]=Number((Math.round(Math.max(lo,Math.min(hi,raw[key]))/step)*step).toFixed(6));
    if([12,18,24,65,120].includes(raw.duration))s.duration=raw.duration;
    if([1,2,4].includes(raw.speed))s.speed=raw.speed;
    if(typeof raw.time==='number'&&Number.isFinite(raw.time))s.time=Math.max(0,Math.min(s.duration,raw.time));
    for(const key of ['parts','envelope','compare'])if(typeof raw[key]==='boolean')s[key]=raw[key];
    return s;
  }
  const M=Object.freeze({fresh,preset,presets,coefficients,resonant,at,sample,envelopeAvailable,envelopeAt,response,peak,matchingPreset,restore});
  if(typeof module==='object'&&module.exports)module.exports=M;else root.ForcedVibrationModel=M;
})(globalThis);
