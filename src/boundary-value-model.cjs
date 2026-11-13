// Exact scalar boundary-value examples and two homogeneous eigenvalue problems.
(function(root){
  const cases=Object.freeze({
    unique:{b:Math.PI/4,target:1,label:'y(π/4) = 1'},
    infinite:{b:Math.PI/2,target:0,label:'y(π/2) = 0'},
    none:{b:Math.PI/2,target:1,label:'y(π/2) = 1'}
  });
  const lengths=Object.freeze({pi:Math.PI,two:2,'two-pi':2*Math.PI});
  const fresh=()=>({view:'endpoints',caseId:'unique',C:.6,length:'pi',bc:'values',lambda:2,A:1,poly:false,family:true,reveal:false,progress:1});
  const L=s=>lengths[s.length];
  function endpointAt(s,x){return {x,y:s.C*Math.sin(2*x),dy:2*s.C*Math.cos(2*x),out:4*s.C*Math.sin(2*x)};}
  function classification(s){
    const {b,target}=cases[s.caseId],factor=Math.sin(2*b);
    if(Math.abs(factor)<1e-12)return {kind:target===0?'infinite':'none',C:null,factor:0};
    return {kind:'unique',C:target/factor,factor};
  }
  function candidate(s,x,lambda=s.lambda){
    const ell=L(s),A=s.A;
    let y,dy;
    if(lambda===0){y=s.bc==='values'?A*x/ell:A;dy=s.bc==='values'?A/ell:0;}
    else if(lambda>0){
      const q=Math.sqrt(lambda);
      if(s.bc==='values'){
        const den=q*ell<Math.PI/2?Math.sin(q*ell):1;
        y=A*Math.sin(q*x)/den;dy=A*q*Math.cos(q*x)/den;
      }else{y=A*Math.cos(q*x);dy=-A*q*Math.sin(q*x);}
    }else{
      const q=Math.sqrt(-lambda);
      if(s.bc==='values'){
        const den=Math.sinh(q*ell);
        y=A*Math.sinh(q*x)/den;dy=A*q*Math.cosh(q*x)/den;
      }else{
        const den=Math.cosh(q*ell);
        y=A*Math.cosh(q*x)/den;dy=A*q*Math.sinh(q*x)/den;
      }
    }
    return {x,y,dy,out:lambda*y,scaled:lambda*y};
  }
  function polynomial(s,x){
    const ell=L(s),A=s.A;
    const y=s.bc==='values'?4*A*x*(ell-x)/(ell*ell):A*(1-x*x/(ell*ell));
    const dy=s.bc==='values'?4*A*(ell-2*x)/(ell*ell):-2*A*x/(ell*ell);
    const out=(s.bc==='values'?8:2)*A/(ell*ell);
    return {x,y,dy,out,scaled:s.lambda*y};
  }
  const eigenvalue=(s,n)=>((s.bc==='values'?n:n-.5)*Math.PI/L(s))**2;
  function mode(s){
    if(s.lambda<=0)return null;
    const n=Math.round(Math.sqrt(s.lambda)*L(s)/Math.PI+(s.bc==='mixed'?.5:0));
    return n>=1&&Math.abs(s.lambda-eigenvalue(s,n))<=1e-10*Math.max(1,s.lambda)?n:null;
  }
  function spectrum(s,max=30){
    const result=[];for(let n=1;eigenvalue(s,n)<=max+1e-10;n++)result.push({n,lambda:eigenvalue(s,n)});return result;
  }
  const next=s=>spectrum(s).find(p=>p.lambda>s.lambda+1e-9)||null;
  // Scaling every candidate to max |y| = 1 preserves the zeros of this residual.
  const residual=(s,lambda=s.lambda)=>candidate({...s,A:1},L(s),lambda).y;
  function sample(s,count=601){
    const end=s.view==='endpoints'?cases[s.caseId].b:L(s);
    const at=s.view==='endpoints'?endpointAt:s.view==='operator'&&s.poly?polynomial:candidate;
    return Array.from({length:count},(_,i)=>at(s,end*i/(count-1)));
  }
  function operatorCheck(s){
    const at=s.poly?polynomial:candidate,left=at(s,0),right=at(s,L(s));
    let error=0;for(let i=0;i<=80;i++){const p=at(s,L(s)*i/80);error=Math.max(error,Math.abs(p.out-p.scaled));}
    return {left:s.bc==='values'?left.y:left.dy,right:right.y,error,relation:!s.poly,eigenfunction:!s.poly&&mode(s)!==null};
  }
  function restore(raw){
    const s=fresh();if(raw?.version!==1)return s;
    if(['endpoints','eigenvalues','operator'].includes(raw.view))s.view=raw.view;
    if(Object.hasOwn(cases,raw.caseId))s.caseId=raw.caseId;
    if(Object.hasOwn(lengths,raw.length))s.length=raw.length;
    if(['values','mixed'].includes(raw.bc))s.bc=raw.bc;
    for(const [key,lo,hi] of [['C',-2,2],['lambda',-4,30],['A',.25,2],['progress',0,1]])
      if(typeof raw[key]==='number'&&Number.isFinite(raw[key]))s[key]=Math.max(lo,Math.min(hi,raw[key]));
    for(const key of ['poly','family','reveal'])if(typeof raw[key]==='boolean')s[key]=raw[key];
    return s;
  }
  const M=Object.freeze({fresh,cases,lengths,L,endpointAt,classification,candidate,polynomial,eigenvalue,mode,spectrum,next,residual,sample,operatorCheck,restore});
  if(typeof module==='object'&&module.exports)module.exports=M;else root.BoundaryValueModel=M;
})(globalThis);
