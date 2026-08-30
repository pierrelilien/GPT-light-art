'use strict';
const C={
 W:20,H:20,N:6,nAir:1.0003,nFilm:1.75,nSlab:1.49,filmNm:185,th:.20,
 internalBounces:3,sceneInteractions:1,minEnergy:.0007,absorptionPerCm:.015,wallZ:-1.5,
 refBeamFull:28,beamProfileExp:1.25,cam:{t:{x:10,y:10,z:1.5},d:48,z:1.16},
 pieces:[{x:6.5,y:13.5,s:4,a:25},{x:12,y:10.5,s:3,a:105},{x:10.5,y:5.5,s:2.5,a:-20}],
 wavelengths:[630,540,450]
};
const S={view:3,debug:0,auto:0,ang:-180,rad:27,z:12,beam:17,yaw:-.75,pitch:.72,drag:0,px:0,py:0};
const c=document.querySelector('#c'),m=document.querySelector('#m'),g=c.getContext('2d'),h=m.getContext('2d'),st=document.querySelector('#st');
const V=(x,y,z)=>({x,y,z}),add=(a,b)=>V(a.x+b.x,a.y+b.y,a.z+b.z),sub=(a,b)=>V(a.x-b.x,a.y-b.y,a.z-b.z),mul=(a,k)=>V(a.x*k,a.y*k,a.z*k),dot=(a,b)=>a.x*b.x+a.y*b.y+a.z*b.z,cross=(a,b)=>V(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x),len=a=>Math.hypot(a.x,a.y,a.z),norm=a=>mul(a,1/(len(a)||1)),cl=(x,a,b)=>Math.max(a,Math.min(b,x));
const specMul=(a,b)=>[a[0]*b[0],a[1]*b[1],a[2]*b[2]],specScale=(a,k)=>[a[0]*k,a[1]*k,a[2]*k],specEnergy=a=>.2126*a[0]+.7152*a[1]+.0722*a[2],specColor=a=>{let mx=Math.max(...a,1e-9),q=a.map(v=>Math.pow(cl(v/mx,0,1),.72));let mn=Math.min(...q),sat=.10;q=q.map(v=>cl((v-mn)*(1+sat)+mn,0,1));return `rgb(${Math.round(q[0]*255)} ${Math.round(q[1]*255)} ${Math.round(q[2]*255)})`};
const src=()=>{let a=S.ang*Math.PI/180;return V(10+S.rad*Math.cos(a),10+S.rad*Math.sin(a),S.z)},aim=()=>V(10,10,0);
function frame(p){let a=p.a*Math.PI/180,u=V(Math.cos(a),Math.sin(a),0);return{u,n:norm(V(-Math.sin(a),Math.cos(a),0)),b:V(p.x,p.y,0)}}
function point(p,u,z){let f=frame(p);return add(add(f.b,mul(f.u,(u-.5)*p.s)),V(0,0,z*p.s))}
function corners(p){let f=frame(p),a=add(f.b,mul(f.u,-p.s/2)),b=add(f.b,mul(f.u,p.s/2));return[a,b,add(b,V(0,0,p.s)),add(a,V(0,0,p.s))]}
function insideFace(q,p,eps=.002){let f=frame(p),r=sub(q,f.b),u=dot(r,f.u);return Math.abs(u)<=p.s/2+eps&&q.z>=-eps&&q.z<=p.s+eps}
function planeT(o,d,p,n){let den=dot(d,n);if(Math.abs(den)<1e-10)return null;let t=dot(sub(p,o),n)/den;return t>1e-5?t:null}
function at(o,d,t){return add(o,mul(d,t))}
function reflect(i,n){return norm(sub(i,mul(n,2*dot(i,n))))}
function refract(i,n,n1,n2){let ci=cl(-dot(n,i),0,1),eta=n1/n2,k=1-eta*eta*(1-ci*ci);return k<0?null:norm(add(mul(i,eta),mul(n,eta*ci-Math.sqrt(k))))}
function thinFilmRgb(nIn,nOut,ci){
 ci=cl(ci,0,1);let si=Math.sqrt(Math.max(0,1-ci*ci)),sf=nIn/C.nFilm*si,so=nIn/nOut*si;
 if(so>=1-1e-9)return[1,1,1];
 let cf=Math.sqrt(Math.max(0,1-sf*sf)),co=Math.sqrt(Math.max(0,1-so*so)),out=[];
 for(let lam of C.wavelengths){
  let r01s=(nIn*ci-C.nFilm*cf)/(nIn*ci+C.nFilm*cf),r12s=(C.nFilm*cf-nOut*co)/(C.nFilm*cf+nOut*co);
  let r01p=(C.nFilm*ci-nIn*cf)/(C.nFilm*ci+nIn*cf),r12p=(nOut*cf-C.nFilm*co)/(nOut*cf+C.nFilm*co);
  let ph=4*Math.PI*C.nFilm*C.filmNm*cf/lam,cp=Math.cos(ph);
  let Rs=(r01s*r01s+r12s*r12s+2*r01s*r12s*cp)/(1+(r01s*r12s)**2+2*r01s*r12s*cp);
  let Rp=(r01p*r01p+r12p*r12p+2*r01p*r12p*cp)/(1+(r01p*r12p)**2+2*r01p*r12p*cp);
  out.push(cl(.5*(Rs+Rp),0,1));
 }
 return out;
}
function beamHalfRad(){return S.beam*Math.PI/360}
function solidAngle(fullDeg){let a=fullDeg*Math.PI/360;return 2*Math.PI*(1-Math.cos(a))}
function powerGain(){return solidAngle(C.refBeamFull)/solidAngle(S.beam)}
function distanceGain(q){let d0=Math.hypot(27,12),d=len(sub(q,src()));return(d0/Math.max(d,1))**2}
function beamIntensity(q){let a=norm(sub(aim(),src())),d=norm(sub(q,src())),co=Math.cos(beamHalfRad()),mu=dot(a,d);if(mu<=co)return 0;let profile=((mu-co)/(1-co))**C.beamProfileExp;return profile*powerGain()*distanceGain(q)}
function absorption(dist){return Math.exp(-C.absorptionPerCm*dist)}
function boardCandidate(o,d){let t=planeT(o,d,V(0,0,0),V(0,0,1));if(t==null)return null;let p=at(o,d,t);return p.x>=0&&p.x<=C.W&&p.y>=0&&p.y<=C.H?{t,p,type:'board'}:null}
function wallCandidate(o,d){let t=planeT(o,d,V(0,0,C.wallZ),V(0,0,1));return t==null?null:{t,p:at(o,d,t),type:'wall'}}
function pieceCandidate(o,d,idx){let p=C.pieces[idx],f=frame(p),side=dot(sub(o,f.b),f.n)>=0?1:-1,n=mul(f.n,side),face=add(f.b,mul(n,C.th/2)),t=planeT(o,d,face,n);if(t==null)return null;let q=at(o,d,t);return insideFace(q,p)?{t,p:q,type:'piece',idx,n}:null}
function sceneHit(o,d,exclude=-1,allowPieces=true){let best=null;for(let q of [boardCandidate(o,d),wallCandidate(o,d)])if(q&&(!best||q.t<best.t))best=q;if(allowPieces)for(let i=0;i<C.pieces.length;i++)if(i!==exclude){let q=pieceCandidate(o,d,i);if(q&&(!best||q.t<best.t))best=q}return best}
function emitReceiver(events,o,d,spec,path,depth,exclude){
 if(specEnergy(spec)<C.minEnergy)return;
 let hit=sceneHit(o,d,exclude,depth<=C.sceneInteractions);
 if(!hit)return;
 if(hit.type==='piece'&&depth<=C.sceneInteractions){interactPiece(events,hit.idx,hit.p,d,hit.n,spec,path+`>P${hit.idx}`,depth+1);return}
 if(hit.type==='piece'){hit=sceneHit(o,d,exclude,false);if(!hit)return}
 events.push({p:hit.p,surface:hit.type,spec,path:path+':'+hit.type,energy:specEnergy(spec),origin:o});
}
function interactPiece(events,idx,entry,inc,n,spec,path,depth){
 let p=C.pieces[idx],ci=cl(-dot(n,inc),0,1);if(ci<=1e-5||!insideFace(entry,p))return;
 let Rrgb=thinFilmRgb(C.nAir,C.nSlab,ci),Trgb=Rrgb.map(x=>1-x);
 let sR=specMul(spec,Rrgb),sIn=specMul(spec,Trgb);
 emitReceiver(events,entry,reflect(inc,n),sR,path+':R0',depth,idx);
 let inside=refract(inc,n,C.nAir,C.nSlab);if(!inside||specEnergy(sIn)<C.minEnergy)return;
 let f=frame(p),opp=add(f.b,mul(n,-C.th/2)),backT=planeT(entry,inside,opp,n);if(backT==null)return;
 let back=at(entry,inside,backT);if(!insideFace(back,p))return;
 sIn=specScale(sIn,absorption(len(sub(back,entry))));let dir=inside,cur=back;
 for(let order=0;order<=C.internalBounces;order++){
  let cib=cl(-dot(n,dir),0,1);let Rb=thinFilmRgb(C.nSlab,C.nAir,cib),Tb=Rb.map(x=>1-x),outB=refract(dir,n,C.nSlab,C.nAir);
  let sOut=specMul(sIn,Tb);if(outB)emitReceiver(events,cur,outB,sOut,path+`:T${order}`,depth,idx);
  sIn=specMul(sIn,Rb);if(specEnergy(sIn)<C.minEnergy)break;
  dir=reflect(dir,n);let frontPlane=add(f.b,mul(n,C.th/2)),tf=planeT(cur,dir,frontPlane,mul(n,-1));if(tf==null)break;let front=at(cur,dir,tf);if(!insideFace(front,p))break;
  sIn=specScale(sIn,absorption(len(sub(front,cur))));let nf=mul(n,-1),cif=cl(-dot(nf,dir),0,1),Rf=thinFilmRgb(C.nSlab,C.nAir,cif),Tf=Rf.map(x=>1-x),outF=refract(dir,nf,C.nSlab,C.nAir),sF=specMul(sIn,Tf);
  if(outF)emitReceiver(events,front,outF,sF,path+`:R${order+1}`,depth,idx);
  sIn=specMul(sIn,Rf);if(specEnergy(sIn)<C.minEnergy)break;
  dir=reflect(dir,nf);let tb=planeT(front,dir,opp,n);if(tb==null)break;cur=at(front,dir,tb);if(!insideFace(cur,p))break;
  sIn=specScale(sIn,absorption(len(sub(cur,front))));
 }
}
function tracePiece(idx){
 let p=C.pieces[idx],A=[],s=src(),f=frame(p);
 for(let i=0;i<C.N;i++)for(let j=0;j<C.N;j++){
  let q=point(p,(i+.5)/C.N,(j+.5)/C.N),inc=norm(sub(q,s)),side=dot(sub(s,f.b),f.n)>=0?1:-1,n=mul(f.n,side),face=add(f.b,mul(n,C.th/2)),t=planeT(s,inc,face,n);if(t==null)continue;
  let entry=at(s,inc,t);if(!insideFace(entry,p))continue;let blocker=sceneHit(s,inc,idx,true);if(blocker&&blocker.type==='piece'&&blocker.t<t-1e-4)continue;let ci=cl(-dot(n,inc),0,1),B=beamIntensity(entry),geom=ci,E=B*geom;if(E<=C.minEnergy)continue;
  let rec={events:[],incident:E};interactPiece(rec.events,idx,entry,inc,n,[E,E,E],`P${idx}`,0);A.push(rec)
 }
 return A;
}
let cache={key:'',sets:[]};function traces(){let k=`${S.ang}|${S.rad}|${S.z}|${S.beam}`;if(k!==cache.key){cache.key=k;cache.sets=C.pieces.map((_,i)=>tracePiece(i))}return cache.sets}
function hull(P){P=P.filter(Boolean).slice().sort((a,b)=>a.x-b.x||a.y-b.y);if(P.length<3)return P;let X=(o,a,b)=>(a.x-o.x)*(b.y-o.y)-(a.y-o.y)*(b.x-o.x),L=[],U=[];for(let p of P){while(L.length>1&&X(L.at(-2),L.at(-1),p)<=0)L.pop();L.push(p)}for(let i=P.length-1;i>=0;i--){let p=P[i];while(U.length>1&&X(U.at(-2),U.at(-1),p)<=0)U.pop();U.push(p)}L.pop();U.pop();return L.concat(U)}
function groups(records){let M=new Map;for(let r of records)for(let e of r.events){let a=M.get(e.path)||[];a.push(e);M.set(e.path,a)}let out=[];for(let [key,a] of M){if(a.length<2)continue;let s=[0,0,0];for(let e of a){s[0]+=e.spec[0];s[1]+=e.spec[1];s[2]+=e.spec[2]}s=s.map(v=>v/a.length);out.push({key,P:hull(a.map(e=>e.p)),spec:s,color:specColor(s),power:a.reduce((z,e)=>z+e.energy,0)/a.length,rep:a[Math.floor(a.length/2)]})}return out.sort((a,b)=>b.power-a.power)}