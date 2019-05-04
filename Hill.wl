(* ::Package:: *)

BeginPackage["Hill`"]
eccentricAnomaly::usage="eccentricAnomaly[ma, e]--Converts mean anomaly (ma) to eccentric anomaly for an orbit with eccentricity e";
trueAnomaly::usage="trueAnomaly[ma, e]--Converts mean anomaly (ma) to true anomaly for an orbit with eccentricity e";
meanAnomaly;
f0;
t0;
\[Omega];
DD;
v0;
tt;
f1;
ftidal;
ttidal::usage="ttidal[a]: Time at which the binary crosses the tidal sphere. a=Log[rp/rt]";
initConditions;
initEccentric;
sol::usage="sol[a, \[Phi], S, EVEC->{0,0,0}]: Solve parabolic Hill's equations for binary separation (in units 
of (m/M)^1/3 rp) as a function of time (in units of rp^1.5/(GM)^0.5))--See Sari et al 2010. a=Log[rp/rt], \[Phi]
is the phase of the binary (mean anomaly in the case of an eccentric binary), and S is the orientation of the binary with respect to the 
center of mass (1 is prograde -1 is retrograde). EVEC is the eccentricity vector of the binary.";
pos1::usage="pos1[a, \[Phi], S, EVEC->{0,0,0}, MR->10^6, Q->1]: Position of star 1--Plug binary separation from 
sol into the equation of motion for star 1. a=Log[rp/rt], \[Phi] is the phase of the binary (mean anomaly in the case of an eccentric binary), and S is the orientation of the binary with respect to the 
center of mass (1 is prograde -1 is retrograde). EVEC is the eccentricity vector of the binary. MR is the ratio between the mass of the central 
object and the binary. Q=m1/m2 (where m1 (m2) is the mass of star 1 (star 2).";
pos2::usage="pos2[a, \[Phi], S, EVEC->{0,0,0}, MR->10^6, Q->1]: Position of star 2 (see documentation for pos1.";
posBound::usage="posBound[a, \[Phi], S, EVEC->{0,0,0}, MR->10^6, Q->1]. Returns the position of the bound star.";
posUnbound;
en1;
enf::usage="enf[a, \[Phi], S, EVEC->{0,0,0}]: Energy of the star 1 after encounter with central object in units of G m1 m2/abin (M/m)^(1/3))(Negative energies mean the star is bound).  a=Log[rp/rt], \[Phi]
is the phase of the binary (mean anomaly in the case of an eccentric binary), and S is the orientation of the binary with respect to the 
center of mass (1 is prograde -1 is retrograde). EVEC is the eccentricity vector of the binary. Return \[ImaginaryI] if binary is not separated, and -\[ImaginaryI] if the 
the integration did not finish...";
enfColl::usage="enfColl[a, \[Phi], S, EVEC->{0,0,0}]: See enf. Also returns the minimum separation of the binary";
getEcc::usage="getEcc[a, \[Phi], S, EVEC->0, MR->10^6, Q->1]: Get orbital elements of the bound star";
enf2::usage="enf2[a, \[Phi], S, EVEC->{0,0,0}, MR->10^6, Q->1]";
(*getEcc2*)
ppPlot;
ppPlot2;
pComp;
pComp2;
DisruptFracGrid::usage="DisruptFracGrid[n, S, EVEC->{0,0,0}], calculates dimensionless energy and minimum separation 
for a range of binary phases and penetration factors. n is the number of grid points to use in \[Phi]. (The grid in pericenter
is hard-coded). S is the orientation of the binary with respect to the 
center of mass (1 is prograde -1 is retrograde). EVEC is the eccentricity vector of the binary.";
Begin["`Private`"]
eccentricAnomaly[ma_,e_,\[Epsilon]_:10^-16]:= Module[{xOld=ma+1000,x=ma,err=0,count=0, countMax=100, ini=1,x1,x2},
While[Abs[x-xOld]>\[Epsilon]&& count<countMax,
xOld=x;
If[count<ini,
x1=ma+e*Sin[x]; x2=ma+e*Sin[x1];x=(x1+x2)*0.5,
x+=(ma-x+e*Sin[x])/(1.-e*Cos[x])
];
count++
];
err=x-ma-e*Sin[x];
x
(*{x,count,err}*)]
(*AnomalyH[ma_,e_,\[Epsilon]_:10^-16]:= Module[{xOld=ma+1000,x=ma,err=0,count=0, countMax=100, ini=1,x1,x2},
While[Abs[x-xOld]>\[Epsilon]&& count<countMax,
xOld=x;
If[count<ini,
x1=e*Sinh[x]-ma; x2=e*Sinh[x1]-ma;x=(x1+x2)*0.5,
x+=(ma+x-e*Sinh[x])/(e*Cosh[x]-1)
];
count++
];
err=e*Sinh[x]-x-ma;
x
(*{x,count,err}*)]*)

trueAnomaly[ma_, e_,\[Epsilon]_:10^-16]:=Module[{eccAnomaly, ss, cc},
eccAnomaly=eccentricAnomaly[ma, e, \[Epsilon]];
ss=Sqrt[1-e^2]Sin[eccAnomaly]/(1-e Cos[eccAnomaly]);
cc=(Cos[eccAnomaly]-e)/(1-e Cos[eccAnomaly]);
ArcTan[cc, ss]
]
(*trueAnomalyH[ma_, e_,\[Epsilon]_:10^-16]:=Module[{eccAnomaly, ss, cc},
eccAnomaly=AnomalyH[ma, e, \[Epsilon]];
2 ArcTan[((e+1)/(e-1))^0.5Tanh[eccAnomaly/2]]
]*)

meanAnomaly[nu_, e_]:=Module[{eccAnomaly, ss, cc},
ss=Sqrt[1-e^2]Sin[nu]/(1+e Cos[nu]);
cc=(e+Cos[nu])/(1+e Cos[nu]);
eccAnomaly=ArcTan[cc, ss];
eccAnomaly-e Sin[eccAnomaly]
]

(*DD=0.1;*)
QDEF=1;
EVECDEF={0,0,0};
MRDEF=10^6;
DMINDEF=0;

DD[a_]:=10^a;
\[Omega][a_]:=DD[a]^(3/2);
v0[a_]:=DD[a]^(1/2);
f0[a_]:=-ArcCos[1/5*DD[a]-1];
t0[a_]:=Sqrt[2]/3 Tan[f0[a]/2](3+Tan[f0[a]/2]^2);
f1[t_,a_]:=(x/.FindRoot[t==Sqrt[2]/3 Tan[x/2](3+Tan[x/2]^2), {x,0}])
tt[f_]:=Sqrt[2]/3 Tan[f/2](3+Tan[f/2]^2);
ftidal[a_]:=-ArcCos[2 DD[a]-1];

ttidal[a_]:=Sqrt[2]/3 Tan[ftidal[a]/2](3+Tan[ftidal[a]/2]^2);

(*Initial phase for a circular binary.*)
\[Phi]0[a_,\[Phi]_, S_]:=\[Phi]+S \[Omega][a] t0[a];
(*Initial conditions for separation vector of an eccentric binary.*)
Options[initEccentric]={EVEC->EVECDEF};
initEccentric[a_, \[Phi]_, S_, opts :OptionsPattern[]]:=Module[{evec, ecc, w, \[Phi]1, rr, per,sma,vx, vy, jj, test},
evec=OptionValue[EVEC];
(*Period (code units) and eccentricity of the binary.*)
per=2 \[Pi]/\[Omega][a];
ecc=Norm[evec];
(*Initial true anomaly*)
\[Phi]1= trueAnomaly[\[Phi]+2 \[Pi] S (t0[a])/per, ecc];
(*Initial pomega.*)
w=If[evec[[1]]==0 && evec[[2]]==0, 0,ArcTan[evec[[1]], evec[[2]]]];
(*SMA of binary*)
sma=1/DD[a];
rr=(sma(1-ecc^2))/(1+ecc Cos[\[Phi]1]);
(*Initial velocity --this is the time derivative of r2-r1*)
{vx, vy}={-(( evec[[2]]+Sin[\[Phi]1+w])/(S Sqrt[(1-ecc^2)]))v0[a], ( evec[[1]]+Cos[\[Phi]1+w])/(S Sqrt[(1-ecc^2)]) v0[a]};
jj=Sqrt[ (1-ecc^2)];
test=Cross[{vx, vy, 0}, {0, 0, jj}]/v0[a]-{Cos[\[Phi]1+w], Sin[\[Phi]1+w], 0};
{rr*Cos[\[Phi]1+w],rr Sin[\[Phi]1+w],vx, vy}
]

(*Parabolic Hill's equations for a binary with a non-zero initial eccentricity; The eccentricity vector of the binary can be specificied 
with the EVEC optional argument.*)
Options[eqnsEccentric]={EVEC->EVECDEF};
eqnsEccentric[a_, \[Phi]_, S_, opts :OptionsPattern[]]:=Module[{initCond},
initCond=initEccentric[a,\[Phi], S, opts];
{D[x[t],t,t]==(1+Cos[f[t]])^3/8 (-x[t]+3(x[t] Cos[f[t]]+y[t] Sin[f[t]])Cos[f[t]])-x[t]/(x[t]^2+y[t]^2)^(3/2),
D[y[t],t,t]==(1+Cos[f[t]])^3/8 (-y[t]+3(x[t] Cos[f[t]]+y[t] Sin[f[t]])Sin[f[t]])-y[t]/(x[t]^2+y[t]^2)^(3/2), 
f'[t]==Sqrt[2 ](1+Cos[f[t]])^2/4,
x'[t0[a]]==initCond[[3]], y'[t0[a]]==initCond[[4]], f[t0[a]]==f0[a], x[t0[a]]==initCond[[1]], y[t0[a]]==initCond[[2]]}
];
Options[sol]={EVEC->EVECDEF};
sol[aa1_, \[Phi]1_, S_,opts :OptionsPattern[]]:=
Module[{ dmin},
(*dmin=OptionValue[DMIN];*)
(*Solving parabolic Hill's equations. Detect extrema in distance--useful for detecting the minimum binary separation.*)
NDSolve[{eqnsEccentric[aa1, \[Phi]1, S, opts], WhenEvent[(x'[t]x[t]+y'[t] y[t])Sqrt[x[t]^2+y[t]^2]==0, Sow[t]]}, {x, x', y, y',  f}, {t,t0[aa1], -xend t0[aa1]}, Method->{"StiffnessSwitching", Method->{"ExplicitRungeKutta", Automatic},  Method->{"SymbolicProcessing"->0}}, MaxSteps->STEPS, PrecisionGoal->PGOAL, AccuracyGoal->AGOAL][[1]]
];
Options[initConditions]={Q->QDEF, EVEC->EVECDEF, MR->MRDEF};
initConditions[aa1_,\[Phi]1_, S_, opts:OptionsPattern[]]:=Module[{ii, mr, q, r1, r2, ff0, xm1, ym1, ff01},
mr=OptionValue[MR];
q=OptionValue[Q];
ii=initEccentric[aa1, \[Phi]1, S, FilterRules[{opts}, Options[initEccentric]]];
ff0=f0[aa1];
ff01=Sqrt[2](1+Cos[ff0])^2/4;
xm1=2mr^(1/3) ((- Sin[ff0])/(1+Cos[ff0])+(Cos[ff0]Sin[ff0]/(1+Cos[ff0])^2)) ff01;
ym1=2mr^(1/3) ((Cos[ff0])/(1+Cos[ff0])+(Sin[ff0]^2/(1+Cos[ff0])^2)) ff01;
r1={2 mr^(1/3) Cos[ff0]/(1+Cos[ff0])-q/(1+q) ii[[1]], 2 mr^(1/3) Sin[ff0]/(1+Cos[ff0])-q/(1+q) ii[[2]],0, xm1-q/(1+q) ii[[3]],ym1-q/(1+q) ii[[4]], 0};
r2={2 mr^(1/3) Cos[ff0]/(1+Cos[ff0])+1/(1+q) ii[[1]], 2 mr^(1/3) Sin[ff0]/(1+Cos[ff0])+1/(1+q) ii[[2]],0, xm1+1/(1+q) ii[[3]],ym1+1/(1+q) ii[[4]], 0};
{r1, r2}
]
Options[pos1]={Q->QDEF, EVEC->EVECDEF, MR->MRDEF};
pos1[aa1_, \[Phi]1_, S_,opts :OptionsPattern[]]:=Module[{mr, q, ss, xx, yy, ff0, ff01, xx1, yy1, xm1, ym1, ii, tend},
ss=sol[aa1, \[Phi]1, S, FilterRules[{opts}, Options[sol]]];
tend=((x/.ss)[[1,1,2]]);
mr=OptionValue[MR];
q=OptionValue[Q];
{xx, yy, xx1, yy1}={x, y, x', y'}/.ss;
ff0=f0[aa1];
ii=initConditions[aa1, \[Phi]1, S, opts][[1]];

NDSolve[{D[r1x[t], t,t]==-mr/(r1x[t]^2+r1y[t]^2)^(3/2) r1x[t]+(q/(1+q))/(xx[t]^2+yy[t]^2)^(3/2) xx[t],D[r1y[t], t,t]==-mr/(r1x[t]^2+r1y[t]^2)^(3/2) r1y[t]+(q/(1+q))/(xx[t]^2+yy[t]^2)^(3/2) yy[t], 
r1x[t0[aa1]]==ii[[1]], r1y[t0[aa1]]==ii[[2]], r1x'[t0[aa1]]==ii[[4]],  r1y'[t0[aa1]]==ii[[5]]}, 
{r1x, r1y, r1x', r1y'}, {t, t0[aa1], tend}]
]
Options[pos2]={Q->QDEF, EVEC->EVECDEF, MR->MRDEF};
pos2[aa1_, \[Phi]1_, S_,opts :OptionsPattern[]]:=Module[{mr, q, ss, xx, yy, ff0, ff01, xx1, yy1, xm1, ym1,ii, tend},
(*Print[FilterRules[{opts}, Options[sol]]];*)
ss=sol[aa1, \[Phi]1, S, FilterRules[{opts}, Options[sol]]];
tend=((x/.ss)[[1,1,2]]);
mr=OptionValue[MR];
q=OptionValue[Q];
{xx, yy, xx1, yy1}={x, y, x', y'}/.ss;
ff0=f0[ aa1];
(*Print[ff0//N];*)
ii=initConditions[aa1, \[Phi]1, S, opts][[2]];
NDSolve[{D[r1x[t], t,t]==-mr/(r1x[t]^2+r1y[t]^2)^(3/2) r1x[t]-(1/(1+q))/(xx[t]^2+yy[t]^2)^(3/2) xx[t],
D[r1y[t], t,t]==-mr/(r1x[t]^2+r1y[t]^2)^(3/2) r1y[t]-(1/(1+q))/(xx[t]^2+yy[t]^2)^(3/2) yy[t], 
r1x[t0[aa1]]==ii[[1]], r1y[t0[aa1]]==ii[[2]],
r1x'[t0[aa1]]==ii[[4]],  r1y'[t0[aa1]]==ii[[5]]}, {r1x, r1y, r1x', r1y'}, {t, t0[aa1], tend}]
]

Options[posBound]={MR->MRDEF, Q->QDEF, EVEC->EVECDEF};
posBound[aa1_, \[Phi]1_, S_, opts :OptionsPattern[]]:=Module[{q, mr, ss, xx, yy, ff0, ff01, xx1, yy1, xm1, ym1, ef,bound,pp},
(*Print[FilterRules[opts, Options[sol]]];*)
ef =enf[aa1, \[Phi]1, S, FilterRules[{opts}, Options[enf]]];
pp=If[ef<0, pos1[aa1, \[Phi]1, S, opts], pos2[aa1, \[Phi]1, S, opts]];
pp
]
Options[posUnbound]={MR->MRDEF, Q->QDEF, EVEC->EVECDEF};
posUnbound[aa1_, \[Phi]1_, S_, opts :OptionsPattern[]]:=Module[{q, mr, ss, xx, yy, ff0, ff01, xx1, yy1, xm1, ym1, ef,bound,pp},
(*Print[FilterRules[opts, Options[sol]]];*)
ef =enf[aa1, \[Phi]1, S, FilterRules[{opts}, Options[enf]]];
pp=If[ef<0, pos2[aa1, \[Phi]1, S, opts], pos1[aa1, \[Phi]1, S, opts]];
pp
]

Options[en1]={EVEC->EVECDEF};
en1[t_, ss_, aa1_?NumericQ, \[Phi]1_?NumericQ, S_, opts :OptionsPattern[]]:=Module[{}, 
(-1/DD[aa1] ((1+Cos[f[t]])^2/4 (x[t] Cos[f[t]]+y[t] Sin[f[t]])+(-Sin[f[t]]x'[t]+(1+Cos[f[t]])y'[t])/Sqrt[2])/.ss)
]

Options[enf]={EVEC->EVECDEF}
enf[aa1_?NumericQ,\[Phi]1_?NumericQ, S_, opts:OptionsPattern[]]:=
Module[{ss, filt, tt, tend, e1, e2,minSep, times, tmp},
ss=sol[aa1, \[Phi]1, S, opts];
tt=-xend t0[aa1];
tend=((x/.ss)[[1,1]][[2]]);
If[tend<tt, Return[-I, Module]];

e1=en1[tend, ss, aa1, \[Phi]1, S, opts];
tend=0.9 tend;
e2=en1[tend, ss, aa1, \[Phi]1, S, opts];
tend=tend/0.9;
filt=(((Sqrt[x[tt]^2+y[tt]^2]>10/DD[aa1]/.ss)) &&(Abs[(e1-e2)/e2]<1.2));

Piecewise[{{I, Not[filt]}}, en1[tt, ss, aa1, \[Phi]1, S, opts]]
]

Options[enfColl]={EVEC->EVECDEF}
enfColl[aa1_?NumericQ,\[Phi]1_?NumericQ, S_, opts:OptionsPattern[]]:=
Module[{ss, filt, tt, tend, e1, e2,minSep, times, tmp, opts2},
opts2=FilterRules[{opts}, Options[sol]];
tmp=Reap[sol[aa1, \[Phi]1, S, opts2]];
ss=tmp[[1]];
times=tmp[[-1]][[1]];
minSep=Min[(x[#]^2+y[#]^2/.ss)^0.5 &/@times];
tt=-xend t0[aa1];
tend=((x/.ss)[[1,1]][[2]]);
If[tend<tt, Return[{-I, -I}, Module]];

e1=en1[tend, ss, aa1, \[Phi]1, S, opts2];
tend=0.9 tend;
e2=en1[tend, ss, aa1, \[Phi]1, S, opts2];
tend=tend/0.9;
filt=(((Sqrt[x[tt]^2+y[tt]^2]>10/DD[aa1]/.ss)) &&(Abs[(e1-e2)/e2]<1.2));
Piecewise[{{{I, I}, Not[filt]}}, {en1[tt, ss, aa1, \[Phi]1, S, opts], minSep DD[aa1]}]

]
Options[enf2]={EVEC->EVECDEF, MR->MRDEF, Q->QDEF};
enf2[aa1_, \[Phi]_, S_, opts :OptionsPattern[]]:=Module[{mr, sol1, sol2, solp1, solp2, teval, p1}, 
mr=OptionValue[MR];
teval=-xend t0[aa1];
p1=pos1[aa1, \[Phi] , S, opts][[1]];
sol1=r1x/.p1;
sol2=r1y/.p1;
solp1=r1x'/.p1;
solp2=r1y'/.p1;
Print[{sol1[t0[aa1]], sol2[t0[aa1]], solp1[t0[aa1]] , solp2[t0[aa1]]}];
1/DD[aa1] (-mr/(sol1[teval]^2+sol2[teval]^2)^0.5+1/2 (solp1[teval]^2+solp2[teval]^2))(2 mr^(-1/3))
]

(*Getting eccentricity vector...*)
Options[getEcc]={DMIN->DMINDEF, EVEC->EVECDEF, MR->MRDEF, Q->QDEF};
getEcc[aa1_, \[Phi]1_, S_, opts:OptionsPattern[]]:=Module[{teval, mr,ss,  pos, vel, kk, jc, rhat,eVec, nu, ecc,pomega, q, as, Sq, ef, abin, minSep},
{kk, minSep}=enfColl[aa1, \[Phi]1, S, FilterRules[{opts}, Options[enfColl]]];

mr=OptionValue[MR];
q=OptionValue[Q];
Sq=q^((1-Sign[kk])0.5);
If[((kk==I)||(kk==-I)), Return[{{kk, kk, kk}, {kk,kk, kk}, kk, kk, kk, kk, kk, kk, kk}, Module]];
If[minSep<OptionValue[DMIN],  Return[{{"Coll", "Coll","Coll"}, {"Coll", "Coll","Coll"}, "Coll", "Coll","Coll", "Coll", "Coll","Coll", "Coll"}, Module]];
teval=-xend t0[aa1];
ss=posBound[aa1, \[Phi]1, S, FilterRules[{opts}, Options[posBound]]][[1]];
pos={(r1x/.ss)[teval], (r1y/.ss)[teval],0};
vel={(r1x'/.ss)[teval], (r1y'/.ss)[teval],0};
(*Print[kk];*)
(*Return eccentricity;*)
jc=Cross[pos, vel];
rhat=pos/Norm[pos];
eVec=mr^-1 Cross[vel ,jc]-rhat;
as=1/DD[aa1] mr^(2./3.) (1+q)/Sq  1/(2 Abs[kk]);
abin=(1/DD[aa1]);
ecc=Sqrt[1-(Cross[pos, vel][[3]] DD[aa1]^(1/2) mr^(-1/2) (abin/as)^(1/2))^2];
pomega=ArcTan[eVec[[1]], eVec[[2]]];
nu=ArcTan[ pos[[1]], pos[[2]]]-pomega;

{eVec, jc, nu, as, as/abin, ecc, pomega,  meanAnomaly[nu, ecc], Sign[kk]}
]
(*Orbital elements for the unbound orbit*)
Options[getEcc2]={DMIN->DMINDEF, EVEC->EVECDEF, MR->MRDEF, Q->QDEF};
getEcc2[aa1_, \[Phi]1_, S_, opts:OptionsPattern[]]:=Module[{teval, mr,ss,  pos, vel, kk, jc, rhat,eVec, nu, ecc,pomega, q, as, Sq, ef, FF, abin},
kk=enf[aa1, \[Phi]1, S, FilterRules[{opts}, Options[enf]]];
Print[kk];
mr=OptionValue[MR];
q=OptionValue[Q];
Sq=q^((1-Sign[kk])0.5);
If[((kk==I)||(kk==-I)), Return[{{kk, kk, kk}, {kk,kk, kk}, kk, kk, kk, kk, kk, kk}, Module]];
teval=-xend t0[aa1];
ss=posUnbound[aa1, \[Phi]1, S,opts][[1]];
pos={(r1x/.ss)[teval], (r1y/.ss)[teval],0};
vel={(r1x'/.ss)[teval], (r1y'/.ss)[teval],0};
(*Print[kk];*)
(*Return eccentricity; Implicitly assumes equal mass ratios for now.*)
jc=Cross[pos, vel];
rhat=pos/Norm[pos];
eVec=mr^-1 Cross[vel ,jc]-rhat;
as=1/DD[aa1] mr^(2./3.) (1+q)/Sq  1/(2 Abs[kk]);
abin=(1/DD[aa1]);
ecc=Sqrt[1+(Cross[pos, vel][[3]] DD[aa1]^(1/2) mr^(-1/2) (abin/as)^(1/2))^2];
pomega=ArcTan[eVec[[1]], eVec[[2]]];
nu=ArcTan[ pos[[1]], pos[[2]]]-pomega;
FF=2ArcTanh[((ecc-1)/(ecc+1))^(0.5) Tan[nu/2]];

{eVec, jc, nu, -as, ecc, pomega, ecc*Sinh[FF]-FF}
]

(*Randomly a, \[Phi] parameter space for a given combination of S, e*)
Options[RandSample]={EVEC->EVECDEF};
RandSample[trials_,S_, opts: OptionsPattern[]]:=Module[{DDs, \[Phi]s, enDat, enDatFilt,evec, tag},
DDs=RandomVariate[UniformDistribution[{0.0001, 4}], trials];
\[Phi]s=\[Pi] RandomVariate[UniformDistribution[{-1, 1}], trials];
enDat=Table[enf[Log10[DDs[[i]]], \[Phi]s[[i]],S, opts], {i, 1, trials}];
enDatFilt=Position[enDat,x_/;(( x!=I  ) &&(Abs[x]<100) )];
evec=ToString/@OptionValue[EVEC];
tag=evec[[1]]<>"_"<>evec[[2]]<>"_"<>evec[[3]];
Export[NotebookDirectory[]<>"enSariFull_S"<>ToString[S]<>"_e"<>tag<>".csv" , Transpose[{DDs, \[Phi]s,enDat}]];
Export[NotebookDirectory[]<>"enSari_S"<>ToString[S]<>"_e"<>tag<>".csv", Transpose[{Extract[DDs, enDatFilt], Extract[\[Phi]s, enDatFilt],Extract[enDat,enDatFilt]}]];
]
Options[DisruptFrac]={EVEC->EVECDEF};
DisruptFrac[trials_,aa1_, S_, opts: OptionsPattern[]]:=Module[{DDs, \[Phi]s, enDat, enDatFilt,evec, tag},
\[Phi]s=\[Pi] Range[-1, 1, 2/trials];


enDat=Table[enf[aa1, \[Phi]s[[i]],S, opts], {i, 1, trials}];
enDatFilt=Position[enDat,x_/;(( x!=I  ) &&(Abs[x]<100) )];
evec=ToString/@OptionValue[EVEC];
tag=evec[[1]]<>"_"<>evec[[2]]<>"_"<>evec[[3]];
{aa1, Length[enDatFilt], trials}
(*Export[NotebookDirectory[]<>"enSariFull_S"<>ToString[S]<>"_e"<>tag<>".csv" , Transpose[{DDs, \[Phi]s,enDat}]];
Export[NotebookDirectory[]<>"enSari_S"<>ToString[S]<>"_e"<>tag<>".csv", Transpose[{Extract[DDs, enDatFilt], Extract[\[Phi]s, enDatFilt],Extract[enDat,enDatFilt]}]];*)
]
Options[DisruptFracGrid]={EVEC->EVECDEF};
(*Energy after binary disruption*)
DisruptFracGrid[trials_, S_, opts: OptionsPattern[]]:=Module[{aas, \[Phi]s,evec, tag, tab, ff},
\[Phi]s=\[Pi] Range[-1, 1+0.1*2.0/(trials-1), 2.0/(trials-1)];
aas=Join[Range[-3, -1.19, 0.2], Range[-1, 0.5, 0.05]];

evec=ToString/@OptionValue[EVEC];
tag=evec[[1]]<>"_"<>evec[[2]]<>"_"<>evec[[3]];

tab=Table[{aas[[i]], \[Phi]s[[j]],enfColl[aas[[i]], \[Phi]s[[j]], S, opts]}, {i, 1, Length[aas]}, {j,1, Length[\[Phi]s]}]//Flatten[#,1]&;
tab=Flatten/@tab;
ff=NotebookDirectory[]<>"enGrid_S"<>ToString[S]<>"_e"<>tag<>".csv";
Export[ff, tab, "Table"];
OpenAppend[ff];
WriteString[ff, "\n"];
Close[ff];
]

ppPlot[aa1_, ph1_, S_, opts:OptionsPattern[]]:=Module[{p0,ann, p10,p1, p2, p3, p4,ss, pRange}, 
ss=sol[aa1,  ph1, S, opts];
p0=ListPlot[{{x[t0[aa1]], y[t0[aa1]]}}/.ss, PlotStyle->{Red, PointSize[Large]}];
ann=Graphics[Text[Style["\[Phi]="<>ToString[ph1/\[Pi]]<>" \[Pi]",20]]];
p10=ParametricPlot[{x[-10^(t1)], y[-10^(t1)]}/.ss, {t1, Log10[-ttidal[aa1]],Log10[-t0[aa1]]}, PlotRange->All, Frame->True, FrameLabel->{"X", "Y"}, FrameStyle->Directive[FontSize->20]];
pRange=PlotRange[p10];
p1=ParametricPlot[{x[-10^(t1)], y[-10^(t1)]}/.ss, {t1, Log10[-ttidal[aa1]],Log10[-t0[aa1]]}, PlotRange->All, Frame->True, FrameLabel->{"X", "Y"}, FrameStyle->Directive[FontSize->20],
 Epilog->Inset[ann,{pRange[[1,2]], pRange[[2,2]]}]];

(*Print[pRange];*)
p2=ParametricPlot[{x[-10^(t1)], y[-10^(t1)]}/.ss, {t1, Log10[-ttidal[aa1]], -2}, PlotRange->All, PlotStyle->{Purple, Dashed}];
p3=ParametricPlot[{x[10^(t1)], y[10^(t1)]}/.ss, {t1, -2,Log10[- ttidal[aa1]]}, PlotRange->All, PlotStyle->{Red,Dashed}];

Show[ p1, p2, p3]
]

ppPlot2[aa1_, ph1_, S_, opts:OptionsPattern[]]:=Module[{p0,ann, p10,p1, p2, p3, p4,ss, pRange}, 
ss=sol[aa1,  ph1, S, opts];
p0=ListPlot[{{x[t0[aa1]], y[t0[aa1]]}}/.ss, PlotStyle->{Red, PointSize[Large]}];
ann=Graphics[Text[Style["\[Phi]="<>ToString[ph1/\[Pi]]<>" \[Pi]",20]]];
p10=ParametricPlot[{-10^t1, (x[-10^(t1)]^2+y[-10^(t1)]^2)/.ss}, {t1, Log10[-ttidal[aa1]],Log10[-t0[aa1]]}, PlotRange->All, Frame->True, FrameLabel->{"X", "Y"}, FrameStyle->Directive[FontSize->20]];
pRange=PlotRange[p10];
p1=ParametricPlot[{-10^t1,(x[-10^(t1)]^2+y[-10^(t1)]^2)^0.5/.ss}, {t1, Log10[-ttidal[aa1]],Log10[-t0[aa1]]}, PlotRange->All, Frame->True, FrameLabel->{"X", "Y"}, FrameStyle->Directive[FontSize->20],
 Epilog->Inset[ann,{pRange[[1,2]], pRange[[2,2]]}]];

(*Print[pRange];*)
p2=ParametricPlot[{-10^t1, (x[-10^(t1)]^2+y[-10^(t1)]^2)^0.5/.ss}, {t1, Log10[-ttidal[aa1]], -2}, PlotRange->All, PlotStyle->{Purple, Dashed}];
p3=ParametricPlot[{10^t1, (x[10^(t1)]^2+y[10^(t1)]^2)^0.5/.ss}, {t1, -2,Log10[- ttidal[aa1]]}, PlotRange->All, PlotStyle->{Red,Dashed}];

Show[ p2, p3]
]

Options[pComp]={MR->MRDEF, Q->QDEF, COL->Green};
pComp[aa1_, \[Phi]1_,  S_, opts:OptionsPattern[]]:=Module[{ sol1, sol2, sol1b, sol2b, pa, pb, pc,pd, q, mr, col},
mr=OptionValue[MR];
q=OptionValue[Q];
col=OptionValue[COL];
sol1=r1x/.pos2[aa1,\[Phi]1, S, MR->mr, Q->q][[1]];
sol2=r1y/.pos2[aa1,\[Phi]1 , S, MR->mr, Q->q][[1]];

sol1b=r1x/.pos1[aa1,\[Phi]1, S, MR->mr, Q->q][[1]];
sol2b=r1y/.pos1[aa1,\[Phi]1 , S, MR->mr, Q->q][[1]];

pa=ParametricPlot[{sol1[tt]-sol1b[tt],sol2[tt]-sol2b[tt]}, {tt, t0[aa1], ttidal[aa1]}, PlotRange->All, PlotStyle->{col, Dashed}];
pb=ParametricPlot[{sol1[tt]-sol1b[tt],sol2[tt]-sol2b[tt]}, {tt, ttidal[aa1], 0}, PlotStyle->{col, Dashed}];
pc=ParametricPlot[{sol1[tt]-sol1b[tt],sol2[tt]-sol2b[tt]}, {tt, 0, -ttidal[aa1]}, PlotStyle->{col, Dashed}];
pd=ParametricPlot[{sol1[tt]-sol1b[tt],sol2[tt]-sol2b[tt]}, {tt, -ttidal[aa1], -3 ttidal[aa1]}, PlotStyle->{col, Dashed}];
Show[pa, pb, pc, pd, PlotRange->All]

]

Options[pComp2]={MR->MRDEF, Q->QDEF};
pComp2[aa1_, \[Phi]1_,  S_, opts:OptionsPattern[]]:=Module[{ sol1, sol2, sol1b, sol2b, pa, pb, pc,pd, q, mr, col},
mr=OptionValue[MR];
q=OptionValue[Q];
sol1=r1x/.pos2[aa1,\[Phi]1, S, MR->mr, Q->q][[1]];
sol2=r1y/.pos2[aa1,\[Phi]1 , S, MR->mr, Q->q][[1]];

sol1b=r1x/.pos1[aa1,\[Phi]1, S, MR->mr, Q->q][[1]];
sol2b=r1y/.pos1[aa1,\[Phi]1 , S, MR->mr, Q->q][[1]];

pb=ParametricPlot[{-10^tt, Sqrt[(sol1[-10^tt]-sol1b[-10^tt])^2+(sol2[-10^tt]-sol2b[-10^tt])^2]}, {tt, Log10[-ttidal[aa1]], -2}, PlotStyle->{Purple}];
pc=ParametricPlot[{10^tt, Sqrt[(sol1[10^tt]-sol1b[10^tt])^2+(sol2[10^tt]-sol2b[10^tt])^2]}, {tt, -2,  Log10[-ttidal[aa1]]}, PlotStyle->{Red}];
(*pd=Plot[Sqrt[(sol1[tt]-sol1b[tt])^2+(sol2[tt]-sol2b[tt])^2], {tt, -ttidal[aa1], -3 ttidal[aa1]}, PlotStyle\[Rule]{Red, Dashed}];*)
Show[pb, pc]
]

PGOAL=Automatic;
AGOAL=Automatic; 
STEPS=3 10^5;
xend=50;
End[]
EndPackage[]







