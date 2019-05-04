(* ::Package:: *)

Needs["Hill`"]
BeginPackage["HardSoft`"]
MC::usage="MC[Ntrials,  opts]--Run Ntrials binary disruption experiments with central black hole of mass Mbh"


Needs["Hill`"]


Begin["`Private`"]

pc=3.08568*10^18;
Msun = 1.98892*10^33;
year = 3.15569 10^7;
G = 6.674*10^-8;
Rsun=6.957*10.^10;
au=1.495978707*10.^13;
c=2.99792458 10^10;

Mbh=4*10^6 Msun;
to=4*10^6  year;
no=1.5*10^6;
ro=0.1 pc;
hr=1;
MMIN=0.1;

loc="/home/aleksey/software/MIST_v1.1_feh_p0.00_afe_p0.0_vvcrit0.4_EEPS/";
rdat=Import[loc<>"r_zams.csv"];
tdat=Import[loc<>"mt0.csv"];
rinterp=rdat//Interpolation;
Minterp=Log10[tdat]//Interpolation[#, InterpolationOrder->1]&;
Mt0[tt_]:=10^Minterp[Log10[tt]]
(*muKroupa[M_]:=Piecewise[{{(M)^(-23/10), M>=1/2},{(1/2)^-1 (M)^(-13/10),(8/100<=M)&&(M<1/2)}},(1/2)^-1  (8/100)^-1 (M)^(-3/10)];*)
muKroupa1[M_]:= Piecewise[{{(M)^-2.3, M>=0.5},{0.5^-1. (M)^-1.3,(0.08<=M)&&(M<0.5)}},0.5^-1 0.08^-1 ( M)^-0.3];
muKroupa[M_]:=muKroupa1[M]/NIntegrate[muKroupa1[m1],{m1, MMIN, 100}]
IMF:=ProbabilityDistribution[muKroupa[M], {M, MMIN, 100}];
PDMF:=ProbabilityDistribution[muKroupa[M], {M, MMIN, Mt0[to/year]}];

(*Mean mass and RMS mass should evolv self-consistently--see below.!*)
norm[mmax_]:=Piecewise[{{(250*mmax^(7/10))/7,mmax<=2/25},{(100*2^(7/10)*5^(3/5))/21-20/(3*mmax^(3/10)),Inequality[2/25,Less,mmax,LessEqual,1/2]}},(-200*2^(3/10))/39+(100*2^(7/10)*5^(3/5))/21-10/(13*mmax^(13/10))]-7.016379835666912
mbarf[mmax_]:=(Piecewise[{{(250*mmax^(17/10))/17,mmax<=2/25},{(-8*2^(7/10)*5^(3/5))/119+(20*mmax^(7/10))/7,Inequality[2/25,Less,mmax,LessEqual,1/2]}},(100*2^(3/10))/21-(8*2^(7/10)*5^(3/5))/119-10/(3*mmax^(3/10))]-0.28323077863336193)/norm[mmax]
M2f[mmax_]:=(Piecewise[{{(250*mmax^(27/10))/27,mmax<=2/25},{(-16*2^(7/10))/(2295*5^(2/5))+(20*mmax^(17/10))/17,Inequality[2/25,Less,mmax,LessEqual,1/2]}},(-50*2^(3/10))/119-(16*2^(7/10))/(2295*5^(2/5))+(10*mmax^(7/10))/7]-0.017524313761887154)/norm[mmax]
M2:=M2f[Mt0[to/year]]Msun^2;
mbar:=mbarf[Mt0[to/year]] Msun;


(*M2=0.13 Msun^2;
mbar=0.3 Msun;*)
vkep:=Sqrt[G Mbh/ro];
sig:=hr vkep;
n:=no/pc^3 hr^-1;
ahs[m_,q_]:=G (q m^2)/( (1+q)^2 sig^2 mbar);
agw[m_,  q_, e_]:= (5/(64*4) (c^5 (1+q)^2)/(G^3 m^3 to q) (1-e^2)^(7/2)/(1+73/24 e^2+37/96 e^4))^(-1/4)
v12[m_, a_, q_]:=Sqrt[G  m/a]
lam12[m_, a_, q_]:=Log[2 sig^2/v12[m,a, q]^2]
aevapRoot[m_?NumericQ,  q_]:=NSolve[0.07 (v12[m, a, q]^2 sig)/(G^2  n M2 lam12[m,a, q])==to, a]
aevap[m_?NumericQ, q_]:=(a/.aevapRoot[m,  q])
rRoche[q_]:=(0.49 q^(2/3)/(0.6 q^(2/3)+Log[1+q^(1/3)]));
acontact[m_,q_]:=Max[rinterp[m/Msun/(1+q)]Rsun/rRoche[q^-1], rinterp[m/Msun q/(1+q)]Rsun/rRoche[q]];


fbin[a1_, amin_, amax_]:=Piecewise[{{a1^-1/Log[amax/amin], a1<=amax && a1>=amin}},0];
powSamp[p_, amin_, amax_]:=Module[{r},
r=RandomVariate[UniformDistribution[{0,1}]];
Piecewise[{{amin*Exp[r*Log[amax/amin]], p==1}}, (r*(amax^(1-p)-amin^(1-p))+amin^(1-p))^(1./(1-p))]
]
Cstring[x_, p_:2]:=ToString[CForm[SetPrecision[x,p]]]
InBetween[ords_]:=Table[0.5(ords[[i+1]]+ords[[i]]), {i,1, Length[ords]-1}];

Options[smaEmpirical]={COMPACT->False, COLL->False, ECC->0};
Options[initConditionsDim]={Q->Hill`Private`QDEF, EVEC->Hill`Private`EVECDEF, MR->Hill`Private`MRDEF};
initConditionsDim[a1_, \[Phi]_, S_, abin_,m_, opts:OptionsPattern[]]:=Module[{m0, abinAU, toDim, i1, rscale, tscale}, 
i1=Hill`initConditions[a1, \[Phi], S, opts];
abinAU=abin/au;
m0=m/Msun;
rscale=(10^a1 abinAU);
tscale=10^(1.5 a1) Sqrt[abinAU^3/(m0)];
toDim= rscale{1, 1, 1, tscale^-1, tscale^-1, tscale^-1};
{i1[[1]] toDim, i1[[2]] toDim}
]
Options[smaEmpirical]={ECC->0, COLL->False, COMPACT->False};
smaEmpirical[opts: OptionsPattern[]]:=Block[{\[Phi]1, a1, q,m, m1, m2,rt, smas, kpro, kretro, ks, a0, filt, smas2, amin, amax, comp, out,aa1,
delta, dmin, ecc, S1, i1, ss, endState, cut, rscale, abinAU, tscale},
m1=RandomVariate[PDMF] Msun;
m=m1(1+q);
m2=q m1;
q=RandomVariate[UniformDistribution[{MMIN Msun/( m1),1}]];

ecc=OptionValue[ECC];
If[OptionValue[ECC]=="Thermal", ecc=RandomVariate[UniformDistribution[{0,1}]]^0.5];
If[OptionValue[ECC]=="Uniform", ecc=RandomVariate[UniformDistribution[{0,1}]]];
(*Cap sma distribution at 100 au...*)
{amin, amax}={Max[acontact[m, q]/(1-ecc), agw[m,  q,ecc]],Min[Max[(aevap[m, q]), ahs[m, q]], 100 au]};
comp=OptionValue[COMPACT];
If[comp, {amin, amax}={agw[m, q, ecc],Min[Max[(aevap[m, q]), ahs[m, q]], 100 au]}];
If[amin>=amax, Return["ERR", Module]];
smas=powSamp[1, amin, amax];
dmin=0;
If[And[OptionValue[COLL], Not[comp]], dmin=acontact[m, q]/smas];
rt=(Mbh/m)^(1/3) smas;
cut=Max[(Mbh/m1)^(1/3) rinterp[m1/Msun] Rsun/(rt), (Mbh/m2)^(1/3) rinterp[m2/Msun] Rsun/(rt), 6 G Mbh/c^2/rt];

aa1=RandomVariate[UniformDistribution[{cut, 3.2}]]//Log10;
\[Phi]1=RandomVariate[UniformDistribution[{-\[Pi], \[Pi]}]];
S1= RandomChoice[{-1,1}];

ss=Hill`sol[aa1, \[Phi]1, S1, EVEC->{ecc,0,0}];
endState={Hill`Private`x[-(Hill`Private`xend)Hill`t0[aa1]], Hill`Private`y[-(Hill`Private`xend)Hill`t0[aa1]], Hill`Private`x'[-(Hill`Private`xend)Hill`t0[aa1]], Hill`Private`y'[-(Hill`Private`xend)Hill`t0[aa1]]}/.ss;
abinAU=smas/au;
rscale=10^aa1 abinAU;
endState=endState rscale;
tscale=10^(1.5 aa1) Sqrt[abinAU^3/(m/Msun)];
endState[[3;;]]=endState[[3;;]]/tscale;

i1={{m/Msun, q, -t0[aa1] tscale},initConditionsDim[aa1, \[Phi]1, S1, smas, m, EVEC->{ecc,0,0}, Q->q, MR->Mbh/m]}//Flatten;
out= Hill`getEcc[aa1, \[Phi]1,S1, EVEC->{ecc,0,0}, Q->q, DMIN->dmin, MR->Mbh/m ][[5;;]];
{i1, Flatten[{10^aa1, \[Phi]1, S1, m/Msun, smas/au, ecc, q, out}], endState}
]

Options[MC]=Options[smaEmpirical];
MC[Ntrials_, opts:OptionsPattern[]]:=Module[{dat, dat1, dat2, dat3, names, base, out},
dat=Table[smaEmpirical[ opts], {i, 1, Ntrials}];
Print[Position[dat, x_/;x=="ERR"]//Flatten];
dat=dat[[Complement[Range[1,Length[dat]],Position[dat, x_/;x=="ERR"]//Flatten]]];
names= { "m", "q", "t0", "x", "y" , "z","vx" , "vy", "vz", "x2", "y2" , "z2","vx2" , "vy2", "vz2"};
dat1=dat[[;;,1]];
dat1=Prepend[dat1, names];

names= { "rp/rt", "phi", "S", "m" , "abin","e", "q", "as/abin", "es", "os", "Ms",  "Sk"};
dat2=dat[[;;,2]];
dat2=Prepend[dat2, names];
names= { "x", "y", "x'", "y'"};
dat3=dat[[;;,3]];
dat3=Prepend[dat3, names];
base="Mbh"<>Cstring[Mbh/Msun]<>"_t"<>Cstring[to/year] <>"_e"<>ToString[OptionValue[ECC]]<>"_coll" <>ToString[OptionValue[COLL]]<>".csv";
out=NotebookDirectory[]<>"trial_init_"<>ToString[base];
Export[out, dat1];
out=NotebookDirectory[]<>"trial_"<>ToString[base];
Export[out, dat2];
out=NotebookDirectory[]<>"trial_end_"<>ToString[base];
Export[out, dat3];
]


End[]
EndPackage[]
