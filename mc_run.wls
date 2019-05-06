#!/usr/bin/env wolframscript
(* ::Package:: *)

Needs["HardSoft`"]
$HistoryLength=0;

getParam[param_?StringQ]:=Module[{},
	p=Select[$ScriptCommandLine, StringMatchQ[#, RegularExpression[param<>"=.*"]] &];
	p=StringSplit[p, "="]//Flatten//#[[-1]]&;
	p=Read[StringToStream[p], Number];
	p
]
to1=getParam["to"];
Mbh1=getParam["Mbh"];
no1=getParam["no"];
ro1=getParam["ro"];

HardSoft`Private`to=to1 HardSoft`Private`year
HardSoft`Private`ro=ro1 HardSoft`Private`pc
HardSoft`Private`no=no1;
HardSoft`Private`Mbh=Mbh1 HardSoft`Private`Msun;
Print[HardSoft`Private`to]
Print[HardSoft`Private`ro]
Print[HardSoft`Private`no]
Print[HardSoft`Private`Mbh]

Ntrials=2*10^4;
HardSoft`MC[Ntrials, ECC->"Thermal"]
HardSoft`MC[Ntrials, ECC->"Thermal", COLL->True]
HardSoft`MC[Ntrials, ECC->"Uniform"]
HardSoft`MC[Ntrials, ECC->"Uniform", COLL->True]
HardSoft`MC[Ntrials, ECC->0]
HardSoft`MC[Ntrials, ECC->0, COLL->True]
