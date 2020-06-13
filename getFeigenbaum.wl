(* ::Package:: *)

BeginPackage["getFeigenbaum`"]
(* calculate Feigenbaum constant and point *)

getFeigenbaum::usage = "get the n first superstable points and increasingly accurate Feigenbaum constant"
getFeigenbaumPoint::usage = "get an approximation of the Feigenbaum point"
getFeigenbaumConstant::usage = "get an approximation of the Feigenbaum constant"
getSuperstablePoints::usage = "get the n first superstable points"
getScalingFactor::usage = "get an approximation of the scaling factor"
getIterate::usage = "get value and derivative of [2^(n-1)]th iterate of f"
getMaxPoint::usage = "get the point where f achieves its maximum"

Begin["`Private`"]
(* approximate nth bifurcation point from previous ones using estimate of delta *)
approximateS[sl_,n_]:=Module[{},

N[sl[[n-1]]+(sl[[n-1]]-sl[[n-2]])/dl[[n-1]]]
]

(* get value and derivative of [2^(n-1)]th iterate of f *)
getIterate[f_,r_,n_]:=Module[{x,xMax,dr=0},

xMax=getMaxPoint[f];
x=xMax;

For[i=1,i<=2^(n - 1),i++,

(* here we are calculating the derivative iteratively *)
df[t_,y_]=Dt[f[t,y],t]/.Dt[y,t]->dr;
dr=N[Re[ComplexExpand[df[r,x]]]];

x=N[f[r,x]];
];

{x-xMax,dr}
]

(* get the point where f achieves its maximum *)
getMaxPoint[f_]:=Module[{},
x/.Last[Quiet[FindMaximum[f[1,x],{x,0}]]]
]

(* get an estimate of Feigenbaum constant from the superstable points *)
getDelta[sl_,n_]:=Module[{},
(sl[[n-1]]-sl[[n-2]])/(sl[[n]]-sl[[n-1]])
]

(* get the n first superstable points and increasingly accurate Feigenbaum constant
	f: the unimodal function
	{rMin_,rMax_}: the interval of f
	d0: an estimate of the feigenbaum constant
	snMax: the number of superstable points used for the approximation
*)
getFeigenbaum[f_,{rMin_,rMax_},d0_:0,snMax_:12]:=Module[{},

xMax=getMaxPoint[f];

f4[r_]:=Nest[f[r,#]&,xMax,2^2];

(* solve for the first 3 superstable points analytically to obtain a first estimate for delta *)
sn3=Quiet[Solve[f4[r]==xMax&&rMin<=r<=rMax,r,Reals]];

(* if the maximum of f is at 0 we cannot determine enough superstable points analytically as to approximate delta. Thus an estimate of delta needs to be passed as argument *)
If[xMax==0&&d0==0,
Throw["Not enough intial superstable points could be determined analytically."]
];

(* buffer the list of superstable points with 0 if not enough  could be determined *)
If[xMax==0,
sl=Prepend[Take[Values[Flatten[sn3]],2],0],
sl=Take[Values[Flatten[sn3]],3]
];

If[d0==0,
dl=ConstantArray[getDelta[sl,3],3],
dl=ConstantArray[d0,3]
];

For[k=4,k<=snMax,k++,

sn={approximateS[sl,k]};

For[j=1,j<1000,j++,

{xn,drn}=getIterate[f,sn[[j]],k];

(* avoid division by 0 *)
If[drn==0,Break[]];

(* this is the actual Newton iteration *)
AppendTo[sn,sn[[j]]-xn/drn];

If[sn[[j+1]]==sn[[j]],Break[]];
];

AppendTo[sl,Last[sn]];

AppendTo[dl,getDelta[sl,k]];
];

Transpose[{sl,dl}]
]

(* get an approximation of the Feigenbaum point
	f: the unimodal function
	{rMin_,rMax_}: the interval of f
	d0: an estimate of the feigenbaum constant
	snMax: the number of superstable points used for the approximation
*)
getFeigenbaumPoint[f_,{rMin_,rMax_},d0_:0,snMax_:12]:=Module[{},
Last[getFeigenbaum[f,{rMin,rMax},d0,snMax]][[1]]
]

(* get an approximation of the Feigenbaum constant
	f: the unimodal function
	{rMin_,rMax_}: the interval of f
	d0: an estimate of the feigenbaum constant
	snMax: the number of superstable points used for the approximation
*)
getFeigenbaumConstant[f_,{rMin_,rMax_},d0_:0,snMax_:12]:=Module[{},
Last[getFeigenbaum[f,{rMin,rMax},d0,snMax]][[2]]
]

(* get the snMax first superstable points
	f: the unimodal function
	{rMin_,rMax_}: the interval of f
	d0: an estimate of the feigenbaum constant
	snMax: the number of superstable points
*)
getSuperstablePoints[f_,{rMin_,rMax_},d0_:0,snMax_:12]:=Module[{},
getFeigenbaum[f,{rMin,rMax},d0,snMax][[All,1]]
]

(* get an approximation of the scaling factor
	f: the unimodal function
	{rMin_,rMax_}: the interval of f
	d0: an estimate of the feigenbaum constant
	snMax: the number of superstable points used for the approximation
*)
getScalingFactor[f_,{rMin_,rMax_},d0_:0,snMax_:12]:=Module[{mu,xMax,d,a,l},
	mu=getSuperstablePoints[f,{rMin,rMax},d0,snMax];
	xMax=getMaxPoint[f];
	d={0};
	a={};
	
	For[i=2,i<=Length[mu],i++,
	l=NestList[N[f[mu[[i]],#]]&,xMax,2^(i-1)-1];
		AppendTo[d,Min[Abs[Drop[l,1]-xMax]]];
		AppendTo[a,d[[i-1]]/d[[i]]];
	];
	Last[a]
];

End[]
EndPackage[]



