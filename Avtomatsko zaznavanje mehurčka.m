(* ::Package:: *)

(* ::Subsubsection:: *)
(*Funkcija Houghove transformacije*)


(* ::Input:: *)
(*HoughCircleDetection[image_Image,radiusmin_Integer: 1,radiusmax_Integer: 40,edgedetectradius_Integer: 10,minfitvalue_Real: .25,radiusstep_Integer: 1,minhoughvoxels_Integer: 4]:=Module[{edgeimage,hough3dbin,hough3dbinlabels,coords,arraydim},edgeimage=SelectComponents[DeleteBorderComponents[EdgeDetect[image,edgedetectradius,Method->"Canny"]],"EnclosingComponentCount",#==0&];*)
(*hough3dbin=DeleteSmallComponents[Image3D[ParallelMap[Binarize[Image@Divide[ListConvolve[#,ImageData@edgeimage,Ceiling[(Length@#)/2]],Total[#,2]],minfitvalue]&,Map[Function[{r},DiskMatrix[r]-ArrayPad[DiskMatrix[r-1],1]][#]&,Range[radiusmin,radiusmax,radiusstep]]]],minhoughvoxels];*)
(*hough3dbinlabels=MorphologicalComponents[hough3dbin];*)
(*coords=ParallelMap[Round[Mean[Position[hough3dbinlabels,#]]]&,Sort[Rest@Tally@Flatten@hough3dbinlabels,#1[[2]]>#2[[2]]&][[All,1]]];*)
(*arraydim=Rest@Dimensions[hough3dbinlabels];*)
(*(*Print["Radii: ",Max[radiusmin+coords[[All,1]]-2]];*)*)
(*{Max[radiusmin+coords[[All,1]]-2],*)
(*ParallelMap[Function[{level,offx,offy},ImageMultiply[image,Image@ArrayPad[DiskMatrix[radiusmin+level-1],{{offx-radiusmin-level,First@arraydim-offx-radiusmin-level+1},{offy-radiusmin-level,Last@arraydim-offy-radiusmin-level+1}}]]][Sequence@@#]&,coords]}]*)


(* ::Text:: *)
(*Uvoz slik*)


(* ::Input:: *)
(*{FileNameSetter[Dynamic[f], "Directory"], Dynamic[f]}*)


(* ::Input:: *)
(*files=FileNames["*.jpg",f];*)
(*all=Import[#]&/@files;*)
(*frejmi=Table[ImageCrop[all[[i]]],{i,Length[all]}];*)


(* ::Input:: *)
(*ms=frejmi[[30]];*)
(*markschulzeedges=EdgeDetect[ms,10,Method->"Canny"]*)
(*ParallelMap[Image@Divide[ListConvolve[#,ImageData@markschulzeedges,Ceiling[(Length@#)/2]],Total[#,2]]&,Map[Function[{r},DiskMatrix[r]-ArrayPad[DiskMatrix[r-1],1]][#]&,Range[14,18,1]]];*)
(**)
(*j=Show[ImageApply[Plus,HoughCircleDetection[#,8,50,7,.30][[2,1]]],ImageSize->ImageDimensions[#]]&[ms]*)
(**)
(*X=ComponentMeasurements[j,{"Centroid","EquivalentDiskRadius"},#AdjacentBorderCount==0&&400<#Area<50000&];*)
(*HighlightImage[ms,Circle@@@X[[All,2]]]*)


(* ::Subsubsection:: *)
(*Preverjanje*)


(* ::Input:: *)
(*img2=ms;*)
(*measurementCircle[Dynamic[center_],Dynamic[r_]]:={Red,Thick,Dynamic@Circle[center,r],*)
(*Dynamic@*)
(*Text[Style[StringJoin[ToString@Round@r,"px"],FontSize->18],Scaled@{0.1,0.9},Background->RGBColor[0,0,0,0.5]]};*)
(**)
(*DynamicModule[*)
(*{center=ImageDimensions[img2]/2,r},Manipulate[*)
(*EventHandler[*)
(*DynamicImage[*)
(*img2,*)
(*Epilog->measurementCircle[Dynamic[center],Dynamic[r]]*)
(*],{"MouseDown":>If[CurrentValue["OptionKey"],center=MousePosition["Graphics"]]},PassEventsDown->Dynamic[Not[CurrentValue["OptionKey"]]]*)
(*],*)
(*{{r,32},8,800},*)
(*FrameMargins->0*)
(*]*)
(*]*)


(* ::Input:: *)
(*korak=15;*)
(*Timing[R=Table[If[HoughCircleDetection[frejmi[[i]],8,50,7,.30][[1]]==-\[Infinity],0,HoughCircleDetection[frejmi[[i]],8,50,7,.30][[1]]],{i,1,Length[frejmi],korak}];]*)


(* ::Input:: *)
(*pxmmraz=0.18/(2Sqrt[2]);*)
(**)
(*Rt=Table[{10^-6 (i-1)3.3333 korak ,R[[i]] pxmmraz},{i,Length[R]}];*)
(*Rp={0,21,28,33,37,39,40,41,42,43,43,41,41,40,37,33,28,21,8,8,9,10} ;*)
(*Rtp=Table[{10^-6 (i-1)3.333 10 ,Rp[[i]]pxmmraz},{i,Length[Rp]}];*)
(*Rpyt=Import["C:\\Users\\Jure Podobnikar\\Desktop\\Diploma\\Algoritmi\\Python\\test.txt","Table"]//Flatten;*)
(*Rpyt[[{1,2,3,4,5,6,7}]]={0,5,7,9,13,15,18};*)
(*Rtpyt=Table[{10^-6 (i-1)3.3333 5 ,Rpyt[[i]]pxmmraz},{i,Length[Rpyt]}]*)


(* ::Input:: *)
(*p1=ListLinePlot[Rt,AxesLabel->{"t[\[Mu]s]","R[mm]"},PlotStyle->Orange,Mesh->All,PlotLegends->{"Mathematica"}];*)
(*p2=ListLinePlot[Rtp,AxesLabel->{"t[\[Mu]s]","R[mm]"},PlotStyle->Blue,Mesh->All,PlotLegends->{"Ro\[CHacek]no"}];*)
(*p3=ListLinePlot[Rtpyt,AxesLabel->{"t[\[Mu]s]","R[mm]"},PlotStyle->Green,Mesh->All,PlotLegends->{"Python"}];*)
(*Show[p1,p3,PlotRange->All,ImageSize->Full]*)


(* ::Input:: *)
(*Print["Odstopanje v Mathematici je ",N[Total[Table[{Abs[R[[i]]-Rp [[i]]]},{i,1,Length[Rp]}]//Flatten]/Length[Rp] pxmmraz,1],"mm."]*)
(*Print["Odstopanje v Pythonu je ",NumberForm[N[Total[Table[{Abs[Rpyt[[i]]-Rp [[i]]]},{i,1,Length[Rp]}]//Flatten]/Length[Rp] pxmmraz],1],"mm."]*)


(* ::Input:: *)
(*Vhodni podatki*)


(* ::Input:: *)
(*Tinf1=Quantity[293.15,"Kelvins"] (*Temperatura vode*)*)
(*Tinf=QuantityMagnitude[Tinf1];*)
(**)
(*pv1=Quantity[Exp[20.386-5132/Tinf]133.322368,"Pascals"](*Tlak nasi\[CHacek]enja*)*)
(*pv=QuantityMagnitude[pv1];*)
(**)
(*\[Sigma]1=Quantity[7.28 10^-2,"Newtons"/"Meters"](*Povr\[SHacek]inska napetost*)*)
(*\[Sigma]=QuantityMagnitude[\[Sigma]1];*)
(**)
(*R01=Quantity[2 10^-4,"Meters"](*Za\[CHacek]etni polmer mehur\[CHacek]ka*)*)
(*R0=QuantityMagnitude[R01];*)
(**)
(*\[Rho]l1=ThermodynamicData["Water","Density",{"Pressure"->Quantity[1000,"Millibars"],(*Gostota vode v kapljevinastem stanju*)"Temperature"->Quantity[QuantityMagnitude[Tinf],"Kelvins"]}]*)
(*\[Rho]l=QuantityMagnitude[\[Rho]l1];*)
(**)
(*\[Rho]v1=ThermodynamicData["Water","Density",{"Temperature"->Quantity[101,"DegreesCelsius"]}](*Gostota vodne pare*)*)
(*\[Rho]v=QuantityMagnitude[1.2];*)
(**)
(*\[Eta]1=ThermodynamicData["Water","Viscosity",{"Pressure"->Quantity[1000,"Millibars"],"Temperature"->Quantity[Tinf,"Kelvins"]}](*Viskoznost vode*)*)
(*\[Eta]=QuantityMagnitude[\[Eta]1];*)
(**)
(*cpl1=ThermodynamicData["Water","IsobaricHeatCapacity",{"Pressure"->Quantity[1000,"Millibars"],"Temperature"->Quantity[QuantityMagnitude[Tinf],"Kelvins"]}]*)
(*(*Specifi\[CHacek]na izobarna toplotna kapacitivnost*)*)
(*cpl=QuantityMagnitude[cpl1];*)
(**)
(*\[Gamma]=1.4(*Razmerje specifi\[CHacek]nih toplot*)*)
(**)
(*\[Lambda]l1=ThermodynamicData["Water","ThermalConductivity",{"Pressure"->Quantity[1000,"Millibars"],"Temperature"->Quantity[Tinf,"Kelvins"]}]*)
(*\[Lambda]l=QuantityMagnitude[\[Lambda]l1];*)
(**)
(*L1=Quantity[2256.9 10^3,"J/kg"]*)
(*L=QuantityMagnitude[L1];*)
(**)
(*t01=Quantity[2 10^-3,"Seconds"]*)
(*t0=QuantityMagnitude[t01];*)
(**)
(*Thermo=0;*)


(* ::Input:: *)
(*Robni pogoji*)


(* ::Input:: *)
(*\[Omega]1=Quantity[23.4 10^3,"Hertz"](*Frekvenca*)*)
(*\[Omega]=QuantityMagnitude[\[Omega]1];*)
(**)
(*p01=Quantity[1.01 10^5,"Pascals"](*Za\[CHacek]etni okoli\[SHacek]ki tlak*)*)
(*p0=QuantityMagnitude[p01];*)
(**)
(*pf1=Quantity[1 1.28 10^5,"Pascals"](*Vzbujalni tlak*)*)
(*pf=QuantityMagnitude[pf1];*)
(**)
(*pg01=p01-pv1+(2\[Sigma]1)/R01*)
(*pg0=QuantityMagnitude[pg01];*)
(**)
(*\[Alpha]l1=\[Lambda]l1/(\[Rho]l 1cpl1)*)
(*\[Alpha]l=QuantityMagnitude[\[Alpha]l1];*)
(**)
(*sum1=(\[Rho]v1 L1)^2/(\[Rho]l1^2 cpl1 Tinf1 Sqrt[\[Alpha]l1]) Thermo;*)
(*sum=(\[Rho]v L)^2/(\[Rho]l^2 cpl Tinf Sqrt[\[Alpha]l]) Thermo*)
(**)
(*pinf[t_]=p0-pf Sin[2 \[Pi] \[Omega] t]*)
(**)


(* ::Input:: *)
(*DOPRIamat={{1/5},{3/40,9/40},{44/45,-56/15,32/9},{19372/6561,-25360/2187,64448/6561,-212/729},{9017/3168,-355/33,46732/5247,49/176,-5103/18656},{35/384,0,500/1113,125/192,-2187/6784,11/84}};*)
(*DOPRIbvec={35/384,0,500/1113,125/192,-2187/6784,11/84,0};*)
(*DOPRIcvec={1/5,3/10,4/5,8/9,1,1};*)
(*DOPRIevec={71/57600,0,-71/16695,71/1920,-17253/339200,22/525,-1/40};*)
(*DOPRICoefficients[5,p_]:=N[{DOPRIamat,DOPRIbvec,DOPRIcvec,DOPRIevec},p];*)


(* ::Input:: *)
(*Clear[pf,\[Omega]]*)
(*solns=ParametricNDSolveValue[{r''[t]==1/(\[Rho]l r[t]) (-((3(r'[t])^2)/2)\[Rho]l-r'[t] sum Sqrt[t] \[Rho]l+pv-p0+pf Sin[2 \[Pi] \[Omega] t]+pg0 (R0^3/r[t]^3)^\[Gamma]-2 \[Sigma]/r[t]-4 \[Eta] r'[t]/r[t]),r[0]==R0,r'[0]==0},{r},{t,0,t0},{pf,\[Omega]},Method->{"ExplicitRungeKutta","DifferenceOrder"->5,"Coefficients"->DOPRICoefficients,"StiffnessTest"->False},MaxStepSize->2 10^-8];*)


(* ::Input:: *)
(*Manipulate[Plot[Evaluate@Through[solns[pf,\[Omega]][t]],{t,0,a t0}],{{pf,0,"Vzbujevalni tlak"},0,10 10^6,Appearance->"Labeled"},{{\[Omega],0,"Vzbujevalna frekvenca"},0,2 10^3,Appearance->"Labeled"},{{a,0.1,"Prikaz"},0,1.,Appearance->"Labeled"}]*)


(* ::Input:: *)
(*s3=NDSolve[{r''[t]==1/(\[Rho]l r[t]) (-((3(r'[t])^2)/2)\[Rho]l-r'[t] sum Sqrt[t] \[Rho]l+pv-p\[Infinity][t]+pg0 (R0/r[t])^(3\[Gamma])-2 \[Sigma]/r[t]-4 \[Eta] r'[t]/r[t]),r[0]==R0,r'[0]==0},r,{t,0,t0},Method->{"ExplicitRungeKutta","DifferenceOrder"->5,"Coefficients"->DOPRICoefficients,"StiffnessTest"->False},MaxStepSize->2 10^-8];*)
(**)
(*r3[t_]=r[t]/.s3[[1]];*)
(*p4=Plot[r3[t],{t,0,0.1t0},AxesLabel->{"t[s]","R[m]"},AxesLabel->{"t[\[Mu]s]","R[px]"},PlotStyle->Red,PlotLegends->{"RPE"}];*)
(*Show[p4,PlotRange->All,ImageSize->Full]*)


(* ::Input:: *)
(**)
(*Clear[R0,pg0]*)
(*pg0=p0-pv+(2\[Sigma])/R0*)
(*pinf[t_]=p0-pf Sin[2 \[Pi] \[Omega] t]*)
(*Rmax=6.08112 10^-3*)
(*tmax=(0.00074997+0.0008333)/2*)


(* ::Input:: *)
(*NSolve[0==1/(\[Rho]l Rmax) (pv-pinf[tmax]+pg0 (R0^3/Rmax^3)^\[Gamma]-(2\[Sigma])/Rmax),R0]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Solve[101000-2338==0.837 \[Rho]l x^2/0.0026^2]*)


(* ::Input:: *)
(* ((2\[Pi])/0.000150)*)


(* ::Input:: *)
(*(3 5500 10^-6 60^2)/(8 \[Pi](2 10^-6)) 0.30*)


(* ::Input:: *)
(*((1.83 10.88 10^-3 Sqrt[\[Rho]l/p0])/(2 \[Pi]))^-1*)


(* ::Input:: *)
(*sin=Interpolation[Rtpyt]*)


(* ::Input:: *)
(*Fit[Rtpyt,{1,x,x^2,x^3,x^4,x^5},x]*)


(* ::Input:: *)
(*p5=Plot[%,{x,0,0.001}]*)
(*Show[p3,p4,p5,PlotRange->All,ImageSize->Large,PlotRange->All]*)


(* ::Input:: *)
(*p\[Infinity][t_]=p0-3/(8 \[Pi]) (6800 (10^-6) 60 Exp[-(t/(0.139 6800 10^-6))] )/(2 10^-6)*)


(* ::Input:: *)
(*Plot[p\[Infinity][t],{t,0,t0}]*)


(* ::Input:: *)
(*Solve[0.001==0.915 r0 Sqrt[\[Rho]l/(p0-pv)]]*)


(* ::Input:: *)
(*\[Omega]0=Sqrt[3\[Gamma] pg0/(\[Rho]l R0^2)-(2\[Sigma])/(\[Rho]l R0^3)]*)
