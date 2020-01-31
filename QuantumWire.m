(* ::Package:: *)

(* ::Title:: *)
(*Quantum Wire*)


(* ::Section:: *)
(*Package Information*)


(* ::Text:: *)
(*This package uses the Fourier-transformed version of the Schr\[ODoubleDot]dinger equation to compute the eigenvalues and eigenfunctions over rectangular regions, with random or ordered delta-function potentials\[LongDash]corresponding to a 2D quantum wire\[LongDash]in Cartesian coordinates.*)


(* ::Text:: *)
(*A more efficient program might use the built-in discrete Fourier Sin Transform.*)


(* ::Section:: *)
(*Package Beginnings*)


BeginPackage["QuantumWire`"];


(* ::Text:: *)
(*Clear pre-existing definitions:*)


Unprotect["QuantumWire`*"];


ClearAll["QuantumWire`*"];


(* ::Subsection:: *)
(*Usage statements*)


MatrixDisplay::usage="MatrixDisplay[m] displays the matrix m as a GridBox.";


LabeledMatrixDisplay::usage="LabeledMatrixDisplay[l,m] displays the matrix m with corresponding index labels l."; 


PotentialDiagram::usage="PotentialDiagram[{Lx,Ly}][spikes] displays the locations of the set of \[Delta]-function spikes {x,y} in the rectangle Lx \[Times] Ly."; 


Indices::usage="IndexSet[n] constructs the set of pairs of indices from 1 to n.";


\[CapitalDelta]Term::usage="\[CapitalDelta]Term[{Lx,Ly}][{x,y}][{n,m}] gives a term in the potential matrix for a \[Delta]-function spike at {x,y} in the rectangle Lx \[Times] Ly.";


\[CapitalDelta]Vector::usage="\[CapitalDelta]Vector[{Lx,Ly}][spikes][indices] constructs the indexed vector for all indices of the \[Delta]-function spike locations in the rectangle Lx \[Times] Ly.";


PotentialMatrix::usage="PotentialMatrix[{Lx,Ly}][spikes][indices] computes the potential matrix for the set of \[Delta]-function spike locations in the rectangle Lx \[Times] Ly using an Outer product involving Dot products of vectors."; 


HamiltonianMatrix::usage="HamiltonianMatrix[{Lx,Ly},\[Alpha]][spikes][indices] computes the Hamiltonian matrix with coupling constant \[Alpha] and \[Delta]-function spikes in the rectangle Lx \[Times] Ly."; 


EnergyOrdered::usage="EnergyOrdered[es] re-orders the eigensystem in order of increasing energy."; 


QuantumWireSystem::usage="QuantumWireSystem[{Lx,Ly},\[Alpha]][spikes][indices,n] determines the coefficient eigensystem for the n lowest energy eigenvalues sorted into increasing order."; 


QuantumWireEnergies::usage="QuantumWireEnergies[{Lx,Ly},\[Alpha]][spikes][indices,n] returns the n lowest energy eigenvalues."; 


QuantumWireStates::usage="QuantumWireStates[{Lx,Ly},\[Alpha]][spikes][indices,n][{x,y}] returns the n lowest energy quantum states as functions of {x,y}, computed from the eigenvalues and vector basis functions."; 


QuantumWireStatesPlot::usage="QuantumWireStatesPlot[{Lx,Ly},\[Alpha]][spikes][indices,n] displays the n lowest energy quantum states, along with the corresponding energy."; 


(* ::Subsection:: *)
(*Begin Private context*)


Begin["`Private`"];


(* ::Section:: *)
(*Matrix Visualization *)


(* ::Text:: *)
(*Display matrices with the corresponding indices:*)


MatrixDisplay[m_List]:=
DisplayForm[
	StyleBox[
		GridBox[m,GridFrame->2,RowLines->{1,0},ColumnLines->{1,0}],
		Background->GrayLevel[0.8]
	]
]


LabeledMatrixDisplay[l_,mat_]:=MatrixDisplay[ArrayFlatten[{{{{""}},{l}},{Transpose[{l}],mat}}]]


(* ::Text:: *)
(*Show the rectangular region, along with the location of the \[Delta]-function spikes*)


PotentialDiagram[{Lx_,Ly_}][spikes:{x_,y_}]:=
Graphics[{
	Opacity[0.2],Rectangle[{0,0},{Lx,Ly}],
	Opacity[1],Red,AbsolutePointSize[5],Point@Transpose@spikes
	},
	Axes->True
]


(* ::Section:: *)
(*System Matrices*)


(* ::Subsection:: *)
(*Indices*)


(* ::Text:: *)
(*Use Tuples to construct the set of pairs of indices:*)


Indices[n_]:=Tuples[Range[n],2]


(* ::Text:: *)
(*To improve convergence, the optimal range of indices will depend upon \[Alpha] and the states of interest, so one can select a different set of indices, which should be optimised.*)


(* ::Subsection:: *)
(*Potential Matrix*)


(* ::Text:: *)
(*Entries in the potential matrix involve summation over a set of \[Delta]-function spikes located at points in the rectangle.*)


(* ::Text:: *)
(*Define each term in the summation:*)


\[CapitalDelta]Term[{Lx_,Ly_}][{x_,y_}][{n_,m_}]:=Sin[(n \[Pi] x)/Lx]Sin[(m \[Pi] y)/Ly]


(* ::Text:: *)
(*Construct the vector of indices for all the points using Listability of the terms in the summation:*)


\[CapitalDelta]Vector[{Lx_,Ly_}][spikes_][indices_]:=\[CapitalDelta]Term[{Lx,Ly}][spikes]/@indices


(* ::Text:: *)
(*Compute the potential matrix \[CapitalDelta] and the sum over the \[Delta]-function spike locations using an Outer product involving Dot products of vectors:*)


PotentialMatrix[{Lx_,Ly_}][spikes_][indices_]:=
Module[{vec=\[CapitalDelta]Vector[{Lx,Ly}][spikes][indices]}, 
	Outer[Dot,vec,vec,1]
]


(* ::Subsection:: *)
(*Hamiltonian matrix*)


(* ::Text:: *)
(*Compute the total Hamiltonian matrix:*)


HamiltonianMatrix[{Lx_,Ly_},\[Alpha]_:1][spikes_][indices_]:=
	((4\[Alpha])/(Lx Ly))PotentialMatrix[{Lx,Ly}][spikes][indices]+
	DiagonalMatrix[indices^2 . {\[Pi]/Lx,\[Pi]/Ly}^2]


(* ::Section:: *)
(*Energies and Eigenfunctions*)


(* ::Subsection:: *)
(*Energy-ordered eigensystems*)


(* ::Text:: *)
(*Re-order the eigensystem in order of increasing energy:*)


EnergyOrdered[es_]:= Transpose[Sort[Transpose[Chop[es]]]]


(* ::Subsection:: *)
(*Quantum Wire Eigensystem*)


(* ::Text:: *)
(*Compute and store the coefficient eigensystem for the n lowest energy eigenvalues sorted into increasing order:*)


QuantumWireSystem[{Lx_,Ly_},\[Alpha]_:1][spikes_][indices_,n_]:= 
QuantumWireSystem[{Lx,Ly},\[Alpha]][spikes][indices,n]=
	EnergyOrdered[
		Eigensystem[HamiltonianMatrix[{Lx,Ly},\[Alpha]][spikes][indices],-n]
	]


(* ::Subsection:: *)
(*Quantum Wire Energies*)


(* ::Text:: *)
(*The  energy eigenvalues in increasing order are the first component of QuantumWireSystem:*)


QuantumWireEnergies[{Lx_,Ly_},\[Alpha]_:1][spikes_][indices_,n_]:=
	First@QuantumWireSystem[{Lx,Ly},\[Alpha]][spikes][indices,n]


(* ::Subsection:: *)
(*Quantum Wire Eigenfunctions*)


(* ::Text:: *)
(*The energy-ordered eigenfunctions are computed from the last component of QuantumWireSystem via a dot product:*)


QuantumWireStates[{Lx_,Ly_},\[Alpha]_:1][spikes_][indices_,n_][{x_,y_}]:=
	 Last@QuantumWireSystem[{Lx,Ly},\[Alpha]][spikes][indices,n] . \[CapitalDelta]Vector[{Lx,Ly}][{x,y}][indices]


(* ::Subsection:: *)
(*State Visualization*)


(* ::Text:: *)
(*Display the eigenfunctions\[LongDash]computed from the eigenvalues and vector basis functions\[LongDash]along with the corresponding energy:*)


QuantumWireStatesPlot[{Lx_,Ly_},\[Alpha]_:1][spikes_][indices_,n_]:=
Module[{\[ScriptCapitalE],\[Psi],x,y},
	\[ScriptCapitalE]=QuantumWireEnergies[{Lx,Ly},\[Alpha]][spikes][indices,n];
	\[Psi]=QuantumWireStates[{Lx,Ly},\[Alpha]][spikes][indices,n][{x,y}];
	TabView@Table[\[ScriptCapitalE][[i]]->Plot3D[\[Psi][[i]],{x,0,Lx},{y,0,Ly},
		AxesLabel->{"x","y"},
		Mesh->None,
		PlotRange->2,
		BoxRatios->{Lx,Ly,1.5}],{i,Length[\[ScriptCapitalE]]}
	]
]


(* ::Section:: *)
(*Package Ending*)


End[]; 


(* ::Subsection:: *)
(*Setting Attributes*)


(*SetAttributes[{...}, {Listable, Protected, ReadProtected, NumericFunction}]; *)


(* ::Subsection:: *)
(*End package*)


EndPackage[];
