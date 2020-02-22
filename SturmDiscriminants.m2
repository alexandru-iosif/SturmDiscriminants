-- -*- coding: utf-8 -*-
newPackage(
	"SturmDiscriminants",
	Version => "0.1",
	Date => "October 2018",
	Authors => {{
		  Name => "Alexandru Iosif",
		  Email => "alexandru.iosif@ovgu.de",
		  HomePage => "https://alexandru-iosif.github.io"}},
    	Headline => "Computation of Sturm Discriminants",
	AuxiliaryFiles => false,    
    	PackageImports => {"Elimination"},
	DebuggingMode => false
)

export {
     -- 'Official' functions
     "SturmDiscriminant",
     "SturmSequence"

     -- Not in the interface:
--   "numeratorMatrix"
--   "denominatorMatrix"
--   "radicalMatrix"
--   "factorsMatrix"
}


SturmSequence = f -> (
-- by Dan Grayson & Frank Sottile:    
     assert( isPolynomialRing ring f );
     assert( numgens ring f === 1 );
     R := ring f;
     assert( char R == 0 );    
     x := R_0;
     n := first degree f;
     c := new MutableList from toList (0 .. n);
     if n >= 0 then (
     	  c#0 = f;
          if n >= 1 then (
               c#1 = diff(x,f);
               scan(2 .. n, i -> c#i = - c#(i-2) % c#(i-1));
               ));
       toList c)

-- the following function takes as input a matrix M and returns
-- another matrix whose entries are the numerator of M:
numeratorMatrix = M ->(
    matrix apply(entries M, i -> apply ( i, j -> numerator j) )
    )

-- the following function takes as input a matrix M and returns
-- another matrix whose entries are the denominator of M:
denominatorMatrix = M ->(
    matrix apply(entries M, i -> apply ( i, j -> denominator j) )
    )

-- the following function takes as input a matrix M and returns
-- another matrix whose entries are the the radicals of the entries of
-- M:
radicalMatrix = M ->(
    R := ring M;
    matrix apply(entries M, i -> apply ( i, j ->(flatten entries gens radical ideal j |{1_R})_0))
    )

-- the following function factors the entries of a matrix:
factorsMatrix = M ->(
    R := ring M;
    M = radicalMatrix M;
    factors := toList set flatten apply (entries M, i -> toList set flatten apply ( i, j -> ((toList factor j)|{1_R})) );
    apply( factors , i -> value i )
    )

-- SturmDiscriminant; the ideal I should be an ideal of
-- Rcoef[variables], where Rcoef=ring of coefficients
SturmDiscriminant = I -> (
    --add an option here; it seems that it takes a lot of time to compute the reduced discriminant; add an option Reduced = true or false
    R := ring I;
    Rcoef := coefficientRing R;
    assert( dim sub(I,frac(Rcoef)[flatten entries vars R]) == 0 );
    assert( char R == 0 );
    assert( not (class Rcoef) === FractionField );
    K := frac(Rcoef);
    Rflat := (flattenRing R)#0;
    J := symbol J;
    fgenerator := symbol fgenerator;
    eliminationVariables := symbol eliminationVariables;
    fsturm := symbol fsturm;
    LeadCoefficientsSturm := symbol LeadCoefficientsSturm;
    ConstantTermsSturm := symbol ConstantTermsSturm;
    LCCTMatrix := symbol LCCTMatrix;
    sturmdiscriminant := set {};
    UnivariateRing := symbol UnivariateRing;
    for i from 0 to numgens R - 1 do(
    	eliminationVariables = toList set flatten entries  vars R - set{(flatten entries vars R)_i};
    	eliminationVariables = flatten entries sub(matrix{eliminationVariables},Rflat);
    	J_i = eliminate(sub(I,Rflat),eliminationVariables);
    	assert (numgens J_i == 1 and not J_i == 0);
	UnivariateRing = K[(vars R)_i_0];
    	fgenerator = sub((gens J_i)_0_0,UnivariateRing);
    	fsturm = SturmSequence fgenerator;
	LeadCoefficientsSturm = sub(matrix {apply (fsturm, i -> leadCoefficient i)},K);
        ConstantTermsSturm = sub(matrix {apply (fsturm, i -> coefficient(1_UnivariateRing,i))},K);
	LCCTMatrix = LeadCoefficientsSturm | ConstantTermsSturm;
--    	sturmdiscriminant = sturmdiscriminant + set factorsMatrix sub(numeratorMatrix LCCTMatrix,Rcoef) + set factorsMatrix sub(denominatorMatrix LCCTMatrix,Rcoef) ;-- Uncomment this line if you want this function to compute a reduced discriminant
   	sturmdiscriminant = sturmdiscriminant + set flatten entries sub(numeratorMatrix LCCTMatrix,Rcoef) + set flatten entries sub(denominatorMatrix LCCTMatrix,Rcoef) ;-- Uncomment this line if you want this function to compute a non reduced discriminant
    	);
--    sturmdiscriminant = flatten entries sub (matrix{toList sturmdiscriminant},Rcoef);
    use R;
    toList (sturmdiscriminant - set{1_Rcoef} - set{0_Rcoef} - set{-1_Rcoef})
    )


beginDocumentation()

document {
    Key => SturmDiscriminants,
    Headline => "a package for computing Sturm Discriminants",
    
    EM "SturmDiscriminants", " is a package that makes use of
    Sturm sequences to compute discriminants of systems with
    positive roots",
    }

document {
     Key => {SturmDiscriminant},
     Headline => "Sturm Discriminant",
     Usage => "SturmDiscriminant I",
     Inputs => {
          "I" => { "a zero dimensional ideal of a polynomial ring of
          the form k[parameters][variables], where are the rational
          numbers "} },
     Outputs => {
          {"a list of polynomials in the parameters of I whose zero
          locus divide the positive orthant into cells with constant
          number of positive roots"} },
     EXAMPLE {
          "R = QQ[a,b,c,d,e,f][x,y]",
          "I = ideal (a*x+b*y+c,d*x+e*y+f)",
          "SturmDiscriminant I",
          },
     Caveat => {"This routines does not work if, after eliminiating
     all the variables except one, the resulting ideal is
     non-principal or zero."}}

end

restart
installPackage "SturmDiscriminants"


R = QQ[a1,a2,a3,c1,c2,c3][x1,x2,x5];
SturmDiscriminant ideal(x1+a1*x1*x2-c1, x5+a3*x1*x2-c2,
x2*x5+a1*x1*x2*x5+a2*x1*x2+a3*x1*x2*x5-c3*x5)

monomials (toList o3)_0

toString (toList o3)_0

monomials (toList o3)_1
monomials (toList o3)_2
monomials (toList o3)_3
monomials (toList o3)_4
monomials (toList o3)_5
monomials (toList o3)_6
monomials (toList o3)_7
monomials (toList o3)_8
monomials (toList o3)_9
monomials (toList o3)_10
monomials (toList o3)_11
monomials (toList o3)_12
monomials (toList o3)_13
monomials (toList o3)_14
monomials (toList o3)_15
monomials (toList o3)_16



For 2-sites phosphorylation
R=QQ[b1,b2,b3,b4,b5,b6,b7,s1,s2,s3]
I=ideal(b1*s1^2*s2*s3+b2*s1*s2*s3+b3*s1*s2-b3*c1,b4*s1^2*s2*s3+b5*s1*s2*s3+b3*s2-b3*c2, (b4+b1)*s1^2*s2*s3+b6*s1^2*s3+(b2+b5)*s1*s2*s3+b3*s3+b7*s1*s3-b3*c3)

In Mathematica:
E3=GroebnerBasis[{b1*s1^2*s2*s3+b2*s1*s2*s3+b3*s1*s2-b3*c1,b4*s1^2*s2*s3+b5*s1*s2*s3+b3*s2-b3*c2, (b4+b1)*s1^2*s2*s3+b6*s1^2*s3+(b2+b5)*s1*s2*s3+b3*s3+b7*s1*s3-b3*c3},{b1,b2,b3,b4,b5,b6,b7,s1,s2,s3},{s1,s2}];

R=QQ[c,]
I=ideal()
