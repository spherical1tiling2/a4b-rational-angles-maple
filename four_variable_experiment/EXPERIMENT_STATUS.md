# Four-variable experiment for equation (2.6)

This is an experimental extension. It is not a complete classification or a theorem of the manuscript. The progress counts below describe the archived report; additional checkpoint files may reflect later partial runs.

## System and phase coordinates

The input consists of the pentagon angle sum, the vertex relation $\alpha+\beta+\gamma=2$, and equation (2.6), with unknown even $f\ge12$. All angles in the formulas are normalized by $\pi$. Eliminate $\alpha=2-\beta-\gamma$ and $\varepsilon=1+4/f-\delta$.

After half-angle substitution and support-lattice reduction, use $X=e^{i\pi\beta}$, $Y=e^{i\pi\gamma}$, $W=e^{2i\pi\delta}$, and $T=-e^{4i\pi/f}$. The resulting polynomial has 32 terms, multidegree $(3,3,3,2)$ in $(T,X,Y,W)$, is reported squarefree and irreducible over the rationals, and has full support-difference lattice $\mathbb Z^4$ of index one. Only nonzero constants and invertible Laurent monomials are cleared.

## First elimination layer

Eliminating $W$ with 15 nonidentity sign transformations and 16 square transformations produces 31 nonzero resultants with no common-factor branches. There are 86 factor occurrences and 42 distinct three-variable factors up to nonzero scalar multiples: eight rank-deficient projections, 15 full-rank factors of lattice index two, and 19 already full-lattice factors.

Before support-lattice reduction, intermediate factors can coincide with transformations and yield gcd branches and artificial degree growth beyond 200. Hermite support-lattice reduction removed those gcd branches in the processed full-lattice cases. The reduction is required at each layer.

## Back-substitution of deficient-rank projections

Under $0<\alpha,\beta,\gamma,\delta,\varepsilon<1$ and even $f\ge12$, the projections $Y=1$, $X=1$, $T=1$, $T=-1$, and $XY=1$ are excluded by angle boundaries or the range of $f$. The additional two-variable projection has only $X,Y=\pm1$ torsion candidates, also on the boundary. The substantive remaining projections are $X=Y$ and $TXY+1=0$.

For $X=Y$, exact factorization gives $F(T,X,X,W)=X(X-1)^2(TX^2+1)(T-W)^2$. The two resulting rational families are

\[
(a,1-a/2,1-a/2,1/2+2/f,1/2+2/f),\quad a\in\mathbb Q,\ 0<a<1,
\]

and

\[
(4/f,1-2/f,1-2/f,d,1+4/f-d),\quad d\in\mathbb Q,\ 4/f<d<1,
\]

for even $f\ge12$. They intersect at $a=4/f$, $d=1/2+2/f$. In the first family, both terms of equation (2.6) vanish because $\beta=\gamma$ and $\delta=\varepsilon$.

The projection $TXY+1=0$ is equivalent to $\alpha=4/f$ in the allowed angle range. Its remaining factor is $Q=W^2X^2Y^2+W(Y-X)(XY+1)-1$. The 15 three-to-two branches yielded 32 factor occurrences and 14 distinct support-lattice projection records, with no gcd branches. Further back-substitution is needed to decide whether this produces additional isolated solutions or families.

## Scope of the archived progress

At the time of the report, eight deficient-rank projections and 18 full-rank factors with at most 200 terms had been processed. Sixteen large factors of 582--2287 terms from the square transformations still required deeper elimination. The two-to-one torsion candidates, full back-substitution, and complete list of isolated solutions remained unfinished.

The two families above are confirmed solutions of the input equations. Their existence does not establish completeness or geometric realizability as tilings.

The implementation must distinguish no solution, isolated angle tuples, and positive-dimensional torsion cosets or parameter families. Each layer requires factorization, support-lattice reduction, gcd splitting, resultants, refactorization, rank detection, and back-substitution. A rank-deficient branch retains free parameters rather than eliminating them as if the solution set were zero-dimensional.

## Recorded formulas
\[
\alpha+\beta+\gamma+\delta+\varepsilon=3+\frac4f,
\qquad
\alpha+\beta+\gamma=2,
\]

\[
\begin{aligned}
&\Big[(1-\cos\pi\beta)\sin\pi\Big(\delta-\frac\alpha2\Big)
 -(1-\cos\pi\gamma)\sin\pi\Big(\varepsilon-\frac\alpha2\Big)\Big]
 \sin\frac{\pi(\delta-\varepsilon)}2\\
&\qquad -(1-\cos\pi(\beta-\gamma))
 \sin\frac{\pi\alpha}2\sin\frac{\pi(\delta+\varepsilon)}2=0.
\end{aligned}
\]

\[
\alpha=2-\beta-\gamma,
\qquad
\varepsilon=1+\frac4f-\delta.
\]

\[
X=e^{i\pi\beta},\qquad
Y=e^{i\pi\gamma},\qquad
W=e^{2i\pi\delta},\qquad
T=-e^{4i\pi/f}.
\]

\[
2^{4+1}-1=31
\]

\[
Y-1,\quad X-1,\quad T-1,\quad T+1,\quad X-Y,\quad XY-1,
\quad TXY+1,
\]

\[
F(T,X,X,W)=X(X-1)^2(TX^2+1)(T-W)^2.
\]

\[
TX^2+1=0
\quad\text{or}\quad
T-W=0.
\]

\[
\boxed{
(\alpha,\beta,\gamma,\delta,\varepsilon)
=
\left(
a,\ 1-\frac a2,\ 1-\frac a2,\
\frac12+\frac2f,\ \frac12+\frac2f
\right)
}
\]

\[
f\ge12\text{ is even},
\qquad a\in\mathbb Q,\quad 0<a<1.
\]

\[
\boxed{
(\alpha,\beta,\gamma,\delta,\varepsilon)
=
\left(
\frac4f,\ 1-\frac2f,\ 1-\frac2f,\
d,\ 1+\frac4f-d
\right)
}
\]

\[
f\ge12\text{ is even},
\qquad d\in\mathbb Q,\quad \frac4f<d<1.
\]

\[
a=\frac4f,
\qquad d=\frac12+\frac2f.
\]

\[
F\big|_{TXY=-1}
=-(XY-1)(X-Y)Q(X,Y,W),
\]

\[
Q=W^2X^2Y^2+W(Y-X)(XY+1)-1.
\]

\[
\texttt{none}\quad\text{or}\quad
(\alpha,\beta,\gamma,\delta,\varepsilon)
\]
