#import "@preview/physica:0.9.5":*


// PARAMETERS
#set page(margin: (x:25mm, y:auto), numbering : "1")
#set text(font:"Lucida Grande", size:11pt)
#set heading(numbering : "1.")
#set par(justify: true)
#set math.mat(delim: "[")
#set figure.caption(separator: [ --- ])


#show emph: it => {
  text(10pt,font:"Helvetica", style:"italic", it.body)
}

#show link: set text(fill: blue)

#show figure.caption: c => context [
  #text(weight: "bold")[#c.supplement #c.counter.display(c.numbering) #c.separator]	
  #c.body
]

// FILE DESC
#set document(
     title:"Lab project report",
     author:"Erwan Le Doeuff (weld)",
     date: datetime(year:2025, month:05, day:27)
)


#let appendix(body) = {
  set heading(numbering: "A.1.", supplement: [Appendix])
  counter(heading).update(0)
  body
}


#show raw.where(block: false): box.with(
  fill: rgb(0%, 30%, 100%, 10%),
  inset: (x: 3pt, y: 0pt),
  outset: (y: 4pt),
  radius: 3pt,
)


#show raw.where(block: true): box.with(
  width: 100%,
  fill: rgb(0%, 30%, 100%, 20%),
  inset: (x: 3pt, y: 0pt),
  outset: (y: 6pt),
  radius: 3pt,
)


#text(17pt)[*Abstract*]

Sol–gel transitions are fascinating phenomena at the crossroads of chemistry, physics, and materials science. Found in everyday systems — from biofilms and mucus to silica aerogels, hydrogels in medicine, and food emulsions — gels exhibit a unique combination of structure and softness. Their tunable mechanical, thermal, and functional properties make them central to cutting-edge applications in biotechnology, soft robotics, and drug delivery. This report explores the sol–gel transition through two complementary lenses: the mean-field Flory–Stockmayer theory and a numerical implementation of Diffusion-Limited Cluster Aggregation (DLCA). Key observables from the literature are reproduced to validate our DLCA implementation, and the strengths and limitations of both approaches are examined. While Flory–Stockmayer provides an elegant, foundational framework for understanding gelation thresholds, it fails to capture key features such as spatial correlations or kinetic constraints. DLCA, although more computationally involved, offers a richer, dynamic description that opens paths toward realistic modeling of gelation in complex systems. We close with a personal perspective on the challenges and opportunities in modeling polymer gels.
\
\
\
\
\
#outline()

#pagebreak()

= Introduction

Duminil-Copin was awarded the Fields Medal in 2022 for his work on percolation theory. Quoting G. R. Grimmet#footnote[G. R. Grimmet, "The Work of Hugo Duminil Copin", 2022, #link("https://arxiv.org/pdf/2207.02022")[arXiv:2207.02022v1]] word for word: "Duminil-Copin has been instrumental in extending many of the classical percolation techniques from Bernoulli percolation (independent) models to a much wider class of dependent models." Whatever that means. Before it became a thing for mathematicians, percolation theory was born out of practical needs - originally among coal industry engineers trying to understand fluid flow in porous media, and almost simultaneously among chemists like Flory and Stockmayer, who used similar ideas to describe the growth of polymer networks at the microscopic scale.

In essence, percolation refers to the emergence of large-scale connectivity from purely local interactions. Imagine occupying sites or bonds randomly on a lattice: as the fraction of occupied elements increases, isolated clusters begin to form, and beyond a critical threshold, a giant connected component suddenly appears. This marks the percolation threshold, a hallmark of phase transition in disordered systems. Whether in lattices, networks, or physical systems, the central question remains the same: at what point does local connectivity become global?

One of the most iconic examples of such a connectivity-driven process is the sol–gel transition. In this transition, a system evolves from a collection of finite, disconnected clusters to the emergence of a spanning, system-wide network --- the gel. Understanding how this transition occurs, and predicting the onset of gelation, is not only of theoretical interest, but also critical for the design of materials in fields ranging from biopolymers to aerogels.

The theoretical framework proposed by Flory and Stockmayer in the 1940s treats the aggregation of monomers as a branching process and predicts quantities such as the gel fraction and the size distribution of finite clusters. However, it relies on strong simplifications: no closed loops, no spatial embedding, and randomness. Whether these assumptions hold in real systems is an open question.

In this report, we test the validity of the Flory–Stockmayer theory by comparing its predictions with numerical results obtained from a diffusion-limited cluster aggregation (DLCA) model. Unlike the theoretical model, DLCA explicitly accounts for spatial constraints and irreversible bonding due to Brownian motion. DLCA puts the theory to the test: it reveals what happens when gelation is governed not by abstract connectivity, but by actual particle motion and contact --- in silico.

#pagebreak()

= Theoretical approach

== The Big Picture

// Insérer les papiers de Flory et de Stockmayer
// Flory, P.J. (1941). "Molecular Size Distribution in Three Dimensional Polymers I. Gelation". J. Am. Chem. Soc. 63, 3083
// Stockmayer, Walter H.(1944). "Theory of Molecular Size Distribution and Gel Formation in Branched Polymers II. General Cross Linking". Journal of Chemical Physics. 12,4, 125

As mentioned in the introduction, Flory-Stockmayer theory was first formulated by P. Flory in the 1940s to explain the gelation process of hyperbranched polymers. Polymers are often described as long chains composed of repeating chemical units called monomers. It's a bit misleading, as it overlooks the structural diversity of polymer systems. Not all polymers form chains --- as in the case of dendrimers --- nor are they always long --- as in the case of oligomers. A more accurate definition would describe polymers as macromolecules composed of repeating units, which can organize into linear, branched, network-like architectures, ... A first broad classification might distinguish synthetic polymers --- typically made from one or two monomers types --- from naturally occurring polymers such as RNA, which exhibits greater chemical heterogeneity and biological functionalities (see *@Table1*). The number of monomers in a given molecule is usually denoted by $N$, and the corresponding molecule is often referred to as an $N$-mer.

#figure(
table(
	columns:(30%, 30%, 30%),
	inset: 6pt,
	align: top + left,
	table.header(
		[],[*Synthetic polymer*],[*Biopolymer*]
	),

	[Monomer diversity],
	[
		One or two\ monomer types\ (e.g. ethylene, styrene)
	],
	[
		Dozens of chemically distinct monomers\ (e.g. ribonucleotides)
	],

	[Sequence specificity],
	[Mostly random, \ sequence control \ possible],
	[Precisely encoded sequence],

	[Chain topology],
	[Mostly linear\ or branched],
	[Hierarchical,\ P1 $->$ P2 $->$ P3#footnote[Informal notation used to denote primary (P1), secondary (P2), and tertiary (P3) structure levels.]],

	[Structure/Function\ relation],
	[Weak,\ emergent from bulk properties],
	[Strong,\ 3D structure governs biochemical functions],

	[Self-assembly],
	[Rare,\ externally driven],
	[Intrinsic\ (e.g. protein folding)],

	[Typical functions],
	[Mechanical, structural,\  packaging],
	[Catalysis, recognition,\ regulation, transport]
),
	caption:[*Short comparison between synthetic and biological polymers*],
	placement:bottom,

)<Table1>

The topology of a synthetic polymer is inherited from the properties of its constituent monomers. In particular, a monomer can be characterized by its functionality $f$, _i.e._ the number of reactive sites it bears. In short, this determines how many other monomers a given monomer can bond with. The bonding of two monomers is called a polymerization process. This infers the overall structure of a given polymer. The case $f=2$ is ubiquitous in both the popular representation and scientific literature on polymers, as it leads to the formation of the "long chains" mentioned earlier. Such systems have been the subject of extensive theoretical and experimental studies, partly because they are relatively simple to model, but also because they exhibit remarkable physical properties. For instance, conductive polymers --- capable of transporting electric current --- only exist in linear configurations, as this is the only topology compatible with a favorable electronic delocalization.


Once you start considering monomers with increased functionality#footnote[In practice, functionalities greater than f = 3 are relatively rare in typical polymer systems, and referring to such high values as “monomer functionalities” may be something of a misnomer — it often reflects collective branching behavior rather than the chemistry of individual monomer units.] ($f >= 3$), the architecture of the resulting polymer becomes more complex, allowing for branching and eventually the formation of extended networks. In a given sample, these structures are generally not unique. Starting from a pool of monomers in solution, several independent growth fronts can develop in parallel, giving rise to distinct branched structures. Interactions between these structures can vary: the arms (in analogy with starfishes) may entangle, stack or remain isolated. However, once bonds start to form between distinct structures, a macroscopic connected component may emerge. This sudden shift --- from a polymixture of branched clusters to a single spanning network --- marks the gelation transition, _i.e._ the onset of a gel.

#figure(
	image("rsrc/fig-structure-polymers.png", width:80%),
	caption: [*Few polymer architectures*\ Adapted from *Rubinstein & Colby, 2003, p. 6 @rubinstein_polymer_2003*],
)<Figure1>


Although preliminary attempts had been made to describe this process, it was *Flory @flory_molecular_1941* who first provided a microscopically grounded and quantitatively predictive framework for gelation. His theory laid the foundation for our modern understanding of gel formation, even though it was not yet cast in the formal language of percolation theory --- a branch of statistical physics concerned with phase transitions driven by local connectivity rules. In order to make the problem mathematically tractable, Flory had to introduce several simplifying assumptions (to be discussed later), but the core idea remains physically meaningful.

The model assigns a probability (governed by physical parameters such as reaction rates, temperature, and monomer concentration) to the formation of a bond between two monomers. The central question then becomes whether a critical probability exists beyond which a gel appears. One of the most striking results is the prediction of a critical bonding probability,

$
	p_c = 1/(f-1), space f >= 3
$
The initial framework did not account for the possibility of mixing different chemical species in the solution, or for the use of crosslinking agents capable of artificially increasing the effective functionality of monomers — thereby enabling the use of linear polymers as building blocks. This limitation was later addressed by *Stockmayer @stockmayer_theory_1944*, whose work extended Flory’s approach and led to what is now known as the Flory–Stockmayer theory.



== Flory-Stockmayer model & mean-field theory

We will follow the development & notation presented in *Rubinstein et al. @rubinstein_polymer_2003* to derive the Flory–Stockmayer formalism. The situation is the following: we consider a reservoir of monomers in solution, and discard the need to define their functionality for the moment. Two monomers can polymerize with an independent probability $p$, commonly referred to as the extent of reaction in the literature. As $p$ increases, the number of chemical bonds in the system grows, leading to the progressive formation of larger and more complex molecular structures.

Initially, monomers link up into small clusters. These clusters grow as more bonds form, and some of them may merge. We can define a quantity $n(p,N)$ as the density of $N$-mers in the system at extent of reaction $p$. We postulate the existence of a critical threshold $p_c$, known as the gel point, beyond which a system-spanning connect component emerges. This marks the gelation transition, and the macroscopic connected structure is referred to as the gel.

At any moment, one can evaluate the total fraction of monomers that are not part of the gel. We denote this quantity by $P_"sol" (p)$ and call it the sol fraction. It accounts for all monomers that are either isolated or embedded within finite polymer structures. Mathematically, $P_"sol"$ can be written as a sum over the $N$-mers that are not in the gel,

$
P_"sol" (p) = sum_(N=1)^(+oo) N n(p,N)
$

A conjugate quantity of $P_"sol"$ can also be defined: the gel fraction, denoted $P_"gel"$, which corresponds to the fraction of monomers that belong to the infinite, percolating component --- the gel. By construction, the two quantities satisfy

$
P_"sol" (p) + P_"gel" (p) = 1
$

At this point, with our current definitions, we are not yet able to derive an explicit expression for the gel point $p_c$. However, we can already anticipate that obtaining an equation for the gel fraction $P_"gel"$ would be highly informative: since it is expected to remain zero below $p_c$ and grow continuously above it, $P_"gel"$ serves as an order parameter for this phase transition.

This is where the assumptions introduced by Flory come into play. To avoid interrupting the main derivation, we will postpone the discussion of their revelance until further sections. For now, we simply state them,

- Monomer bonding is random and occurs with equal probability $p$ between any pair of compatible sites
- All monomers are of the same type, i.e. have the same functionality $f$
- Cyclic structures (loops), whether within finite clusters or across the system, are neglected.

These assumptions are typical of a mean-field approach#footnote[Which reveals to be exact in dimensions greater than or equal to 6, where spatial fluctuations are negligible.], in which microscopic details are replaced by a simplified, analytically tractable picture. Concretely, this corresponds to replacing the monomer reservoir with an abstract structure known as a Bethe lattice — an infinite loopless graph (a tree in the sense of graph theory), where each node represents a monomer and each edge a chemical bond. Each node has a fixed degree equal to the monomer functionality $f$. This topology eliminates cycles and renders the derivation of both the gel point $p_c$ and the gel fraction $P_"gel"$ mathematically tractable through self-consistency arguments.

From the point of view of statistics, the Bethe lattice is translationally invariant. We can therefore begin our computations by focusing on a randomly chosen monomer, which we will call the root of the tree. Due to the absence of loops, the structure that unfolds from this root forms a hierarchy of generations: the monomer is bonded to a maximum of $f - 1$ neighbors in the first generation (excluding the parent link), each of which connects to a maximum of $f - 1$ new monomers in the next generation, and so on. Therefore, per generation, a monomer has in average $mu = p(f-1)$ children and $mu^2$ children plus grandchildren. So, after $n$ generations, the expected number of nodes is $mu^n$. Two cases then arise,

- if $mu < 1$, the expected number of nodes decreases with each generation. The tree dies out rapidly, and only finite clusters of bonded monomers can form.

- if $mu > 1$, the number of descendants grows exponentially. Some lineages persist indefinitely, resulting in an infinite, system-spanning cluster: the gel.

Remark that $mu > 1$ never happens for $f = 2$. The gelation transition occurs precisely when $mu = 1$, i.e. when the average number of descendants per monomer is exactly one. This yields the gel point,

$
p_c = 1/(f-1)
$

This critical value separates the regime of finite clusters from the onset of macroscopic connectivity. What about $P_"gel"$ now? To proceed, we introduce a new quantity $u$ which represents the probability that, when we follow a bond from a given monomer, we do not end up reaching the gel. In other words, $u$ is the probability that a bond leads only to a finite cluster. Let's think recursively and suppose we are at a given monomer, and we want to follow one of its bond.

- With probability $1-p$, no bond is formed at all --- so we certainly do not reach the gel.
- With probability $p$, a bond is formed --- but to avoid reaching the gel, all the other branches from the connected monomer must also lead only to finite clusters. Since these branches behave independently, the total probability is $u^(f-1)$.

Putting this together, we obtain a self-consistent equation for $u$,

$
u = 1 - p + p u^(f-1)
$

This equation tells us that a bond avoids the gel either because it does not exist (first term), or because it leads to a monomer whose remaining branches all avoid the gel (second term). Now, let’s go back to the original monomer --- the one we picked at random at the start. It has $f$ possible neighbors, meaning it can potentially form up to $f$ bonds. To determine whether this monomer belongs to the gel, we ask the following: do any of its bonds eventually connect it to the infinite cluster?

If none of the f bonds lead to the gel --- meaning all of them either don’t exist or only lead to finite clusters --- then this monomer is not part of the gel. Since all bonds behave independently, and each avoids the gel with probability $u$, the total probability that none of the $f$ bonds connect the monomer to the gel is $u^f$. But this is precisely the probability that the monomer itself is not part of the gel, i.e. the sol fraction. In other words,


$
P_"sol" = u^f space <=> space P_"sol"^(1\/f) &= u
$

Plugging it into the self-consistent equation yields,

$
 P_"sol"^(1\/f) = 1 - p + p  P_"sol"^((f-1)\/f)
$

We easily check that $P_"sol" = 1$ is always a solution, which means no gel, for all $p$. However, when $p > p_c$, a second solution appears. To see this, we rewrite the self-consistent equation as the problem of finding the zeros of the function:

$
h(x) = x^(1\/f) - 1 + p - p x^((f-1)\/f)
$


Graphically (see *@Fig2 (left)*), this means identifying the intersection points between the two curves $x^(1\/f)$ and $1 - p + p x^((f-1)\/f)$. Below the gel point, they intersect only once at $x = 1$. However, beyond a critical value $p_c$, they intersect twice: once at $x = 1$, and once at some $x < 1$  ($<=> P_"sol" < 1$). Since, $P_"sol" + P_"gel" = 1$, then $P_"gel"$ is necessarily non-zero. An analytical formula can be derived for $f = 3$, namely,

$
P_"gel" (p) = cases(
 1 - ((1-p)/p)^3  "if" p > p_c,
 0 "otherwise"
 )
$

With $p_c = 1\/(f-1) = 0.5$. This leads to *@Fig2 (right)*, where both $P_"sol"$ and $P_"gel"$ are shown. We can clearly observe that before $p_c$, $P_"gel"$ remains flat at zero. Then, it rises continuously with $p$, until it eventually reaches 1, at which point the entire system becomes a gel.

#figure(
    grid(
        columns: (auto, auto),
        rows:    (auto, auto),
        gutter: 1em,
        [ #image("rsrc/fig-h.png",   width: 92%) ],
        [ #image("rsrc/fig-phase-transition.png", width: 90%) ],
    ),
    caption: [*(Left)* Graphical resolution of the self-consistent equation for $x = P_"sol"$ at different values of $p$. All curves satisfy $h(x) = 0$ at $x = 1$, which corresponds to the absence of gel. For $p > p_c$, a second solution appears with $x < 1$, indicating the formation of a nonzero gel fraction. *(Right)* Evolution of the sol and gel fractions ($P_"sol", P_"gel"$) with the extent of reaction $p$ and the functionality $f = 3$. The system remains entirely in the sol phase ($P_"gel" = 0$) below the threshold $p_c$, and transitions continuously into the gel phase for $p > p_c$, where $P_"gel"$ grows smoothly from 0 to 1.],
) <Fig2>


== Critical exponents

Let's continue the discussion and derive more expressions for the Flory-Stockmayer model. The question we will ask ourselves is the following: what happens for $P_"gel"$ when we approach the gel point? We recall that just before the gel point,  $P_"sol"$ was 1. Right after, $P_"sol" < 1$. We can adopt a perturbative approach to the self-consistent equation near $p_c$, by setting $P_"sol" = 1 - delta space$ with $(P_"gel" = delta << 1)$.
We begin by plugging the perturbed expression into the self-consistent equation,

$
	P_"sol"^(1\/f) &= 1 - p + p dot P_"sol"^((f-1)\/f) \
	=> (1 - delta)^(1\/f) &= 1 - p + p dot (1 - delta)^((f-1)\/f)
$

Since $delta << 1$, both sides of the equation can be Taylor expanded to second order. We write,

$
(1 - delta / f - (f-1)/(2f^2) delta^2) &= 1 - p + p dot (1 - (f-1)/f delta - (f-1)/(2f^2) delta^2)
$

We now simplify both sides and isolate terms. Rearranging things around,

$
delta / f + (f-1)/(2f^2) delta^2  &=  p dot(f-1)/f delta + p dot (f-1)/(2f^2) delta^2 \
<=> 1 + (f-1)/(2f) delta &= p (f-1) + p dot (f-1)/(2f) delta \
$

It's time now to introduce a perturbation in $p$. Since we are near the gel point, we can perturb $p$ around its critical value, as $p = p_c + epsilon$ with $epsilon << 1$. Given $1 - p dot (f-1) = - epsilon (f-1)$, it follows that,

$
1 + (f-1)/(2f) delta &= p (f-1) + p dot (f-1)/(2f) delta \
<=> (f-1)/(2f) delta &= epsilon (f-1) + p dot (f-1)/(2f) delta \
<=> delta (1 - p) (f-1)/(2f) &= epsilon (f-1)
$

But since $1 - p = 1 - p_c - epsilon = [(f-2) dot (f-1)^(-1)] - epsilon$, we have,

$
 delta ((f-2)/(f-1) - epsilon)  (f-1)/(2f)  = epsilon (f-1)
$

Truncating to first order, the term in $epsilon delta$ can be safely neglected, yielding,

$
delta approx (2f(f-1))/(f-2) space epsilon
$

This justifies the classical result for mean-field bond percolation,

$
delta prop epsilon space => space P_"gel" prop (p - p_c )^beta "with" beta = 1
$

Such scaling laws are characteristic of second-order phase transitions, and the exponent $beta$ is, indeed, a critical exponent. It describes how rapidly the order parameter (here, the gel fraction) grows near the transition. In the vinicity of the critical point, many physical observables begin to obey similar power laws, each governed by its own critical exponent. For instance, Rubinstein & Colby derived the Fisher exponent $tau = 5\/2$ for the density of $N$-mers, in the form of the power law,

$
n(p_c, N) prop N^(-tau)
$

A proper derivation#footnote[Also called $epsilon$-expansions.] of these exponents should be done under the framework of the renormalisation group, _i.e._ by analyzing how the system's parameters transform under successive coarse-graining operations --- a method for which I admittedly lack the mathematical background to explore in full details. The advised and interested reader may read *@aharony_universal_1980 @stauffer_introduction_2018*.

Still, it seems that the deeper lesson here is that systems may not depend as strongly on their microscopic details as one might expect, due to scale invariance. Instead, they can often be grouped into universality classes --- categories of models that, despite differing at the microscopic level, exhibit the same macroscopic critical behavior. Within each class, certain relationships between critical exponents emerge naturally, known as hyperscaling relations. One may asks then, _do DLCA and Flory-Stockmayer models belong to the same universality class_? In short, it is very difficult to answer that question and it will be raised in our further developments, once we will have introduced the DLCA properly.

= Simulation approach through DLCA

== Diffusion-Limited Cluster Aggregation (DLCA)

The Flory-Stockmayer model and DLCA share a few similarities. Both start by randomly distributing monomers within a system --- typically a lattice --- and both are fundamentally stochastic. That said, the differences between them are far more important than the similarities.

First, DLCA is a dynamical model. There are many ways to implement it, but we will follow the version proposed by *Gimel et al.* @gimel_transition_1995, and adopt their notation. The system is defined on a cubic lattice of size $L^3$. We randomly place $N$ monomers (each of mass $m$) on the sites and define the initial monomer density $phi.alt_0 = N\/L^3$ which will be of a great importance. Unlike in Flory-Stockmayer, monomers are not fixed in place: they are allowed to move randomly across the lattice, mimicking Brownian motion. As time goes on, this motion brings them into contact.

Secondly, bonding is strictly local: when two monomers (or two clusters) touch, they stick together irreversibly, forming a new cluster. This models systems where attractive forces dominate thermal motion, so that aggregation is permanent. The dynamics are therefore purely aggregative: the number of clusters should decreases as they merge through time.


#figure(
	image("rsrc/percolating.png", width:45%),
	caption: [*Percolating cluster (top-bottom) in the DLCA model*, in red#footnote[Generated using VMD Tachyon ray tracer module. Positions of the atoms extracted from our simulations, for $L = 100, phi.alt_0 = 0.1$.].],
	//placement:bottom,
)<Figure3>





This process is implemented via a kinetic Monte Carlo algorithm. At each step, a cluster is chosen uniformly at random#footnote[A monomer is considered as a cluster of mass $m = 1$.], and a direction among the six nearest neighbors is selected. We then attempt to move the entire cluster one step in that direction. If all the destination sites are unoccupied, the move is accepted with probability

$
P = m^(-alpha) "with" alpha = 0.55
$

where $m$ is the mass of the cluster, and alpha#footnote[The choice of $alpha = 0.55$ corresponds to a fractal dimension $d_f approx 1.8$, which seems to be characteristic of three-dimensional DLCA aggregates --- that is, clusters whose mass scales sublinearly with their radius, due to their irregular, branching structure.] is a dynamical exponent related to the fractal dimension $d_f$ via ($alpha = 1\/d_f$). The dependence in $1\/m$ makes the larger clusters less mobile, which reflects how diffusive motion slows down with increasing size in real systems.

The simulation stops when all monomers belong to a single cluster, or when one cluster spans two opposite sides of the simulation box as shown in *@Figure3*. In that case, we call this cluster the gel.

== Main algorithms

The DLCA implementation is written in Fortran and follows a structure typical of kinetic Monte Carlo simulations. A detailed walkthrough of the code is provided in the appendices (see *@Code*). Still, a few key aspects are worth highlighting here:

- The code is organized into separate modules, as parts of the main procedure must be rewritten depending on the physical quantities of interest. Custom data types (similar to algebraic product types) and, as an example,  are used to represents particles: each particle stores its mass, spatial position, a unique ID and the cluster ID it currently belongs to.

- Whenever possible, some parallelisation is achieved using OpenMP --- for instance, when computing a quantity that depends on different values of $phi.alt_0$.

- As described in the previous section, at each attempt, a cluster is selected and possibly moved. This is done by independently moving all particles belonging to the same cluster, simply by checking their cluster ID. However, when two clusters merge, an algorithm is required to propagate one of the cluster IDs to all particles of the resulting cluster. This is handled using a so-called flood fill algorithm.

Basically, the flood fill procedure works by scanning all particles in the system and, for each unvisited one, launching a breadth-first search (BFS) traversal on the 3D grid. In this context, connectivity is defined through the six nearest-neighbor directions ($plus.minus x, plus.minus y, plus.minus z$), meaning that two particles belong to the same cluster if they occupy adjacent sites along one of these directions.

The BFS starts from a "seed" particle and explores all particles connected to it by visiting neighbors layer by layer: first the direct neighbors, then the neighbors of those, and so on. To manage this process, the algorithm uses a queue structure. Each time a neighbor is discovered, it is marked as visited and added to the queue for further exploration. This guarantees that all particles belonging to the same connected component are explored exactly once and assigned to the same cluster ID. The process is repeated until all particles have been visited, ensuring that each cluster is consistently labeled.


== Quality of Code


Unit testing doesn’t really make sense for the type of algorithm used here. Instead, what is needed are practical tools and physical arguments.

The main goal of our debugging tools is to verify the basic steps of the simulation: how do particles move on the lattice? Do they stay aggregated after collision? Can we detect percolation correctly? To inspect our trajectories, we used `VMD` (Visual Molecular Dynamics). Although `VMD` is designed for molecular dynamics, it supports standard file formats commonly used in that field. So our first step was to write a small utility to export simulation trajectories for visualization and debugging. We settled on the `XYZ` file format, which is lightweight and convenient: it only requires specifying the atom type, its coordinates, and the total number of atoms in each frame. The trajectories are stored uncompressed, which can result in fairly large files — but fortunately, debugging doesn’t require large systems, so this remains manageable.

The choice of using `VMD` came naturally, as it offers built-in tools and plugins to monitor a wide range of quantities. For instance, it can compute the pair-correlation function between two atoms or export a trajectory as an animated `GIF`. Its ability to render and customize visual output made it possible to generate *@Figure3*.

#figure(
	image("rsrc/mass-time.png", width:50%),
	caption: [*The maximum mass tracked along time,* for a given system. It serves a debugging purpose, as an healthy-check for the reasoned physicist.],
	placement:bottom,
)<Figure4>


While visualizing trajectories is certainly useful during the algorithmic development phase, it is by no means the only debugging strategy available. Other sanity checks can and should be implemented. For instance, given the irreversible & aggregative nature of the process we simulate, and assuming the flood-fill routine is correctly applied after each trial, the maximum cluster mass should be a monotonically increasing function of time.

*@Figure4* illustrates this expectation. We track the largest cluster mass (normalized by the total number of monomers) as a function of time. The curve shows a staircase-like growth, with long plateaus and sudden jumps. These jumps correspond to merging events, while the plateaus reflect intervals where no aggregation occurs. The regularity of this pattern provides a minimal yet convincing check that the core mechanics of the simulation behave as intended.

Our last word regarding quality of code is to say that we can imagine situations where the simulation cannot converge to one of the stopping criterion mentioned earlier.

In particular, when the density $phi.alt_0$ is low, the number of Monte Carlo attempts needed in order to converge may not be finite. A crude example would be a case where a pretty large not-percolating cluster is formed and we still have one monomer remaining. The large cluster is almost immobile, and the remaining monomer should then explore the whole volume, stochastically, until it finally reaches the cluster. Another example would be the formation of two large, immobile, non-percolating clusters, isolated from each other. Similar scenarii can be found, up to imagination.

The article makes no mention of how to handle these edge cases, leaving the practical implementation somewhat underspecified. It is likely that the authors introduced a pragmatic cutoff --- for instance, halting the simulation after a given number of successful moves. We imposed a criterion such that if the maximum cluster size is always the same after a certain number $X$ of succesful Monte Carlo attempts, the simulation is terminated. The issue is that, depending on the value of $X$, some of our reported quantities might be influenced by this arbitrary choice. 

In *@Figure5*, the probability $P$ to obtain a gel (a percolating cluster) is measured as a function of $phi.alt_0$. This figure should be directly compared with Figure 2 of the original article. Both exhibit the characteristics $S$-shaped curves. Nonetheless, increasing the cutoff value $X$ systematically shifts the transition to the left --- a trend consistently observed across various values of $(L, phi.alt_0)$.


#figure(
	image("rsrc/max-aggreg.png", width:50%),
	caption: [*Probability of gel formation as a function of $phi.alt_0$, for $L=20$.* \
Each curve corresponds to a different cutoff value for the maximum number of accepted Monte Carlo steps without aggregation. Increasing this cutoff shifts the transition to lower densities.],
	placement:top,
)<Figure5>

As a practical rule of thumb, we take $X$ to be equal to the volume of the simulation box, which corresponds to the number of steps a monomer would typically need to explore the entire system via random motion.


= Physical discussion

== Results of DLCA

This section serves a dual purpose. Naturally, we aim to discuss the physics underlying our simulation results. But more importantly, we seek to reproduce the findings of *@gimel_transition_1995*, since the authors did not make their implementation available. Put differently, to what extent are their result implementation-dependent?

#figure(
	image("rsrc/P.png", width:50%),
	caption: [*Probability of gelation $P(phi.alt_0, L)$ as a function of initial density.* \  Results shown for different system sizes $L$, each averaged over 100 independent simulations.]
)<Figure6>


//#figure(
  //grid(columns: 2, row-gutter: 2mm, column-gutter: 1mm,
  //image("rsrc/P.png"), image("rsrc/tgel.png")),
  //caption: [*(Left)* Probability of gel formation as a function of $phi.alt_0$. *(Right)* Gel time.]
//)

We begin by revisiting one of the earliest figures discussed in the original article. *@Figure6* shows the probability of gel formation as a function of the initial monomer density $phi.alt_0$, for various system sizes $L$. For small systems (e.g., $L = 5$), the transition is broad and gradual. As $L$ increases, the curve becomes progressively steeper, suggesting the development of a sharp threshold in the thermodynamic limit.

However, this steepening is accompanied by a noticeable shift of the transition point towards lower densities. This observation must be carefully interpreted: an extrapolation to ($L -> +oo$) would not only sharpen the transition but also lower the apparent threshold density. This supports the idea, already suggested in the original paper, that gelation can, in principle, occur at arbitrarily low densities — provided the system is large enough. In that sense, the process may no longer qualify as a genuine second-order phase transition, as postulated by Flory-Stockmayer theory, but rather as a purely kinetic, system-size-dependent crossover. The qualitative agreement with the results of *Gimel et al.* is excellent. On the quantitative side, our smaller systems tend to gel slightly later than those reported in the original article.

One of the key advantages of kinetic Monte Carlo methods is that they introduce a notion of dynamics, and with it, the ability to define characteristic timescales. This idea is illustrated in *@Figure7*, where we show the distribution of gelation times for various system sizes. Here, the physical time is defined as $t = 1 \/ N_C$, with $N_C$ the number of clusters remaining in the system. The average gelation time, reported as $t_"gel"$, appears largely independent of the system size $L$, and agrees within a margin of $plus.minus 5$ units of time with the values reported by the original authors. While the distributions seem to peak more sharply as $L$ increases — as mentioned in the article — our dataset remains too limited to draw firm conclusions on that trend.


#figure(
	image("rsrc/tgel.png", width:70%),
	caption: [*Distribution of gelation times for different system sizes $L$.* \ Each histogram is built from 100 independent simulations at fixed density $phi.alt_0 = 10%$. The mean gelation time $t_"gel"$ is shown in the legend. Standard deviation is given as 1$sigma$.],
)<Figure7>


The final set of results we wish to discuss concerns the polydispersity index, denoted by $K$. Polydispersity is defined as the ratio between the second and first moments of the cluster mass distribution, both normalized. Physically, it represents how broad the distribution of cluster sizes is: a low value of $K$ indicates that most clusters have similar sizes, while a high value reflects a wider disparity, with a few very large clusters coexisting alongside many smaller ones. Polydispersity thus serves as a useful indicator of the system’s structural heterogeneity as it evolves toward gelation.

The evolution of the polydispersity *@Figure8* is to be put alongside figure 6 of Gimel. Both figures illustrate the same general behavior: at low concentrations, the polydispersity remains nearly constant, close to the theoretical value $K = 2$ (a broader discussion of this value is available in the following section), indicating a narrow mass distribution. For higher concentrations, $K$ increases dramatically at a certain point in time, signaling the formation of a macroscopic cluster and the onset of gelation.

Quantitatively, the agreement is also satisfactory. The critical growth of $K$ appears sharper at higher $phi.alt_0$ and the divergence occurs earlier. In both cases, the transition from a below, stable $K$ to a steep increase is well captured. The shape and timing of these curves are consistent with the published data, although our $K$ values tend to saturate slightly later.



#figure(
	image("rsrc/polydispersity.png", width:50%),
	caption: [*Time evolution of polydispersity index $K$ at different concentrations.* System size fixed at $L = 200$. Each curve corresponds to a single trajectory for the indicated initial density $phi.alt_0$.],
)<Figure8>


== Comparison, limits of the models

It’s time to turn ourselves  to some qualitative criticisms of the models. DLCA, as presented by *Gimel et al.*, offers rather limited tunability. In practice, the dynamics is governed by very few parameters: the lattice geometry, the exponent $alpha$, and the mass of the monomer units.

If the lattice geometry remains an open direction of exploration, the parameters $alpha$ and $m$ directly govern the dynamics. Recall that the acceptance probability is proportional to $m^(-alpha)$. Changing either of these values would therefore modify how mobility is distributed among clusters, potentially favoring the growth of large aggregates or, on the contrary, hindering it. However, tuning $alpha$ is essentially a phenomenological choice, and there is little guidance from physical arguments to justify a particular _a priori_ value — this raises the question of how predictive the model really is.

More can be — and has been — done by relaxing some of the model’s constraints. In fact, several studies following this foundational work took that route. For instance, *M. Rottereau’s PhD thesis @rottereau_agregation_nodate*, supervised by the original authors, devotes an entire section to comparing how the dynamics and resulting structures change once the lattice constraint is released. This allows for the emergence of more realistic geometries, where clusters can adopt non-cubic shapes and where diffusion is less constrained by artificial boundaries.

Other possible extensions involve enriching the chemistry of the system. For example, one could introduce several types of monomers — say, A, B, and C — with specific rules for reactivity: A can only bind with B, B with C, and so on. Such a model would break the core assumption of the Flory-Stockmayer theory that all reactive sites are equivalent and independent. Similarly, the irreversible nature of the aggregation process in DLCA contrasts with many real-world gels that are thermoreversible. Capturing such behaviors would require introducing a temperature dependence into the model. Additionally, the role of the solvent — its viscosity, interactions, and screening effects — is completely absent from both Flory-Stockmayer and DLCA, though it is known to play a crucial role in the gelation pathways of real systems. All of these extensions could pave the way for more complex percolation scenarios and help explore how gelation thresholds shift in the presence of selective interactions or hierarchical binding rules — phenomena that are particularly relevant in biological systems and modern soft-matter design.

But that’s not all. The theory developed in *@stockmayer_theory_1944* predicts that the polydispersity index should lie between 1.5 and 2.0, given the assumptions of the model. In contrast, what we observe in *@Figure8* shows a clear deviation: even without changing the initial density, the approximation $K approx 2$ only holds at early times or in very low-density conditions — the so-called flocculation regime. Beyond that, $K$ increases significantly until it diverges, a behavior that the Flory-Stockmayer theory is fundamentally unable to predict.

Several other assumptions of the Flory-Stockmayer model also deserve to be challenged. As discussed in the theoretical section, the model assumes that clusters grow as loopless trees --— a mean-field approximation that entirely neglects intramolecular interactions. However, recent studies such as *de Keer et al. @de_keer_going_2021* suggest that the formation of loops and other intramolecular structures plays a significant role in the gelation process. In particular, these interactions tend to shift the percolation threshold away from the classical prediction $p_c = 1\/(f - 1)$, highlighting the limitations of the tree-like assumption when it comes to describing real systems.

Another key assumption — that of equal reactivity — also breaks down in real systems. Flory-Stockmayer theory assumes that all reactive sites are statistically equivalent: each has the same probability of reacting, regardless of its position within a molecule or its chemical environment. But in practice, steric effects can dramatically influence the accessibility and reactivity of a site. These effects introduce correlations between reactions that are entirely absent from the mean-field picture and can profoundly alter the kinetics and topology of gel formation. Moreover, it is fundamentally unable to predict gelation in systems where connectivity emerges from purely topological features, such as polymer knots or entanglements. These more complex mechanisms fall completely outside the scope of the Flory-Stockmayer framework --- which, to the best of my knowledge, has remained the most commonly invoked model to describe sol-gel transitions.


#pagebreak()
= Conclusion

== General conclusion

Overall, the work undertaken by Flory and Stockmayer has proved to be foundational in the way the scientific community understands sol-gel processes#footnote[If you found the topic interesting, I highly recommend T. Vicsek, _Fractal Growth Phenomena (Second Edition)_, World Scientific, June 1992 --- with the preface of B. Mandelbrot. It's so lovely.]. It established a theoretical framework that made the analogy with percolation theory not only possible, but fruitful. Despite its simplifying assumptions — equal reactivity, absence of loops, and a mean-field approach — the theory introduced core concepts such as gelation thresholds and molecular weight distributions. Its grounding in theoretical physics has had a lasting impact, and its conceptual simplicity made it possible to connect with deeper insights from critical phenomena and the renormalization group framework.

Still, the simplifying assumptions of the Flory–Stockmayer model have proven insufficient to capture the full complexity and diversity of gelation processes. Real-world systems exhibit a wide range of behaviors that go far beyond the scope of mean-field theory. The chemical nature of gels can vary dramatically, and the diversity of systems that undergo gelation is astonishing — from synthetic polymers and aerogels to biological systems like DNA and protein networks. Their apparent similarities often mask profoundly different underlying mechanisms. This highlights the need for more flexible, dynamic models capable of incorporating spatial correlations, kinetic constraints, and specific chemical interactions — all of which play a crucial role in the real physics of gelation.

As such, Diffusion-Limited Cluster Aggregation (DLCA) can be viewed as a response to a specific subset of these processes — namely, those where kinetics and spatial constraints dominate the aggregation dynamics. While the precise range of physical systems it accurately describes remains, to the writer, not entirely clear at this stage, DLCA nonetheless provides a valuable framework for exploring irreversible aggregation and the emergence of large-scale connectivity from local interactions.

In the end, both Flory–Stockmayer and DLCA tell part of the story. Neither is complete, yet each sheds light on different aspects of gelation — one through equilibrium statistical reasoning, the other through dynamic irreversibility. As our understanding deepens and experimental systems become more diverse, the dialogue between theory and simulation will remain essential. There is still much to explore.

== A personal note

A few closing words to share my personal take on this work. The physics of polymers is vast, and navigating the many available modeling approaches has been a real challenge. Classical tools like molecular dynamics are often ill-suited for gels: these systems are inherently large, making simulations computationally expensive, and the choice of force field — that is, the interactions considered — can drastically alter the system’s behavior, sometimes leading to states with completely different physical properties.

I initially approached my supervisors with the hope of establishing a dialogue with the ongoing experimental work in their lab. In that respect, the objective wasn’t fully met — the gap between modeling and experimental relevance turned out to be more difficult to bridge than I had anticipated. Still, my desire to deepen my understanding of sol-gel processes and polymer physics has been fully satisfied. I first encountered soft matter in my second undergraduate year, and since then, I’ve been wholly dedicated to learning more about it. I do regret how little visibility the field sometimes receives, especially as attention seems to have shifted toward biological and active systems. I admire my supervisors — Isabelle, Florian, and Boris — for their contributions to soft matter, and I hope to follow in their footsteps, bringing my computational skills to bear and contributing my own small stone to the collective effort.

I have some ideas that could be of interest for future students willing to explore this field. One direction worth investigating lies in the theoretical treatment of the gel itself. In classical mean-field models like Flory–Stockmayer, the gel is never truly part of the ensemble — it’s excluded from sums, inaccessible to physical observables, and often treated more as a conceptual threshold than a real, evolving object. This leaves many key questions unanswered. What does it really mean, microscopically, for a system to gel? How can we develop observables that genuinely probe the gel phase, rather than just its onset? Despite their importance across chemistry, biology, and materials science, we lack a microscopic, dynamic framework that convincingly captures the emergence and structure of the gel itself. In that sense, I believe sol–gel transitions are still theoretically underdeveloped.

As for me, I intend to keep working on DLCA — I’m starting to see it as a useful stepping stone for generating controlled gel structures that could serve as initial conditions for coarse-grained molecular dynamics simulations. To the reader: thank you for your time, and I hope this work has been as enjoyable to read as it was for me to write and explore.

#pagebreak()

#bibliography("rsrc/biblio.bib")

#pagebreak()

#show: appendix


= Code overview<Code>


This appendix presents a structured overview of the Fortran code developed for the DLCA simulations. Rather than describing the implementation narratively, we adopt the format of a reference manual, similar to an Application Programming Interface (API) documentation. Each module, type, and key subroutine is introduced with its purpose, inputs and ouputs, and relevant implementation details. This structure aims to facilitate future reuse or extension of the codbase.

The source can be found on #link("https://github.com/atomicwelding/DLCA")[GitHub].



== Language, dependencies

The codebase is written in Fortran 2008 and makes use of modern features such as derived types and modular programming. It has been tested and compiled using GNU Fortran (`gfortran`), version 14.02.

The project relies on minimal external dependencies. It is mostly self-contained, with the exception of `OpenMP`, which is used to parallelize specific parts of the computation. However, this parallelism is only invoked in the provided main file, which serves as a user template. Depending on the physical quantity one wishes to measure, the user is expected to write or adapt their own main program. As such, `OpenMP` is optional and can be disabled if not needed.

A `build.sh` shell script is provided to compile the project. It may require minor adjustments depending on the operating system, especially to correctly link `OpenMP`. The code has been successfully built and tested on macOS Sequoia 15.15.

Compilation typically produces an executable named `dlca`. It does not take any option/argument.

== Parameters (`params.f90`)

This module defines all global constants and control parameters used throughout the DLCA simulation. It encapsulates system size, runtime settings and file output path in a single location.

```fortran
character(len=*), parameter :: filepath
```

- *filepath* --- Path where the main output will be written.



```fortran
integer, parameter :: L
```

- *L* --- Linear size of the cubic simulation box. The total number of lattice sites is $L^3$.



```fortran
real, parameter :: alpha
```

- *alpha* --- Dynamical exponent governing cluster mobility. Move probability is scaled as $P prop m^(-alpha)$, where $m$ is the cluster mass. Typically, $alpha = 0.55$.

```fortran
  integer, parameter :: Npts
  real, parameter :: MIN_PHI, MAX_PHI
  integer, parameter :: runs 
```

- *MIN_PHI*
- *MAX_PHI*
- *Npts* --- Number of density points to simulate between `MIN_PHI` and `MAX_PHI`.
- *runs* --- Number of runs per point.

```fortran
  integer, parameter :: MAX_STEPS_WITHOUT_AGGREGATION
```

- *MAX_STEPS_WITHOUT_AGGREGATION* --- Stopping criterion. Maximum number of consecutive Monte Carlo steps without the maximum cluster size growing. A good heuristic is $L^3$.

== Types (`types.f90`)

This module defines the fundamental data types used to represent monomers, clusters, and the overall simulation state. It provides the foundational structures on which the entire DLCA simulation is built.

```fortran
type :: Particle 
     integer :: id
     integer :: x, y, z
     integer :: m = 1
     integer :: cluster_id 
end type Particle
```
Represents a monomer in the system. Each particle knows its own position, mass, and which cluster it belongs to.
- *id* --- Unique identifier of the particle.
- *x,y,z* --- Coordinates of the particle on the 3D lattice.
- *m* --- Mass of the particle (default = 1).
- *cluster_id* --- ID of the cluster to which the particle currently belongs.

```fortran
type FaceTouch
     logical :: min_x = .false.
     logical :: max_x = .false.
     logical :: min_y = .false.
     logical :: max_y = .false.
     logical :: min_z = .false.
     logical :: max_z = .false.
end type FaceTouch
```
A utility structure used to check for percolation of a cluster through the simulation box. Each field is a boolean flag indicating whether a cluster touches the corresponding face of the cube.
- *min_x, max_x* --- Contact with the $x=0$ or $x=L-1$ faces.
- *min_y, max_y* --- Contact with the $y=0$ or $y=L-1$ faces.
- *min_z, max_z* --- Contact with the $z=0$ or $z=L-1$ faces.

```fortran
type :: SimulationState
     integer :: N 
     integer, dimension(0:L-1, 0:L-1, 0:L-1) :: grid 

     type(Particle), allocatable :: particles(:)
     integer,        allocatable :: cluster_size(:)

     integer :: no_growth_steps

     logical :: hasPercolated 
     logical :: hasEnded
end type SimulationState
```
Central structure storing the entire state of the simulation at a given moment. This includes the particles, grid occupancy, cluster sizes, and flags for termination or percolation.
- *N* --- Total number of particles in the simulation.
- *grid* --- Each lattice site contains either -1 (empty) or a particle `id`.
- *particles* --- Allocatable array of all particles in the system.
- *cluster_size* --- Size of each cluster (indexed by `cluster_id`).
- *no_growth_steps* --- Counter used to detect stagnation.
- *hasPercolated* --- `.TRUE.` if a percolating cluster has been detected.
- *hasEnded* --- `.TRUE.` if the simulation termination condition has been met.

== Subroutines (`subroutines.f90`)

This module implements the main simulation logic of the DLCA model. It provides a set of core routines to initialize the simulation state, evolve the system through Monte Carlo steps, and evaluate stopping conditions. These routines operate on the `SimulationState` type, defined in the `types.f90` module.

The module is designed to be modular and readable, with clearly separated responsibilities. Each subroutine can be reused independently and modified to fit alternative simulation strategies if needed.

```fortran
subroutine init_particles(SimulationState)
```
*Purpose*: Randomly distributes `N` monomers in the 3D lattice without overlap. Initializes particle IDs, positions, and cluster assignments.
- *Input*:  `SimulationState`, `intent(inout)`.
- *Output*: Populates the `grid`, `particles`, and `cluster_size` arrays.


```fortran
integer function count_active_clusters(SimulationState) result(res)
```
*Purpose*: Returns the number of active clusters (clusters with non-zero size).
- *Input*:  `SimulationState`, `intent(in)`.
- *Output*: `integer :: res` --- Number of distinct, non-empty clusters.


```fortran
subroutine flood_fill_all(SimulationState)
```
*Purpose*: Performs a flood-fill on the grid to update `cluster_id` for all particles. Ensures connected particles share the same cluster ID.
- *Input*:  `SimulationState`, `intent(inout)`.
- *Output*: Updates `SimulationState%particles(:)%cluster_id`.


```fortran
subroutine update_cluster_sizes(SimulationState)
```
*Purpose*:  Recomputes the size (number of particles) of each cluster.
- *Input*:  `SimulationState`, `intent(inout)`.
- *Output*:  Updates `SimulationState%cluster_size(:)` array.


```fortran
subroutine checkEndSimulation(SimulationState)
```
*Purpose*: Checks if the simulation has ended based on three criteria,
1. All particles belong to a single cluster.
2. A cluster spans opposite faces of the box (percolation).
3. A maximum number of MC steps has occurred with no growing events.
- *Input*:  `SimulationState`, `intent(inout)`.
- *Output*: Sets `SimulationState%hasEnded` and `SimulationState%hasPercolated` flags accordingly.

```fortran
subroutine trial(SimulationState)
```
*Purpose*: Performs a single Monte Carlo step.
1. Selects a random cluster.
2. Attempts a move in a random direction.
3. Applies acceptance criterion based on mass and `alpha`.
4. Updates state and checks stopping conditions.
- *Input*:  `SimulationState`, `intent(inout)`.
- *Output*: Updates particle positions, cluster IDs, grid, and termination flags.

== Main program (`main.f90`)

This part of the code is meant to be adapted or rewritten by the user, depending on the physical quantity one wishes to compute. Still, we can illustrate what can typically be done with it through a concrete example.


*Purpose* The provided version performs a systematic DLCA simulation campaign: for a range of initial monomer densities $phi.alt_0$, it runs several independent simulations and measures,

- the probability $P$ to obtain a gel (i.e., a percolating cluster),
- the average gelation time $t_"gel"$ (conditioned on percolation),
- the average gel fraction $P_"gel"$ (conditioned on percolation).

These quantities are averaged over runs independent realizations for each $phi.alt_0$.

*Notes*

- The loop over runs is parallelized with `OpenMP`.
- Random seeds are reinitialized at each run using `system_clock`, to reduce correlations.
- The simulation stops either when a percolating cluster is formed, all particles belong to a single cluster, or the maximum number of allowed steps without aggregation is reached.
- If percolation is impossible for a given $phi.alt_0$ (e.g., if $N < L$), the program skips the simulation and writes trivial outputs.

*Output*
Results are written to the file defined by the filepath parameter. The output file contains a short header followed by one line per $phi.alt_0$, following the format,

```
phi0   P   t_gel   P_gel
```

If no percolation occurred across the runs, we conventionally assign $t_"gel" = -1$. Once there, the data are plotted using `python` scripts.

