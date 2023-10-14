<figure id="fig-intro">
<img src="fig-intro.png" />
<figcaption>Illustration of <em>A15</em> at three different scales. <span id="fig-intro" label="fig-intro"></span></figcaption>
</figure>

# Introduction to Floats and Spatial Partitioning with *A15*

Binary floating-point numbers, commonly known as *floats*, are notorious for giving slightly different, *approximately correct* answers.

$$\dfrac{\mathit{\scriptstyle significand}}{2^{\mathit{exponent}}}$$

Akin to the expression above[^1], floats resemble fractions or ratios. Their integer numerators cycle linearly $0$–$\mathit{significand}$ once per denominator, whereas their log-linear denominators must double or split on strict powers-of-two. This representation approximates the vast majority of rational, base$_{10}$ numbers, and accumulating tiny, order-dependent rounding errors is *expected*. There are ninety-three approximations in the first one-hundred $^1/_n$ reciprocals alone, where denominators are *most* dense.

Floats shed the stability of integers to dramatically boost their reach. Hardware-accelerated for decades in pursuit of ever-more floating-point operations per second (FLOPS), a single, fixed-bit memory format meaningfully ranges from quantum foam to cosmic web. Floats are *inescapably* abundant—the sand of software—their skillful workings endow modern data structures with great strength and timeless clarity.

This research accepts and embraces binary floats as they are. It scales both the *A15* phase structure, also known as $\beta$–$W$, and its corresponding Voronoi honeycomb, the Weaire–Phelan honeycomb, to align with precise, IEEE 754-2008 floating-point specifications. Each honeycomb cell collapses its internal, 3D floating-point coordinates to an integer-encoded *A15* site at its core, and together, cores identify a compact, higher-order space—a well-rounded snapshot of the original, 3D floating-point space—where every high-dimensional *A15*-encoded coordinate mirrors an exact 3D float.

This research applies to IEEE 754-2008 base$_{2}$ floating-point numbers of all bit sizes, and may refer to them as *binary$_{64}$*, *binary$_{32}$*, *base$_{2}$ floats*, or simply *floats*. For coding and analysis purposes, binary$_{64}$ is preferred, due to its large size and prevalence in modern CPUs. For baselines and performance comparisons, binary$_{32}$ is preferred, due to its widespread presence in hardware, software, and network stacks.

The nomenclature and parlance used throughout is primarily that of crystallography, borrowing from other disciplines as necessary.

<figure id="fig-cell2">
<img src="fig-cell2.png" />
<figcaption>Left-handed <span class="math inline">\(^1/_2\)</span> unit cell. <span id="fig-cell2" label="fig-cell2"></span></figcaption>
</figure>

## Foundational Understanding of *A15* and Spatial Partitioning

Partitioning virtual spaces is a matter of fairness. “Fairness," as it applies to transformations on structures in 3D space, is a measure of isometry and isotropy—reducing the bit space *must not* significantly warp distances and angles between any two sites. While isometry on its own is readily achievable, combining it with isotropy is much more difficult. The $SO(3)$ group, or the set of all possible 3D rotations, is spherical—highly-isotropic structures appear “rounder" from the perspective of an individual site—and simply cannot fit nicely inside a cubical lattice “box". This innate tension between translation-preserving symmetries and rotation-preserving symmetries drastically shrinks the pool of N-fold designs available to perfectly isometric 3D lattices. In accordance with the crystallographic restriction theorem, C12 is the maximum coordination number, and the only possible angles are 180$^{\circ}$ (2-fold), 120$^{\circ}$ (3-fold), 90$^{\circ}$ (4-fold), and 60$^{\circ}$ (6-fold)—icosahedral designs (5-fold) with $I_h$ symmetry are not possible.

*A15* eschews perfect isometry. It trades 100% identical sites for a blended mix of exactly 75% C14 major sites—axes-aligned tetradecahedral layers (Weaire–Phelan) or cubes (Tetrastix) with 14 connections each—25% C12 minor sites—pyritohedral voids (Weaire–Phelan) or cubes (Tetrastix) with 12 connections each—and two different site-to-site distance metrics. True 5-fold symmetry appears in the form of alternating left- and right-handed sites with $T_h$ *pyritohedral* symmetry, an isometric subgroup (4-of-10 3-fold axes) of the full icosahedral symmetry group $I_h$. This localized asymmetry drastically increases isotropy (13.5 mean coordination) without impacting long-range isometric order.

<figure id="fig-wp-ts">
<p><img src="fig-wp.png" alt="image" /><img src="fig-ts.png" alt="image" /></p>
<figcaption>The Weaire–Phelan honeycomb (left) and the Tetrastix prism (right).<span id="fig-wp-ts" label="fig-wp-ts"></span></figcaption>
</figure>

In a sense, *A15* is already binary. Its fractional lattice coefficients use nothing but the first three multiples of $2^{-2}$, and all eight basis sites are perfect binary floats. Quadruple its fractional coordinates into the integers, and its two, site-to-site distance metrics become $2$ (major-major) and $\sqrt{5}$ (major-minor). $\sqrt{5}$ is the hypotenuse of a 2:1 right triangle and the crux of the golden ratio. *A15* can be defined at unit scale without stability issues. However, at any scale, *A15* only represents the destination encoding, leaving open the question of *how* higher-density bit spaces should discretize themselves to a valid *A15* site. In other words, finding the nearest site requires a precise definition of *nearest*.

The first definition, depicted in <a href="#fig-wp-ts" data-reference-type="autoref" data-reference="fig-wp-ts">[fig-wp-ts]</a> (left), is available to any 3D point set. Starting from an *A15* integer crystal lattice, identify its Voronoi honeycomb from the set of inflection points between neighboring *A15* sites—edges in this secondary structure have exactly two nearest neighbors in *A15*, and vertices have three or more—and the Weaire–Phelan honeycomb appears. A simpler, less isotropic definition of *nearest* is also available to *A15*. As seen in <a href="#fig-wp-ts" data-reference-type="autoref" data-reference="fig-wp-ts">[fig-wp-ts]</a> (right), when the angle between sites in the secondary structure is fixed to 90°—forcing exactly six nearest neighbors and filling the space with *unit* cubes—the Tetrastix prism emerges instead. The price for this simplicity is reduced spatial accuracy and more directional aliasing.

## Spatial Indexing and Related Spatial Partitioning Methods

Spatial partitioning is closely associated with spatial indexing. In this context, the partitioner is more dynamic and specific—space is split on maximally-coincident hyperplanes, or enclosed within minimally-overlapping polytopes, for the purposes of cataloging sites and facilitating retrieval. Incoming sites are unlikely to adhere to any meaningful symmetry and freely utilize the full range and precision of the ambient space. This ignorance of an implied external structure is critical for spatial indexing, but renders well-known binary space partitioners (BSP), like octrees and KD trees, and bounded-polytope solutions, like R-trees, R$^*$-trees, and its derivatives, less attractive as *implicit*, memory-efficient, interactive virtual space partitioners, because their preferred hyperplanes and polytopes do not maintain the spatial symmetries of the ambient space. However, if they did maintain ambient symmetries, the resultant structures might resemble objects that are comparable to *A15*: space groups and space-filling honeycombs.

## Comparative Analysis of *A15* with Space Groups and Honeycombs

Three-dimensional space admits fourteen crystallographic lattice types known as Bravais lattices—fourteen distinct, prototypical pairings between one-of-seven lattice systems and one-to-four lattice centerings—and every discrete, *periodic* tesselation of 3D space shares its translational isometries with a Bravais lattice. Non-translational isometries, such as reflections and rotoinversions, are known as 3D point groups, and the thirty-two that satisfy the crystallographic restriction theorem are deemed the crystallographic point groups. The complete set of 230 space groups emerges from all isomorphic combinations of the fourteen lattice types with the thirty-two crystallographic point groups, and fully characterizes any periodic tesselation of 3D space.

<figure>
<p><strong>(todo) Table.</strong> <span class="math inline">\(Pm\bar{3}n\)</span> (223) alongside other groups.</p>
</figure>

*A15*’s space group, $Pm\bar{3}n$, pairs the $O_h$ symmetry of the $cP$ Bravais lattice with the $T_h$ pyritohedral symmetry of the $m\bar{3}$ crystallographic point group. $T_h$ pyritohedral symmetry is an isometric subgroup of the *non-crystallographic*, full icosahedral symmetry group, $I_h$. Crystallographic point groups with $T_h$, $O$, and $T_d$ symmetries are all order 24, and second only to order 48, $O_h$ cubic symmetry. Since $T_h$ is the maximal subgroup between $O_h$ and $I_h$—between the existing cubical-octahedral isometries of *A15*’s Bravais lattice and the highly-desirable, *non-crystallographic* icosahedral isometries of $I_h$—any point group with higher order than $m\bar{3}$ is also more cubical. The link from $m\bar{3}$ to $I_h$ symmetry through $T_h$ symmetry is strong evidence that $m\bar{3}$ is isotropically ideal.

<figure>
<p><strong>(todo) Table.</strong> <em>A15</em> alongside other honeycombs (mean coordination, etc).</p>
</figure>

(todo) Tetrahedral-octahedral honeycomb

:   Simplectic, quasiregular honeycomb and face-centered cubic (FCC or $A_3$ or $D_3$) lattice; mean coordination 12; space group $Fm\bar{3}m$ (225); vertex- and edge-transitive; ideal 3-space packing of identical spheres; reciprocal lattice is BCC.

(todo) Tetragonal disphenoid honeycomb

:   Uniform honeycomb and body-centered cubic (BCC or $A_3^*$ or $D_3^*$) lattice; mean coordination 8; space group $Im\bar{3}m$ (229); vertex-, face-, and cell-transitive; reciprocal lattice is FCC; ideal k-space samples in $R^3$.

(todo) Rhombic dodecahedral honeycomb

:   Non-uniform Voronoi honeycomb of FCC lattice; mean coordination 5.5; space group $Fm\bar{3}m$ (225); edge-, face-, and cell-transitive; 3-space parallelohedron.

(todo) Bitruncated cubic honeycomb

:   Uniform Voronoi honeycomb of BCC lattice; mean coordination 4; space group $Im\bar{3}m$ (229); vertex-, edge-, and face-transitive; 3-space permutohedron; best-known ideal foam (Kelvin problem) for a century, then superseded by the Weaire–Phelan honeycomb.

# Experimental Design, Implementation, and Validation

**(todo)** State reasoning for finding the smallest-possible integer representation, if any; show unstable configuration; walk through construction of *A15* and surrounding it with the Weaire–Phelan honeycomb, identifying its minimum prescale factor along the way; relate integer scaling to fraction-like floating-point definition from introduction; add the next layer of lattice and additional prescale due to separation distance; binary splits of this final prescale are “binary" scales, multiples of these splits are “stable" scales, and everything else is “unstable"; bits shuffle cleanly between range and density, facilitating the definition of subspaces; show Tetrastix at Weaire–Phelan’s prescale and frame volume difference as an error domain.

## *A15*.py’s Experimental Design and Exploration

**(todo)** Explain *A15*.py generation process, reasoning, capabilities, and assertions; describe main image and such things as $N_1$ (cell width), $\epsilon_N$, $\epsilon_\delta$, $\epsilon_\Delta$, and $\epsilon$.

<figure id="fig-main">
<p><img src="fig-main.png" alt="image" /><span id="fig-main" label="fig-main"></span></p>
</figure>

## Validation Techniques and Statistical Analysis

<figure id="fig-hist">
<p><img src="fig-histb.png" alt="image" /><img src="fig-hists.png" alt="image" /><img src="fig-histu.png" alt="image" /></p>
<figcaption>Example <em>binary</em> (left) <em>stable</em> (middle) and <em>unstable</em> (right) configurations.<span id="fig-hist" label="fig-hist"></span></figcaption>
</figure>

**(todo)** Explain *A15*.py histogram bar graph and its epsilons; show unstable configurations generating a smattering of epsilons and unused gaps; contrast this with gapless, sequential configurations that always generate a limited number of epsilons (stable) or an exact number (binary).

# Results Interpretation, Insights, and Limitations

*A15* has a long and storied history, from holding the high-temperature superconductor record for decades, to its close association with other interesting structures—such as the Weaire–Phelan honeycomb and the Tetrastix prism—each with their own unique qualities. The Weaire–Phelan honeycomb, in its relaxed, non-polyhedral “bubble" form (combinatorially equivalent to the polyhedral honeycomb), consistently yields highly-isotropic measurements from different physical quantities, including thermal expansion rate, compressional load transfer, photonic wave propagation, and quantum noise distribution. *A15*’s position is further reinforced through its crystallographic space group properties, such as an exceptionally high coordination number (mean of 13.5 connections per site), second-highest symmetry order (24), and maximal intersection with the *non-crystallographic* $I_h$ group. Centered on the Weaire–Phelan honeycomb—the best-known equal-volume partitioner of 3D space and lowest-energy solution to the Kelvin problem—and simultaneously compatible with the Tetrastix prism—an attractive alternative to Weaire–Phelan-based discretization when trading spatial inaccuracy for performance is acceptable or desirable—makes *A15* uniquely qualified for partitioning interactive 3D space.

## Interpretation of Findings and In-Depth Insights

This research presents *A15* as the ideal, coordinated partitioner of shared, interactive virtual space. It highlights *A15*’s innate mapping to hardware floating-point representations, and its compatibility with global measurement systems. It details two, high-dimensional addressing schemes, and explains their differences in terms of lattice centering, chiral balancing, and interlocking extrema. *A15*’s numeric floating-point stability is confirmed, and three distinct classifications—*binary*, *stable*, and *unstable*—are identified. It asserts the binary floating-point stability of the Weaire–Phelan honeycomb, the Tetrastix prism, and *A15* itself, using both a series of mathematical statements and generated empirical evidence. It constructs well-defined, heterogeneous environments, and characterizes how additional bits either double the range or double the density. Original code is shared in full, alongside detailed examples and documentation. This outcome, combined with *A15*’s body of existing materials research, makes *A15* an excellent candidate for *reshaping* the metaverse.

## Noteworthy Limitations and Topics of Concern

Regardless of addressing scheme, most on-lattice symmetry operations require follow-up translations to maintain *A15*’s desirable invariants. *A15*’s primitive unit cell contains eight valid sites—known as its basis, or crystal motif—but only one is also a valid lattice point.

*A15* is a linear space partitioner. Creating massive, open-world virtual environments is less straightforward because *A15* spaces are small, irregularly sized, and end abruptly. Federating virtual real estate in a fifty-story skyscraper might require independent *A15* spaces per floor, or per delegable unit. This matches expectations after reclaiming unused dynamic floating-point range as memory savings, but fast-paced transitions at spatial joins are potential sources of bugs, complexity, and overhead.

This research generates *A15* crystal-lattice cuboids with rectangular faces, and with edge lengths in proportion to each dimension’s allocation of the total bit space. Cuboids best-encode spaces like playing fields and office buildings—spaces with high average utilization, and clear, axes-aligned boundaries—and are less efficient for spaces with arbitrary terrain, irregular boundaries, or internal holes. Such features generate effectively unreachable pockets of addressable space. However, compared to the unused dynamic range of floats, this underutilization is much simpler to quantify.

Cited Weaire–Phelan research uses a relaxed form of the Weaire–Phelan honeycomb to satisfy Plateau’s laws and the constraints of the Kelvin problem. This structure is combinatorially equivalent to the polyhedral form, albeit with softer angles and perfectly equal-volume pyritohedra and tetradecahedra. Results from these studies may not transfer cleanly. However, discrepancies are limited to Weaire–Phelan-based discretization claims and not *A15* itself.

This research facilitates the long-term storage and analysis of human-generated spatial tracking data—personally identifiable information (PII) with both legal and ethical requirements—and *demands* that implementers honor and regard it with the utmost care and respect.

# Immediate Potential and Future Prospects

## Universal Compatibility and Widespread Applicability

**(todo)** L/R balancing for both Unreal and Unity; reliably mutate, consolidate, or replay from the edge; clean mapping into global measurement systems; spatial range limits are implicit anti-cheating mechanisms; range limiting combined with other symmetry-reliant techniques expected to save more memory than *baseline* of 50%.

## Future Investigations and Potential for Immersive Experience

**(todo)** Expecting to train AI bots; expecting to capture high speed games; resultant libraries not limited to metric spaces, could also apply to eg. function spaces; investigate rotating *A15* favorably with respect to expected traversal patterns, eg. Miller index $(111)$, and boost effective isotropy; consider ways to efficiently carve terrain, holes, and non-planar or axes-misaligned boundaries; seek esoteric 3-space objects outside crystallography and convex geometry, eg. stars, gyroids, quasicrystals, plesiohedra, Delone sets, and centroid network of the Laves graph; continue search for higher-dimensional polytopes to fold and project into 3D space.

# Supplementary Materials and Extended Analyses

## Comprehensive Code Guidelines and Replication Protocols

**(todo)** Share commands to rebuild figures in this research; illuminate *A15*.py’s dark, undocumented corners; note the binary scale’s clean steppings compared to stable scaling; in some addressing schemes, the L/R orientation of major sites is independent from that of minor sites, with four or more valid, internal orientations possible; when more than one internal orientation is available, this research arbitrarily chooses one; all parties must know, negotiate, or discover an orientation before transmitting coordinates; some orientations might be more favorable than others, both at the boundary and in the bulk.

## Additional Data, Visualizations, and Mathematical Exploration

**(todo)** <https://infima.space>; A15.py itself; proofs of different assertions made, eg. 5v5 100 kbps per user baseline; tables of vertices; list figures of interesting heterogenous environments; note how surrounding minor sites with a 7/5 pyritohedron generates correct Weaire–Phelan spacing and possible dual space; support and recommend $2^{-6}$ as the preferred default scale; note coincidence of *A15* integer scale, Tetrastix unit scale, and favorable $2^{-6}$ overall scale; demonstrate generating other crystal lattices; unit of least precision; meets meet joins join meets meet joins join...

<div class="acknowledgements">

**(todo)** Acknowledgements.

</div>

<div id="refs" class="references csl-bib-body hanging-indent" entry-spacing="0">

<div id="ref-todo" class="csl-entry">

References. 2023. “(Todo).” *(Todo)*.

</div>

</div>

[^1]: Exponents are effectively signed—*bias* is an implementation detail—and $2^{-n}$ is $^1/_{2^n}$.
