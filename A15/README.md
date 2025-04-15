# Executive Summary

The *A15* encoding framework addresses a fundamental challenge in
spatial computing: achieving deterministic, bandwidth-efficient
coordinate representation across heterogeneous systems. By mapping
continuous 3D coordinates onto a crystallographically-inspired integer
lattice, this approach establishes precise spatial representation within
the inherently imprecise world of floating-point arithmetic.

At its core, the framework leverages the unique properties of the *A15*
phase structure ($`\beta`$–$`W`$), a crystal structure with exceptional
isotropy and binary-compatible coordinates. Through a carefully designed
multi-stage scaling pipeline that prioritizes integer representation
during geometric construction, it guarantees bit-identical spatial
encoding when operating within defined stable scaling regimes
($`\epsilon_\Delta = 0`$). This eliminates the representation errors
that plague floating-point coordinates, ensuring consistent behavior
across different hardware, operating systems, and compiler
implementations.

For networked virtual environments—particularly competitive VR,
collaborative simulations, and distributed computing—this translates
directly to several practical benefits:

- **Network Efficiency:** Reduces coordinate transmission size by
  approximately 50% compared to raw float vectors, enabling higher
  update rates within limited bandwidth.

- **Deterministic Guarantee:** Enables verifiable event sequences,
  accurate replays, and fair competition across heterogeneous client
  devices.

- **Geometric Fidelity:** Preserves spatial relationships with minimal
  directional bias owing to the structure’s high coordination (13.5) and
  near-icosahedral local ordering.

The framework includes robust tools for verifying numerical stability
through histogram analysis, enabling developers to confidently implement
deterministic spatial systems while maintaining flexibility in
coordinate scale selection appropriate to their application domain.

<figure id="fig-intro">
<img src="fig-intro.png" />
<figcaption>Illustration of the <em>A15</em> structure (<span
class="math inline">\(\beta\)</span>–<span
class="math inline">\(W\)</span>) at three different scales, generated
using <code>A15.py</code> (<code>fig-intro.png.txt</code>). The
visualization highlights the constituent pyritohedral and
tetradecahedral elements, demonstrating how the structure maintains its
fundamental <span class="math inline">\(Pm\bar{3}n\)</span> space group
symmetry while revealing increasing detail at finer
resolutions.</figcaption>
</figure>

# Introduction: Floats, Fairness, and the A15 Foundation

Binary floating-point numbers, commonly known as “floats” and largely
governed by the IEEE 754 standard (IEEE 2019), are the computational
workhorses for representing real numbers. They offer a vast dynamic
range essential for graphics and science, yet this reach often comes at
the cost of precision. Floats are notorious for providing *approximately
correct* answers[^1], representing most values inexactly. These tiny,
fundamental representation errors, inherent in their structure
(<a href="#eq-float-representation" data-reference-type="ref+Label"
data-reference="eq-float-representation">[eq-float-representation]</a>),
can accumulate and interact, particularly challenging the deterministic
consistency required for networked virtual reality (VR) or distributed
simulations. This research confronts these challenges head-on, proposing
an alternative spatial representation rooted in the unique
crystallographic properties of the *A15* phase structure
($`\beta`$–$`W`$, space group $`Pm\bar{3}n`$, No. 223 (Aroyo 2016)), as
illustrated in <a href="#fig-intro" data-reference-type="ref+Label"
data-reference="fig-intro">1</a>. Furthermore, establishing a robust
foundation for spatial representation with inherent numerical stability
and structural consistency extends beyond immediate interactive
fidelity. Such a foundation offers compelling advantages for long-term
data archival, application-agnostic interoperability across diverse
platforms, and the efficient decomposition and distribution of spatial
computations in parallel and edge computing environments.

## The Challenge of Floating-Point Precision

Standard normalized binary floating-point numbers derive their value
from three components: a sign bit, a biased exponent, and a significand
(mantissa). Conceptually[^2], their value relates to:
``` math
\label{eq-float-representation}
    \text{value} = (-1)^{\text{sign}} \times (1 + \text{significand}) \times 2^{\text{exponent}}
```
The fixed-size significand offers constant *relative* precision within
an exponent range, while the exponent scales the value logarithmically.
This design excels at exactly representing powers of two and fractions
whose denominators are solely powers of two, but only a finite subset
thereof. Consequently, most decimal values and even many simple
fractions (e.g., $`1/10`$) are merely approximated (Goldberg 1991).

While seemingly minuscule, these representation errors can compound,
critically affecting reproducibility. In distributed systems like
multiplayer games or collaborative VR, this becomes paramount. Identical
logical inputs processed on different hardware architectures, operating
systems, or even with different compiler optimizations can yield subtly
divergent floating-point results. This variance breaks simulation
consistency, undermines fairness in competition (Claypool and Claypool
2006), complicates debugging and verification, and fundamentally opposes
the requirement for determinism.

These challenges manifest most acutely in latency-sensitive,
spatially-complex applications where:

- Small coordinate discrepancies can produce visible jitter or
  desynchronization

- State consistency must be maintained across heterogeneous client
  devices

- Replays and verifiable event sequences require bit-exact reproduction

- Network bandwidth constraints demand efficient coordinate
  representation

Yet, widespread hardware acceleration makes IEEE 754 formats—primarily
binary32 (32, single-precision) and binary64 (64,
double-precision)—ubiquitous in modern CPUs and GPUs. Binary32, in
particular, remains vital for performance-sensitive applications like
game development and VR, establishing a key baseline for efficiency
comparisons. This work, therefore, accepts the practical necessity of
floats but seeks to structure their use, mitigating risks to determinism
via the disciplined partitioning and verifiable scaling offered by the
*A15* framework detailed in
<a href="#sec-framework-design" data-reference-type="ref+Label"
data-reference="sec-framework-design">2</a>.

## The A15 Phase Structure: A Crystallographic Foundation

Effectively partitioning virtual space requires a measure of geometric
“fairness”—a balance between **isometry** (preserving distances between
corresponding points) and **isotropy** (uniformity of properties across
different directions). Discretizing continuous space onto a lattice
while preserving both ensures that spatial relationships remain
consistent and free from directional bias. While regular lattices
achieve perfect isometry through translational symmetry, attaining high
isotropy is constrained by fundamental geometric principles. The
continuous rotational symmetry group in 3D, $`SO(3)`$, implies that
highly isotropic structures should appear locally spherical. However,
the imposition of discrete translational symmetry restricts the allowed
rotational symmetries to only 1-, 2-, 3-, 4-, and 6-fold axes. This is
the essence of the **crystallographic restriction theorem** (Ashcroft
and Mermin 1976), which explicitly forbids the global 5-fold symmetry
characteristic of highly isotropic polyhedra like the icosahedron
(associated with the $`I_h`$ point group symmetry (Coxeter 1973)) in
periodic lattices.

The *A15* structure ($`Pm\bar{3}n`$, No. 223 (Aroyo 2016)) navigates
this restriction with notable elegance. It is based on the simple
primitive cubic ($`cP`$) Bravais lattice but incorporates a complex
basis (or motif) containing two distinct types of crystallographic sites
within its conventional unit cell (Frank and Kasper 1958, 1959), as
illustrated in <a href="#fig-cell2" data-reference-type="ref+Label"
data-reference="fig-cell2">2</a>:

- 25% are C12 sites: 12-coordinated, occupying Wyckoff position 2a
  (e.g., at fractional coordinates $`(0, 0, 0)`$ and
  $`(1/2, 1/2, 1/2)`$). Their local coordination environment exhibits
  characteristics related to pyritohedral symmetry.

- 75% are C14 sites: 14-coordinated, occupying Wyckoff position 6d
  (e.g., at $`(1/4, 1/2, 0)`$ and its symmetry equivalents (Aroyo
  2016)). Their local environment resembles a tetradecahedron (a
  14-faced polyhedron).

This intricate arrangement results in a high average coordination
number, calculated as
$`Z_{avg} = (0.25 \times 12) + (0.75 \times 14) = 13.5`$, suggesting a
densely connected local structure conducive to efficient packing. The
local symmetry around these sites is described by the crystallographic
point group $`T_h`$ ($`m\bar{3}`$, order 24). Significantly, $`T_h`$ is
a maximal subgroup common to both the full cubic symmetry group $`O_h`$
(order 48) and the non-crystallographic icosahedral group $`I_h`$ (order
120) (Coxeter and Moser 1972). While lacking the forbidden 5-fold axes
of $`I_h`$, $`T_h`$ retains key 3-fold rotational symmetries present in
icosahedral structures. This unique symmetry allows the *A15* structure
to incorporate significant near-icosahedral local ordering, boosting
local isotropy considerably compared to simpler cubic structures, all
while maintaining the long-range periodicity required by
crystallography. The structure also features alternating left- and
right-handed local environments around the 6d sites, adding further
complexity relevant to implementation
(<a href="#subsubsec-limits-handedness" data-reference-type="ref+Label"
data-reference="subsubsec-limits-handedness">3.3.6</a>).

Crucially for this work, the canonical definition of the *A15* structure
relies on fractional coordinates inherently suitable for binary
representation: $`0`$, $`1/4`$, and $`1/2`$. These values are exactly
representable in base-2 systems. When these fractional coordinates are
mapped onto an integer lattice (conceptually, by scaling the unit cell
coordinates by 4, aligning with the derivation of the 96-unit baseline
in <a href="#subsubsec-scaling-baseline" data-reference-type="ref+Label"
data-reference="subsubsec-scaling-baseline">2.3.3</a>), two
characteristic squared site-to-site distances emerge within this integer
framework: $`d^2 = 4`$ (yielding distance 2) and $`d^2 = 5`$ (yielding
distance $`\sqrt{5}`$). The presence of $`\sqrt{5}`$, the hypotenuse of
a fundamental 2:1 right triangle, connects the structure’s geometry
directly to the golden ratio $`\phi = (1+\sqrt{5})/2`$, further
reflecting the embedded near-icosahedral characteristics. This intrinsic
numerical simplicity and compatibility with base-2 representation
ensures that the *A15* structure can be defined and scaled with
precision, forming a robust foundation for the numerically stable
encoding framework developed herein
(<a href="#subsec-stability" data-reference-type="ref+Label"
data-reference="subsec-stability">2.4</a>).

<figure id="fig-cell2">
<img src="fig-cell2.png" />
<figcaption>Components illustrating the local arrangement within a
left-handed <span class="math inline">\(^1/_2\)</span> unit cell
fragment of the <em>A15</em> structure. Shows C12 (Wyckoff 2a, center of
pyritohedral environment) and C14 (Wyckoff 6d, center of tetradecahedral
environment) sites demonstrating the local complexity within the overall
<span class="math inline">\(Pm\bar{3}n\)</span> symmetry. Generated via
<code>A15.py</code> (<code>fig-cell2.png.txt</code>).</figcaption>
</figure>

## Local Discretization Methods: WPH and TSP

While the *A15* structure defines the target lattice points for
coordinate encoding, mapping arbitrary continuous coordinates from a
virtual environment onto these discrete identifiers requires a
well-defined *local discretization method*. This method effectively
defines the boundaries of the region (the “cell”) surrounding each *A15*
site; any continuous point falling within a given cell is quantized or
mapped to that specific site’s identifier. This framework primarily
considers two such methods, both intrinsically linked to the *A15*
structure, sharing its $`Pm\bar{3}n`$ space group symmetry, and
implemented within the accompanying `A15.py` software
(<a href="#fig-wp-ts" data-reference-type="ref+Label"
data-reference="fig-wp-ts">3</a>):

<div class="description">

This research utilizes the geometric polyhedral form of the
Weaire–Phelan structure as a partitioning method. It divides space using
cells of two distinct shapes: a 12-faced pyritohedron and a 14-faced
tetradecahedron, arranged in a precise 1:3 ratio (Weaire and Phelan
1994). Significantly, these flat-faced polyhedra constitute the Voronoi
decomposition for the A15 lattice points; the pyritohedron is the
Voronoi cell for the C12 (Wyckoff 2a) site, and the tetradecahedron is
the Voronoi cell for the C14 (Wyckoff 6d) site. This relationship makes
the Weaire–Phelan Honeycomb geometry a canonical choice for
discretization, as each cell naturally defines the region of space
closest to its corresponding A15 lattice site. These cell shapes
correspond topologically to the local coordination environments of the
C12 and C14 sites within the *A15* structure, respectively. This
partitioning is notable for its high degree of isotropy, closely
mirroring the symmetry characteristics of the underlying lattice. (It is
important to distinguish this geometric, space-filling honeycomb from
the related but distinct relaxed, minimal-surface structure which
famously provided a counter-example to Kelvin’s conjecture on foam
partitioning (Thomson 1887; Kusner and Sullivan 1996; Weaire and Hutzler
2001).)

A structurally simpler alternative partitioning method, activated within
`A15.py` via the `-stix` configuration option. This method employs
axis-aligned planar faces, effectively dividing space into cubic blocks
centered on the *A15* lattice sites. While offering computational
advantages for certain operations (such as point-in-cell tests) due to
its simpler geometry, this approach generally exhibits lower isotropy
compared to the Weaire–Phelan Honeycomb method.

</div>

Crucially, both the Weaire–Phelan Honeycomb and the Tetrastix Prism
serve as interchangeable local discretization methods within this
framework. They define the specific geometry used to quantize continuous
space around each *A15* lattice site. The choice influences the precise
shape of these local regions and boundaries, impacting factors like the
resulting partition’s isotropy and the computational cost of the
quantization step itself. However, the fundamental coordinate identifier
that is ultimately transmitted or stored remains based on the position
within the underlying *A15* lattice, regardless of which local method is
used for the mapping.

<figure id="fig-wp-ts">
<p><img src="fig-wp.png" alt="image" /> <img src="fig-ts.png"
alt="image" /></p>
<figcaption>Visualization of local discretization methods related to the
underlying <em>A15</em> structure. Left: The geometric polyhedral
variant of the <strong>Weaire–Phelan Honeycomb</strong>, featuring
pyritohedra (12 faces) and tetradecahedra (14 faces). Right: The simpler
<strong>Tetrastix Prism</strong> partitioning using axis-aligned planar
faces resulting in cubic blocks centered on <em>A15</em> sites. Both
share the <span class="math inline">\(Pm\bar{3}n\)</span> space group
symmetry. Generated using <code>A15.py</code>
(<code>fig-wp.png.txt</code> and
<code>fig-ts.png.txt</code>).</figcaption>
</figure>

## Clarification on Weaire-Phelan Variants

It is essential to emphasize that this research specifically utilizes
the **geometric polyhedral form** of the Weaire-Phelan structure as a
partitioning method for the A15 framework. This geometric variant
consists of flat-faced polyhedra (pyritohedra and tetradecahedra)
arranged in a precise 1:3 ratio, constituting the Voronoi decomposition
for the A15 lattice points.

This geometric WPH variant should be explicitly distinguished from the
related but distinct **minimal-surface foam approximation** of
Weaire-Phelan, which is widely known in the physics literature for
providing a counter-example to Kelvin’s conjecture on foam partitioning
(Thomson 1887; Kusner and Sullivan 1996; Weaire and Hutzler 2001). The
minimal-surface variant features curved interfaces minimizing surface
area under specific physical constraints, resulting in subtle geometric
differences from the polyhedral form used here.

<figure id="fig-wp-variants">

<figcaption>Conceptual comparison between the geometric polyhedral form
of the Weaire-Phelan structure used in this framework (left) and the
minimal-surface foam approximation known from physics literature
(right). The geometric form uses flat-faced polyhedra with exact vertex
coordinates, while the minimal-surface variant features curved
interfaces that minimize surface area.</figcaption>
</figure>

The practical implications of this distinction are significant for
implementation:

- The geometric WPH used in this framework can be precisely represented
  using flat-faced polyhedra with exact vertex coordinates, facilitating
  deterministic computational processing.

- Point-in-cell tests against the geometric WPH can be implemented using
  standard computational geometry techniques for polyhedra, avoiding the
  complexity of curved surface representations.

- The structural rigidity of the geometric WPH ensures consistent
  partitioning regardless of cell contents or external influences, in
  contrast to physical foam structures that might deform based on
  pressure differentials.

In contrast, the minimal-surface WPH would require more complex surface
representations, potentially introducing additional numerical challenges
contradictory to the framework’s determinism goals. Therefore, all
references to the Weaire-Phelan Honeycomb (WPH) within this manuscript
specifically denote the geometric polyhedral form, not the
minimal-surface variant from foam physics.

This clarification is particularly important because the minimal-surface
WPH has gained broader recognition through its architectural application
in structures like the Beijing National Aquatics Center (the "Water
Cube") for the 2008 Olympics, potentially creating confusion about which
variant is employed in this spatial encoding framework.

## Distinction from Traditional Spatial Indexing

The *A15*-based partitioning approach presented here relates to, yet
fundamentally differs from, traditional spatial indexing techniques
commonly employed in databases, geographic information systems (GIS),
and various domains of computer graphics (e.g., octrees (Finkel and
Bentley 1974), R-trees (Guttman 1984), k-d trees (Bentley 1975); see
(Samet 1990) for a comprehensive overview). Conventional indexing
methods typically utilize data-driven, adaptive strategies; their
partitioning boundaries often adjust dynamically based on the
distribution and density of existing spatial objects, with the primary
goal of optimizing query performance (such as range searches or nearest
neighbor lookups) for that specific dataset.

In stark contrast, the *A15* framework imposes a predetermined, regular,
crystallographically-inspired structure onto the virtual space itself,
largely independent of the dynamic content residing within it. This
*structure-first* methodology prioritizes the creation of a universal,
efficient, and geometrically sound fabric for spatial representation and
encoding. The key distinctions are summarized in
<a href="#tab-comparison-indexing" data-reference-type="ref+Label"
data-reference="tab-comparison-indexing">1</a>. The *A15* approach
deliberately trades the dynamic query optimization focus of traditional
indexing for significant advantages in predictability, communication
efficiency (via compact identifiers that implicitly encode the regular
structure), guaranteed numerical determinism (when using stable scaling
regimes), and inherent geometric consistency derived from
crystallographic symmetry—properties particularly valuable for networked
virtual environments demanding shared spatial understanding, robustness
against floating-point discrepancies, and efficient state
synchronization.

<div id="tab-comparison-indexing">

<table>
<caption>Comparison of Spatial Organization Approaches.</caption>
<thead>
<tr>
<th style="text-align: left;"><strong>Characteristic</strong></th>
<th style="text-align: left;"><strong>Traditional Spatial
Indexing</strong><br />
(e.g., Octree <span class="citation" data-cites="Finkel1974">(Finkel and
Bentley 1974)</span>, R-tree <span class="citation"
data-cites="Guttman1984">(Guttman 1984)</span>, k-d tree <span
class="citation" data-cites="Bentley1975">(Bentley 1975)</span>)</th>
<th style="text-align: left;"><strong><em>A15</em> Crystallographic
Partitioning</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">Main Objective</td>
<td style="text-align: left;">Efficient query/retrieval of existing,
often arbitrary, data</td>
<td style="text-align: left;">Uniform, deterministic spatial
discretization and efficient coordinate encoding</td>
</tr>
<tr>
<td style="text-align: left;">Partition Strategy</td>
<td style="text-align: left;">Data-driven; adaptive boundaries; often
hierarchical; structure follows data</td>
<td style="text-align: left;">Structure-driven; regular lattice;
predetermined boundaries (via WPH/TSP local methods); data follows
structure</td>
</tr>
<tr>
<td style="text-align: left;">Symmetry Awareness</td>
<td style="text-align: left;">Generally low; structure adapts to data,
often breaking ambient symmetries</td>
<td style="text-align: left;">High; leverages crystallographic symmetry
(<span class="math inline">\(T_h\)</span>, <span
class="math inline">\(Pm\bar{3}n\)</span>) of the imposed lattice</td>
</tr>
<tr>
<td style="text-align: left;">Coordinate Handling</td>
<td style="text-align: left;">Typically preserves input precision (e.g.,
float coordinates)</td>
<td style="text-align: left;">Quantizes coordinates to discrete
representations (<em>A15</em> identifiers, often integers)</td>
</tr>
<tr>
<td style="text-align: left;">Memory Usage</td>
<td style="text-align: left;">Variable, depends on data
density/distribution; includes tree/node overhead</td>
<td style="text-align: left;">Predictable based on defined scale/volume;
highly efficient storage using integer identifiers</td>
</tr>
<tr>
<td style="text-align: left;">Isotropy</td>
<td style="text-align: left;">Variable; depends heavily on data
distribution and algorithm specifics</td>
<td style="text-align: left;">High local isotropy inherent in the
<em>A15</em> structure, especially when using Weaire–Phelan Honeycomb
partitioning</td>
</tr>
<tr>
<td style="text-align: left;">Neighborhood Info</td>
<td style="text-align: left;">Requires explicit tree traversal or
complex range/proximity queries</td>
<td style="text-align: left;">Implicit, regular neighbor relationships
directly derivable from the lattice structure</td>
</tr>
<tr>
<td style="text-align: left;">Temporal Stability</td>
<td style="text-align: left;">Partitioning structure can change
significantly as data moves or updates</td>
<td style="text-align: left;">Fixed partitioning grid provides temporal
coherence for coordinate mapping (though content moves within)</td>
</tr>
<tr>
<td style="text-align: left;">Determinism Guarantee</td>
<td style="text-align: left;">Generally low; susceptible to float
variance affecting boundary tests/traversal</td>
<td style="text-align: left;">High; guaranteed bit-level consistency
achievable with stable scaling regimes (<span
class="math inline">\(\epsilon_\Delta = 0\)</span>)</td>
</tr>
<tr>
<td style="text-align: left;">Main Use Case</td>
<td style="text-align: left;">Databases, GIS <span class="citation"
data-cites="Samet1990">(Samet 1990)</span>, dynamic collision detection,
view frustum culling</td>
<td style="text-align: left;">Foundational spatial fabric for
deterministic VR/simulations, compact coordinate encoding, networked
state consistency, verifiable replays</td>
</tr>
</tbody>
</table>

</div>

## Comparative Context: Other Symmetric Structures

The 230 crystallographic space groups categorize all possible ways to
combine periodic translational symmetry (defined by the 14 Bravais
lattices) with rotational and reflectional symmetries (defined by the 32
crystallographic point groups) in 3D space (Aroyo 2016). Situating the
*A15* structure (space group $`Pm\bar{3}n`$, No. 223) within this
framework highlights its distinctive characteristics relevant to packing
and partitioning
(<a href="#tab-spacegroups" data-reference-type="ref+Label"
data-reference="tab-spacegroups">2</a>).

<div id="tab-spacegroups">

| **Structure Type** | **No.** | **HM Symbol** | **Point Group** | **Order** | **Bravais Lattice** |
|:---|:--:|:--:|:--:|:--:|:---|
| ***A15* Type** | **223** | $`\boldsymbol{Pm\bar{3}n}`$ | $`\boldsymbol{m\bar{3}}`$ $`\boldsymbol{(T_h)}`$ | **24** | **Primitive (cP)** |
| Simple Cubic (SC) | 221 | $`Pm\bar{3}m`$ | $`m\bar{3}m`$ $`(O_h)`$ | 48 | Primitive (cP) |
| BCC Type (e.g., W) | 229 | $`Im\bar{3}m`$ | $`m\bar{3}m`$ $`(O_h)`$ | 48 | Body-Centered (cI) |
| FCC Type (e.g., Cu) | 225 | $`Fm\bar{3}m`$ | $`m\bar{3}m`$ $`(O_h)`$ | 48 | Face-Centered (cF) |

Comparison of Relevant Cubic Space Groups.

</div>

Data sourced from International Tables for Crystallography, Vol A (Aroyo
2016). HM Symbol: Hermann-Mauguin notation.

As shown, *A15*’s $`Pm\bar{3}n`$ space group utilizes the primitive
cubic ($`cP`$) Bravais lattice but applies the $`T_h`$ point group
symmetry (order 24). This group possesses lower overall symmetry than
the full cubic $`O_h`$ group (order 48) characteristic of the
highest-symmetry simple cubic, body-centered cubic (BCC, $`cI`$), and
face-centered cubic (FCC, $`cF`$) structures. However, as discussed
(<a href="#subsec-intro-a15" data-reference-type="ref+Label"
data-reference="subsec-intro-a15">1.2</a>), the specific $`T_h`$
symmetry is significant because it represents a maximal crystallographic
subgroup linking cubic ($`O_h`$) and non-crystallographic icosahedral
($`I_h`$) symmetries (Coxeter and Moser 1972). This unique combination
allows the *A15* structure to accommodate its complex basis with
distinct C12 and C14 sites, enabling its exceptionally high average
coordination number (13.5), which relates to efficient local packing and
high connectivity.

When compared to common space-filling honeycombs and related sphere
packings, the *A15* structure, particularly when coupled with
the Weaire–Phelan Honeycomb partitioning method, offers further
distinctions relevant to this framework:

<div class="description">

The body-centered cubic lattice (BCC, coordination number Z=8) and its
dual, the bitruncated cubic honeycomb (Kelvin’s structure (Thomson
1887)), have lower fundamental coordination than *A15* (13.5),
suggesting a less densely connected local environment.

The face-centered cubic lattice (FCC, coordination Z=12) represents the
densest packing of identical spheres (Conway and Sloane 1999) and is
related to the rhombic dodecahedral honeycomb. While achieving maximal
coordination for identical points, *A15* surpasses this average by
utilizing two distinct, efficiently arranged site types.

The geometric Weaire–Phelan Honeycomb, used as a local discretization
method, shares the $`Pm\bar{3}n`$ space group and its cell types
correspond topologically to the A15 site environments. This ensures the
partitioning boundaries naturally align with the symmetries and
neighborhood relationships of the underlying A15 encoding lattice, a
unique synergy not offered by partitions derived from simpler lattices.

</div>

In summary, the unique combination offered by the *A15* structure—its
high mean coordination, significant local isotropy via $`T_h`$ symmetry,
a direct relationship to the compatible Weaire–Phelan Honeycomb
partitioning structure, and inherently binary-suitable base
coordinates—makes it exceptionally well-suited for the deterministic,
efficient, and geometrically sound spatial encoding framework proposed
here.

## Comprehensive Lattice Comparison

The selection of *A15* as the foundation for this framework was not
arbitrary but resulted from a systematic evaluation of alternative
crystallographic structures.
Table <a href="#tab-lattice-comparison" data-reference-type="ref"
data-reference="tab-lattice-comparison">3</a> provides a comprehensive
comparison of *A15* against common alternative lattices for spatial
partitioning, highlighting the specific properties that make it
exceptionally well-suited for deterministic spatial encoding.

<div id="tab-lattice-comparison">

| **Lattice** | **Coord. \#** | **Isotropy** | **Binary?** | **Uniform?** | **Complex?** |
|:---|:--:|:--:|:--:|:--:|:--:|
| Simple Cubic (SC) | 6 | Low | High | Perfect | Low |
| Body-Centered Cubic (BCC) | 8 | Medium | Medium | Perfect | Medium |
| Face-Centered Cubic (FCC) | 12 | High | Low | Perfect | Medium |
| ***A15* (A15)** | **13.5 avg** | **High** | **High** | **Dual-type** | **High** |

Comprehensive Comparison of Lattice Structures for Spatial Partitioning

</div>

Each property contributes significantly to the lattice’s suitability for
deterministic spatial encoding:

Coordination Number:  
The average number of nearest neighbors for each lattice point directly
relates to the connectivity and local packing efficiency. *A15*’s
exceptional mean coordination number (13.5) enables a more densely
connected local structure than even the FCC lattice (12), often
considered optimal for uniform spheres.

Isotropy:  
The uniformity of geometric properties across different directions
affects how consistently the lattice represents spatial relationships
regardless of orientation. While SC exhibits significant directional
bias, both FCC and *A15* achieve high isotropy. However, *A15* maintains
this isotropy specifically through its $`T_h`$ point group symmetry,
which incorporates key aspects of icosahedral ordering
(<a href="#subsec-intro-a15" data-reference-type="ref+Label"
data-reference="subsec-intro-a15">1.2</a>).

Binary Compatibility:  
The natural alignment of the structure’s coordinates with binary
representation directly impacts numerical stability. *A15*’s basis site
coordinates (involving only factors of $`0`$, $`1/4`$, and $`1/2`$) are
perfectly compatible with binary representation, unlike the irrational
coordinates required for optimal FCC packing.

Cell Uniformity:  
While SC, BCC, and FCC consist of identical unit cells, *A15*
incorporates two distinct site types (C12 at Wyckoff 2a and C14 at
Wyckoff 6d) in a precise 1:3 ratio. This heterogeneity enables its
exceptional mean coordination without sacrificing periodicity or
computational tractability.

Complexity:  
Implementation complexity ranges from the straightforward SC to the more
involved *A15*.

<figure id="fig-lattice-comparison">

<figcaption>Comparison of crystal lattice structures for spatial
partitioning. From left to right: Simple Cubic (SC) with coordination
number 6, Body-Centered Cubic (BCC) with coordination number 8,
Face-Centered Cubic (FCC) with coordination number 12, and A15 structure
with average coordination number 13.5. The A15 structure combines high
isotropy with binary-friendly coordinates, making it especially suitable
for deterministic spatial encoding.</figcaption>
</figure>

The key insight from this comparison is that *A15* occupies a unique
position in the design space: it achieves exceptionally high
coordination (surpassing even FCC) while maintaining perfect
compatibility with binary representation—a combination not available in
other crystallographic structures. This makes it particularly
well-suited for applications requiring both high isotropy and guaranteed
determinism.

While simpler lattices like SC offer lower implementation complexity,
they suffer from significant directional bias that can adversely affect
the fairness of spatial representation. BCC improves upon isotropy but
still falls short of the near-icosahedral characteristics of *A15*. FCC
achieves excellent isotropy and packing density but lacks the
binary-friendly coordinates essential for deterministic representation.

The *A15* structure’s combination of high coordination, significant
local isotropy, and binary-compatible coordinates provides a robust
foundation for the deterministic, efficient, and geometrically sound
spatial encoding framework proposed in this research.

# The *A15* Encoding Framework: Design and Implementation

The heart of this research lies in translating the idealized
crystallographic description of the *A15* structure
(<a href="#subsec-intro-a15" data-reference-type="ref+Label"
data-reference="subsec-intro-a15">1.2</a>) into a practical, numerically
robust system for partitioning continuous 3D space and encoding
coordinates. This requires establishing an appropriate and consistent
scale, mapping continuous coordinates to discrete lattice identifiers,
and ensuring the process avoids the pitfalls of floating-point
approximation
(<a href="#subsec-intro-floats" data-reference-type="ref+Label"
data-reference="subsec-intro-floats">1.1</a>). The framework achieves
this through a carefully staged approach to coordinate scaling,
primarily operating within an internal integer coordinate system during
geometric construction to preserve fidelity and defer floating-point
conversion until the final output step. This section details this
methodology and introduces the reference implementation, `A15.py`
(Risinger 2024a), used throughout this work for generation, analysis,
and validation.

## A15 Encoding Pipeline: Visual Overview

To provide a clear conceptual overview of the A15 encoding framework,
Figure <a href="#fig-system-diagram" data-reference-type="ref"
data-reference="fig-system-diagram">6</a> illustrates the complete
pipeline from continuous space to deterministic A15 representation. This
visualization complements the detailed technical description of the
multi-stage scaling framework
(<a href="#subsec-scaling-pipeline" data-reference-type="ref+Label"
data-reference="subsec-scaling-pipeline">2.3</a>) and numerical
stability analysis
(<a href="#subsec-stability" data-reference-type="ref+Label"
data-reference="subsec-stability">2.4</a>), offering an accessible entry
point to the framework’s core mechanisms.

<figure id="fig-system-diagram">

<figcaption>Comprehensive system diagram of the A15 encoding framework.
The process begins with continuous 3D coordinates, which are quantized
to the nearest A15 lattice point using either the WPH or TSP
partitioning method. This produces a compact integer A15 identifier
suitable for efficient transmission or storage. For visualization or
further processing, the A15 identifier can be deterministically
reconstructed to floating-point coordinates. When operating in binary or
stable scaling regimes (<span class="math inline">\(\epsilon_\Delta =
0\)</span>), this reconstruction guarantees bit-identical results across
all platforms, creating islands of determinism within the otherwise
approximated float continuum.</figcaption>
</figure>

The pipeline consists of four main stages:

Input: Continuous 3D Space  
The process begins with continuous coordinates from the application’s
native coordinate system, typically represented as floating-point
vectors subject to the approximation challenges detailed in
<a href="#subsec-intro-floats" data-reference-type="ref+Label"
data-reference="subsec-intro-floats">1.1</a>.

Quantization  
Using either the WPH or TSP partitioning method
(<a href="#subsec-intro-partitioning" data-reference-type="ref+Label"
data-reference="subsec-intro-partitioning">1.3</a>), continuous
coordinates are mapped to the nearest A15 lattice point. The
computational efficiency of this step varies based on implementation
strategy.

Compact A15 Identifier  
The output of quantization is an integer identifier that uniquely
specifies the corresponding A15 lattice point. This compact
representation typically requires significantly less storage than the
original floating-point coordinates.

Reconstruction  
When needed for visualization or processing, the A15 identifier can be
deterministically reconstructed to floating-point coordinates. When
operating in stable scaling regimes ($`\epsilon_\Delta = 0`$,
<a href="#subsubsec-stability-regimes" data-reference-type="ref+Label"
data-reference="subsubsec-stability-regimes">2.4.4</a>), this
reconstruction guarantees bit-identical results across all platforms and
environments.

<figure id="fig-stability-analysis">

<figcaption>Analysis of stability regimes based on the relationship
between output scale factor <span
class="math inline">\(\epsilon_\delta\)</span> and inherent base scale
<span class="math inline">\(\epsilon_N\)</span>. The binary regime
(<span class="math inline">\(\epsilon_\delta = \epsilon_N\)</span>)
produces a specific pattern of discrete histogram peaks corresponding to
the intrinsic geometric characteristics of the A15 structure. The stable
regime (<span class="math inline">\(\epsilon_\delta = m \cdot
\epsilon_N, m &gt; 1\)</span>) shows the same pattern plus additional
peaks from the integer multiple scaling. The unstable regime (<span
class="math inline">\(\epsilon_\delta \neq m \cdot \epsilon_N\)</span>)
exhibits a gapped or sparse distribution, indicating approximation
errors. Operating in either binary or stable regimes (<span
class="math inline">\(\epsilon_\Delta = 0\)</span>) is essential for
guaranteeing bit-identical results across platforms.</figcaption>
</figure>

<figure id="fig-hierarchical-representation">

<figcaption>Hierarchical A15 representation using multi-scale relative
addressing. The framework naturally supports efficient encoding across
different spatial scales, from global positioning down to fine-grained
object details. Each level uses an appropriate A15 scale (<span
class="math inline">\(\epsilon_\delta\)</span>) with coordinates defined
relative to its parent context. This enables compact representation
while maintaining precision where needed, with bit allocation tailored
to the requirements of each hierarchical level.</figcaption>
</figure>

The framework’s key innovation lies in the carefully designed
multi-stage scaling pipeline
(<a href="#subsec-scaling-pipeline" data-reference-type="ref+Label"
data-reference="subsec-scaling-pipeline">2.3</a>) that operates within
an internal integer representation during geometric construction. This
approach, coupled with the rigorous stability analysis
(<a href="#subsec-stability" data-reference-type="ref+Label"
data-reference="subsec-stability">2.4</a>), ensures that the A15 lattice
points themselves can be represented exactly in binary floating-point
format, creating islands of perfect determinism within the otherwise
approximated float continuum.

Figure <a href="#fig-stability-analysis" data-reference-type="ref"
data-reference="fig-stability-analysis">7</a> illustrates the three
stability regimes defined by the relationship between output scale
factor $`\epsilon_\delta`$ and inherent base scale $`\epsilon_N`$,
alongside their characteristic histogram patterns. The binary regime
($`\epsilon_\delta = \epsilon_N`$) produces a specific pattern of
discrete peaks corresponding to the inherent geometric characteristics
of the A15 structure. The stable regime
($`\epsilon_\delta = m \cdot \epsilon_N, m > 1`$) shows the same pattern
plus additional peaks, while the unstable regime
($`\epsilon_\delta \neq m \cdot \epsilon_N`$) exhibits a gapped or
sparse distribution indicative of approximation errors.

The hierarchical representation capability shown in
Figure <a href="#fig-hierarchical-representation" data-reference-type="ref"
data-reference="fig-hierarchical-representation">8</a> further extends
the framework’s utility, enabling efficient encoding across different
spatial scales with appropriate precision for each level. This
multi-scale approach allows applications to balance compact global
positioning with high-precision local details, all within the same
coherent A15 framework.

Together, these visualizations provide a comprehensive overview of the
A15 encoding framework, illustrating its pipeline structure, stability
characteristics, and hierarchical representation capabilities. They
serve as a visual complement to the detailed technical descriptions
throughout the manuscript, enhancing accessibility without sacrificing
precision.

## Core Principle: Internal Integer Representation

Applying a discrete lattice like *A15* to represent continuous space
necessitates a rigorous scaling framework. Simply using floating-point
coordinates throughout the geometric construction risks introducing the
very numerical inconsistencies the framework aims to eliminate.
Therefore, the methodology prioritizes calculations within a
well-defined internal integer coordinate system, delaying the mapping to
inexact output formats until the final step. This preserves the precise
geometric relationships inherent in the *A15* structure and forms the
bedrock for numerical stability.

The essential insight is that while floating-point values inherently
approximate most real-world coordinates, integers can exactly represent
discrete positions within a lattice. By constructing the *A15*
structure’s vertices, edges, and cells using only integer coordinates
initially, the framework maintains perfect internal consistency. This
integer-first approach enables deterministic construction of the spatial
partitioning before any potential floating-point approximation occurs
during the final output scaling.

Crucially, this reliance on an internal integer foundation, particularly
when an application deliberately aligns its native coordinate system and
units with a chosen stable A15 scale ($`\epsilon_\Delta = 0`$,
<a href="#subsubsec-stability-diff" data-reference-type="ref+Label"
data-reference="subsubsec-stability-diff">2.4.3</a>), opens pathways for
highly efficient quantization. In such aligned scenarios, mapping
continuous or application-native coordinates to A15 identifiers may
simplify significantly, potentially reducing complex geometric searches
to fast integer arithmetic operations like truncation or bit shifts,
thereby mitigating concerns about quantization overhead through careful
co-design.

## Multi-Stage Scaling Pipeline

The conversion from the abstract *A15* definition to concrete
coordinates involves a systematic, multi-stage pipeline implemented
within `A15.py`. This process uses several interacting parameters and
fixed factors to progressively build the structure within the internal
integer system before final output scaling. The key stages are detailed
below, and the parameters are summarized in
<a href="#tab-scaling-params" data-reference-type="ref+Label"
data-reference="tab-scaling-params">4</a>.

### Primitive Definition and Resolution (`prescale`)

The pipeline originates with the definition of the fundamental geometric
units associated with the *A15* basis sites—typically the pyritohedra
and tetradecahedra related to the Weaire–Phelan Honeycomb method
(<a href="#subsec-intro-partitioning" data-reference-type="ref+Label"
data-reference="subsec-intro-partitioning">1.3</a>), or the cubic blocks
for the Tetrastix Prism method. Functions within `A15.py` (e.g.,
`pyritohedron()`, `tetradecahedra()`) apply an internal integer
multiplier, the `prescale` parameter, to the canonical *dimensionless
fractional* coordinates (like $`0, 1/4, 1/2`$) defining the vertices of
these shapes. This crucial first step converts the fractional values
into *primitive-relative integer coordinates*, establishing the minimum
resolution necessary to accurately represent the individual geometric
primitives without internal approximation. Default `prescale` values in
`A15.py` (typically 20 for standard WPH-derived shapes, 24 for
TSP-derived shapes via `-stix`) are chosen specifically to ensure these
base vertices land precisely on an integer grid relative to the shape’s
center.

### Lattice Placement (Factor 24)

The `lattice()` function in `A15.py` replicates these integer-vertex
polyhedra periodically across 3D space according to the $`Pm\bar{3}n`$
symmetry rules. It determines the positions for the *centers* of these
polyhedra based on *integer lattice vectors* `xyz` relative to an origin
offset `o`. A key element in this placement is the fixed **Lattice
Spacing Factor of 24 units**. Within the `lattice()` function, the
calculation `(xyz + o) * 24` multiplies the integer lattice vector by
this factor. This defines the spacing between the nominal lattice points
(i.e., the centers of the polyhedra) *measured in the integer units
established by `prescale`*. This factor ensures the correct relative
positioning of the repeating units within the overall *A15* lattice
framework.

### Basis Accommodation (Factor 96 Baseline)

The *A15* structure’s complexity arises from its multi-atom basis; it
contains distinct crystallographic sites (Wyckoff 2a and 6d) not just at
the nominal lattice points defined above, but also at specific
fractional offsets, notably involving $`1/4`$ for the 6d sites. To place
*all* required sites onto a single, consistent internal integer grid
requires a resolution finer than the 24-unit spacing factor alone
provides. The smallest denominator involved is 4 (from $`1/4`$).
Therefore, the internal integer grid must be conceptually scaled such
that one unit step along a primitive lattice vector corresponds to
$`4 \times 24 = 96`$ internal integer units along each principal axis.

This derived **Effective Lattice Unit Factor of 96** represents the
fundamental period or effective unit dimension of the comprehensive
internal integer grid required to represent the *complete* *A15*
structure (including its basis) without loss of precision due to these
fractional offsets. This factor of 96 establishes the crucial **baseline
dimension** relative to which the structure’s inherent numerical
precision requirements ($`\epsilon_N`$, see
<a href="#subsubsec-stability-epsilon-n" data-reference-type="ref+Label"
data-reference="subsubsec-stability-epsilon-n">2.4.2</a>) are determined
and against which the final output scale ($`\epsilon_\delta`$) is
compared for stability analysis.

### Optional Binary Rescaling (`rescale`)

Before the lattice generation stage fully populates the internal integer
coordinates, an optional power-of-two scaling factor, `rescale`, can be
applied. This factor is specified via the `-rescale` command-line flag
or implicitly through ‘+/-‘ suffixes appended to shape names in the
input (e.g., `pyritohedra++` implies a rescale factor of $`2^2=4`$ in
`A15.py`). This `rescale` factor multiplies the effective size
represented by the internal integer coordinates relative to the 96-unit
baseline established above. It allows adjustments to the overall size or
relative proportions of different generated components *before* the
final output scaling step, operating purely in powers of two within the
integer domain. Importantly, the choice of `rescale` influences the
resulting inherent base scale $`\epsilon_N`$ of the internal geometry.

### Resultant Internal Integer Coordinates

The cumulative effect of the initial `prescale`
(<a href="#subsubsec-scaling-prescale" data-reference-type="ref+Label"
data-reference="subsubsec-scaling-prescale">2.3.1</a>), the lattice
placement logic
(<a href="#subsubsec-scaling-lattice" data-reference-type="ref+Label"
data-reference="subsubsec-scaling-lattice">2.3.2</a>, establishing the
96-unit effective baseline
<a href="#subsubsec-scaling-baseline" data-reference-type="ref+Label"
data-reference="subsubsec-scaling-baseline">2.3.3</a>), and any applied
`rescale` factor
(<a href="#subsubsec-scaling-rescale" data-reference-type="ref+Label"
data-reference="subsubsec-scaling-rescale">2.3.4</a>) collectively
defines the complete set of vertex coordinates for the entire generated
structure. These coordinates exist within a consistent, potentially
large-valued, *internal integer* system referenced against the 96-unit
effective lattice baseline. By prioritizing integer arithmetic
throughout these construction stages, `A15.py` preserves the precise
geometric relationships inherent in the *A15* structure, deferring any
floating-point approximation until the very last output step. This
internal integer representation forms the bedrock for the numerical
stability analysis detailed next.

Furthermore, the explicit control over scale inherent in this pipeline
naturally supports hierarchical or multi-scale representations via
relative addressing. This allows, for instance, coarse-grained
coordinates identifying an entity’s global position to coexist
efficiently with fine-grained coordinates defining its internal details
(e.g., limb articulations relative to a centroid), each utilizing an
appropriate A15 scale ($`\epsilon_\delta`$) within the same unified
structural logic, optimizing both data density and precision.

<div id="tab-scaling-params">

| **Parameter/Factor** | **Role and Implementation in `A15.py`** |
|:---|:---|
| `prescale` | Internal integer multiplier applied within shape functions (e.g., `pyritohedron()`). Converts base fractional coordinates (0, 1/4, 1/2, etc.) to integer vertices relative to shape center. Establishes primitive resolution (Defaults: 20 WPH, 24 TSP). |
| Lattice Spacing Factor (24) | Fixed factor used in `lattice()`. Multiplies integer lattice vectors `xyz` to determine nominal lattice point positions relative to origin offset `o`, measured in units defined by `prescale`. |
| Effective Lattice Unit Factor (96) | Derived factor ($`4 \times 24 = 96`$). Represents the fundamental period of the internal integer grid along a principal axis needed to accommodate *A15* basis site offsets (like 1/4) precisely. **This is the internal baseline dimension** relative to which inherent precision $`\epsilon_N`$ and output scale $`\epsilon_\delta`$ are compared. |
| `rescale` | Optional pre-generator power-of-two scaling (via `-rescale` flag or `+/-` suffixes). Multiplies internal integer coordinates, adjusting size relative to the 96-unit baseline *before* final output scaling. Affects the resulting $`\epsilon_N`$. |
| `scale` ($`\epsilon_\delta`$) | Final global scaling factor (`-scale=<value>`) applied by `figure()` after internal integers are generated. Maps internal integers (relative to the 96-unit baseline) to output coordinates (typically float). Sets real-world size and determines numerical stability regime ($`\epsilon_\Delta`$) by comparison with $`\epsilon_N`$. Accepts various formats (integer, fraction, power, float). |
| `n` | Controls lattice extent (`-n=<value>`). Integer specifies cuboid dimensions (number of lattice cells defined by the 96-unit baseline along axes); float specifies a spherical radius cutoff based on lattice vector magnitude relative to the origin. |

Summary of Key Scaling Parameters and Factors in `A15.py`.

</div>

## Numerical Stability: Regimes and Validation

Achieving deterministic spatial representation—ensuring that identical
logical operations yield bit-identical results across different
systems—is a primary motivation for this framework. This guarantee
hinges critically on the relationship between the scale chosen for the
final output coordinates ($`\epsilon_\delta`$) and the inherent
precision requirements ($`\epsilon_N`$) of the underlying *A15*
geometry, as captured by the internal integer representation
(<a href="#subsec-scaling-framework" data-reference-type="ref+Label"
data-reference="subsec-scaling-framework">2.2</a>). The framework
defines specific, verifiable conditions under which the mapping from the
precise internal structure to the output coordinate system can be
performed without introducing floating-point approximation errors
relative to that internal structure.

### Global Output Scaling Factor ($`\epsilon_\delta`$)

The primary user control over the final representation’s physical scale
or resolution is the `scale` parameter, specified via the
`-scale=<value>` command-line option in `A15.py`. We denote this crucial
factor as $`\epsilon_\delta`$. It is applied globally by the `figure()`
function *after* all internal integer coordinates (relative to the
96-unit baseline) have been generated. This $`\epsilon_\delta`$ defines
the mapping from the dimensionless internal integer system to the output
coordinate system (typically floating-point). Therefore,
$`\epsilon_\delta`$ effectively sets the real-world size represented by
one unit of the internal 96-unit baseline dimension, and its specific
value is the key determinant of the resulting numerical stability
regime.

It is important to emphasize that $`\epsilon_\delta`$ does not eliminate
the inherent approximation nature of floating-point numbers—most
real-world coordinates will still require approximation when represented
as binary floats. Rather, it establishes a carefully chosen scale at
which the *specific discrete points of the *A15* lattice* can be exactly
represented without such approximation, creating "islands of
determinism" within the otherwise approximated float continuum.

### Inherent Base Scale ($`\epsilon_N`$)

The internal integer geometry produced by the construction process
(<a href="#subsec-scaling-pipeline" data-reference-type="ref+Label"
data-reference="subsec-scaling-pipeline">2.3</a>) possesses an intrinsic
precision limit relative to the 96-unit baseline. This limit arises from
the specific combination of the *A15* basis site coordinates (involving
factors of $`1/4`$), the chosen `prescale` value, and any applied
`rescale` factor. We define $`\epsilon_N`$ as the **inherent base
scale** required to represent this specific internal geometry exactly
when mapped to a base-2 representation. Conceptually, $`\epsilon_N`$ is
the finest scaling factor, expressible as a power of two ($`1/2^N`$ for
some integer $`N`$), that allows all generated internal integer vertex
coordinates to be represented perfectly as rational numbers without
approximation when measured against the 96-unit baseline. It captures
the structure’s intrinsic geometric precision limit within the binary
system used by floating-point numbers. The `A15.py` script
computationally infers $`\epsilon_N`$ by analyzing the power-of-two
denominators required for the exact rational representation of the
generated internal integer geometry relative to the 96-unit baseline.

### Stability Condition and Difference ($`\epsilon_\Delta`$)

Numerical stability is achieved if, and only if, the chosen global
output scale ($`\epsilon_\delta`$) is commensurate with the inherent
base scale ($`\epsilon_N`$). Specifically, the mapping from the internal
integer representation to the output coordinate system is guaranteed to
be exact (free from representation error relative to the internal grid)
if $`\epsilon_\delta`$ is a positive integer multiple ($`m`$) of
$`\epsilon_N`$. The `A15.py` script quantifies this relationship by
calculating the **stability difference**, $`\epsilon_\Delta`$.
Conceptually, $`\epsilon_\Delta`$ measures the mismatch or residual
error when checking if $`\epsilon_\delta`$ aligns perfectly with the
grid defined by $`\epsilon_N`$:
``` math
\label{eq-stability-condition}
\epsilon_\Delta = 0 \quad \iff \quad \epsilon_\delta = m \cdot \epsilon_N \quad \text{for some integer } m \ge 1
```
If this condition holds ($`\epsilon_\Delta = 0`$), the framework
operates in a stable regime. Otherwise ($`\epsilon_\Delta \neq 0`$), the
scaling is unstable, and approximation errors are necessarily introduced
relative to the internal structure.

### Stability Regimes

The stability condition leads to three distinct operational regimes:

<div class="description">

This optimal regime occurs when the chosen output scale precisely
matches the minimum required inherent base scale
($`\epsilon_\delta = \epsilon_N`$). All internal integer coordinates map
directly and exactly onto the binary floating-point grid defined by
$`\epsilon_N = 1/2^N`$ without any approximation relative to the
internal structure.
``` math
\label{eq-scaling-binary}
    \epsilon_\delta = \epsilon_N = \frac{1}{2^N} \quad \implies \quad \epsilon_\Delta = 0
```
This represents the most compact scale that allows exact representation
and guarantees determinism. The validation histogram
(<a href="#fig-hist" data-reference-type="ref+Label"
data-reference="fig-hist">9</a>, left) shows a specific pattern of
discrete peaks at the exponents corresponding to the inherent A15
geometry.

This regime occurs when the output scale is an exact positive integer
multiple ($`m`$) of the inherent base scale
($`\epsilon_\delta = m \cdot \epsilon_N`$). All internal coordinates
still map exactly to representable binary floating-point values without
approximation error relative to this scaled grid.
``` math
\label{eq-scaling-stable}
    \epsilon_\delta = m \cdot \epsilon_N = \frac{m}{2^N}, \quad m \in \mathbb{Z}^+, m > 1 \quad \implies \quad \epsilon_\Delta = 0
```
This regime maintains perfect representability and determinism while
providing flexibility in choosing the overall physical scale. The
smallest resolvable difference (Unit of Least Precision,
<a href="#subsubsec-notes-ulp" data-reference-type="ref+Label"
data-reference="subsubsec-notes-ulp">6.4.5</a>) corresponds to $`m`$
units at the $`\epsilon_N`$ scale, effectively making the ULP equal to
$`\epsilon_\delta`$. The validation histogram
(<a href="#fig-hist" data-reference-type="ref+Label"
data-reference="fig-hist">9</a>, middle) displays the same
characteristic pattern of peaks as the binary regime, plus additional
peaks reflecting the integer scaling factor $`m`$.

This regime occurs whenever the chosen output scale $`\epsilon_\delta`$
is *not* a positive integer multiple of the inherent base scale
$`\epsilon_N`$. Under these conditions, the precise internal integer
geometry cannot be perfectly represented on the binary floating-point
grid implied by $`\epsilon_\delta`$. Rounding errors are necessarily
introduced during the final scaling from internal integers to output
floats.
``` math
\label{eq-scaling-unstable}
    \epsilon_\delta \neq m \cdot \epsilon_N \quad \text{for any} \quad m \in \mathbb{Z}^+ \quad \implies \quad \epsilon_\Delta \neq 0
```
Utilizing unstable scales fundamentally compromises the core benefit of
the framework regarding determinism. It introduces representation
errors, leading to subtle geometric inconsistencies, non-deterministic
outcomes in sensitive calculations, and significant challenges in
achieving reliable state synchronization. The validation histogram
(<a href="#fig-hist" data-reference-type="ref+Label"
data-reference="fig-hist">9</a>, right; also
<a href="#fig-main" data-reference-type="ref+Label"
data-reference="fig-main">10</a>) clearly reveals this instability
through a broad, sparse, or gapped distribution of denominator
exponents.

</div>

Therefore, operating exclusively within the **Binary** or **Stable**
scaling regimes ($`\epsilon_\Delta = 0`$) is essential for numerical
consistency, reproducibility, and guaranteed deterministic behavior.

## Validation via `A15.py` Histogram Analysis

The `A15.py` script provides not only the means to generate structures
based on this framework but also includes a direct, quantitative method
for validating the numerical stability regime resulting from any chosen
set of parameters, particularly the global output scale
$`\epsilon_\delta`$. This validation capability is accessed via the
`-bars` command-line option and provides empirical confirmation of the
stability concepts discussed above
(<a href="#subsec-stability" data-reference-type="ref+Label"
data-reference="subsec-stability">2.4</a>). When enabled, the script
performs the following analysis within its `figure()` function:

1.  It iterates through the components (x, y, z) of output
    floating-point coordinates at the selected `scale` factor
    ($`\epsilon_\delta`$). For each component, it obtains the exact
    rational representation $`(k, d)`$ using Python’s built-in
    `float.as_integer_ratio()` method (Python Software Foundation 2025).
    This captures the precise value representable by the float variable,
    including any approximation introduced if the scaling was unstable
    relative to the internal grid
    (<a href="#subsubsec-scaling-internal" data-reference-type="ref+Label"
    data-reference="subsubsec-scaling-internal">2.3.5</a>).

2.  It then analyzes the denominator $`d`$ of this exact rational
    representation. Since floats operates within a base-2 context, the
    denominator $`d`$ always simplifies to $`2^n`$ for some integer
    exponent $`n`$. The script extracts this effective exponent
    $`n = \log_2 d`$, which conceptually represents $`-\log_2`$ of the
    fractional part’s required precision at the output scale
    $`\epsilon_\delta`$.

3.  Finally, it visualizes the distribution of these calculated
    exponents $`n`$ across all analyzed vertex components as a histogram
    (<a href="#fig-hist" data-reference-type="ref+Label"
    data-reference="fig-hist">9</a>). This histogram provides immediate
    visual feedback on the numerical stability of the chosen
    configuration.

The shape of this generated histogram directly corresponds to the
theoretical stability regimes
(<a href="#subsubsec-stability-regimes" data-reference-type="ref+Label"
data-reference="subsubsec-stability-regimes">2.4.4</a>): a specific
pattern of discrete peaks signifies Binary scaling; the same pattern
plus additional peaks signifies Stable scaling; and a sparse, broad, or
gapped distribution signifies Unstable scaling, revealing the
introduction of approximation errors relative to the internal geometry.
Furthermore, the script explicitly calculates and displays the stability
difference $`\epsilon_\Delta`$ alongside the histogram
(<a href="#fig-main" data-reference-type="ref+Label"
data-reference="fig-main">10</a>), offering direct numerical
confirmation of the regime. This built-in quantitative validation
provides objective criteria for selecting scaling parameters that
correctly balance spatial resolution, memory efficiency, and the
critical requirement for numerical robustness and determinism.

<figure id="fig-hist">
<p><img src="fig-histb.png" alt="image" /> <img src="fig-hists.png"
alt="image" /> <img src="fig-histu.png" alt="image" /></p>
<figcaption>Example histograms generated by <code>A15.py</code>
(<code>-bars</code>) illustrating numerical stability regimes for
different output scaling factors (<span
class="math inline">\(\epsilon_\delta\)</span>). Left
(<code>fig-histb.png.txt</code>, <span
class="math inline">\(\epsilon_\delta=1/64\)</span>):
<strong>Binary</strong> scale (<span
class="math inline">\(\epsilon_\Delta=0\)</span>) showing a distinctive
pattern of discrete denominator exponents (<span class="math inline">\(n
= \log_2 d\)</span>) reflecting the intrinsic geometric characteristics
of the A15 structure. Middle (<code>fig-hists.png.txt</code>, <span
class="math inline">\(\epsilon_\delta=683/2^{16}\)</span>):
<strong>Stable</strong> scale (<span
class="math inline">\(\epsilon_\Delta=0\)</span>) showing the same
pattern plus additional exponents introduced by the integer multiple
scaling. Right (<code>fig-histu.png.txt</code>, <span
class="math inline">\(\epsilon_\delta=1/96\)</span>):
<strong>Unstable</strong> scale (<span
class="math inline">\(\epsilon_\Delta \neq 0\)</span>) showing a widely
spread distribution with gaps, indicating floating-point approximation
errors introduced during final scaling relative to the internal integer
geometry.</figcaption>
</figure>

## `A15.py` Reference Implementation Overview

The Python script `A15.py` (Risinger 2024a) serves as the reference
implementation, exploratory tool, and validation instrument for the
concepts presented in this research. It requires Python 3 and leverages
standard scientific libraries—NumPy (Harris et al. 2020) for efficient
numerical array operations, SciPy (Virtanen et al. 2020) for geometric
computations (such as convex hull calculations used for visualization
via `scipy.spatial.ConvexHull`), and Matplotlib (Hunter 2007) for
versatile 2D and 3D visualization. The script’s workflow systematically
applies the principles of the *A15* encoding framework, incorporating
the crucial multi-stage scaling logic
(<a href="#subsec-scaling-pipeline" data-reference-type="ref+Label"
data-reference="subsec-scaling-pipeline">2.3</a>) and the quantitative
numerical stability analysis
(<a href="#subsec-stability-validation" data-reference-type="ref+Label"
data-reference="subsec-stability-validation">2.5</a>) described
previously.
<a href="#tab-a15py-workflow" data-reference-type="ref+Label"
data-reference="tab-a15py-workflow">5</a> summarizes the core stages of
its operation, from parameter parsing to final visualization and
validation. This structured process allows `A15.py` to generate,
visualize, and analyze configurations, providing concrete examples and
empirical validation of the framework’s properties, such as the
composite output shown in
<a href="#fig-main" data-reference-type="ref+Label"
data-reference="fig-main">10</a>.

<div id="tab-a15py-workflow">

| **Stage** | **Key Functions & Purpose in `A15.py`** |
|:---|:---|
| Configuration Processing | `flags()`, `configuration()` |
|  | Parses command-line arguments/files (`*.txt`). Interprets parameters (e.g., `scale`, `rescale`, `n`, `prescale`, `stix`, visualization flags like `-edges`, `-faces`, `-bars`). Handles defaults, automatic configurations (`-auto`), inheritance (colon notation). Determines final parameter set for generation and visualization. |
| Shape Definition | `pyritohedron()`, `tetradecahedra()` |
|  | Constructs fundamental geometric units based on parameters. Applies internal `prescale` factor converting base fractional coordinates to integer vertices relative to shape center. Handles local symmetry (handedness). |
| Lattice Generation | `lattice()` |
|  | Replicates shape units across 3D grid based on extent `n` and `at` filter (for basis site mapping). Calculates center positions using integer vectors `xyz`, origin offset `o`, and the fixed `* 24` spacing factor (relative to `prescale`d units), implicitly establishing the 96-unit effective baseline. Yields (vertex_array, config_object) pairs representing internal integer geometry. |
| Scaling & Stability Analysis | `figure()` |
|  | Collects internal integer vertex arrays. Applies global `scale` parameter ($`\epsilon_\delta`$) mapping internal integers (relative to the 96-unit baseline) to output floats. Infers inherent base scale $`\epsilon_N`$. Calculates stability difference $`\epsilon_\Delta`$ (verifying if $`\epsilon_\delta = m \cdot \epsilon_N`$ for $`m \in \mathbb{Z}^+ \implies \epsilon_\Delta=0`$). |
| Visualization & Validation | `figure()` |
|  | Renders 3D geometry via Matplotlib (Hunter 2007) using specified options (`-edges`, `-faces`, etc.). Uses SciPy (Virtanen et al. 2020) for hull calculations if needed. If `-bars` requested, performs stability analysis (via `float.as_integer_ratio()` (Python Software Foundation 2025) denominators) and generates histogram (<a href="#fig-hist" data-reference-type="ref+Label"
data-reference="fig-hist">9</a>), visualizing the stability regime and displaying $`\epsilon_\Delta`$. Outputs to screen (`-interactive` option implies pop-up) or file (`-savefig`). Includes annotations (<a href="#fig-main" data-reference-type="ref+Label"
data-reference="fig-main">10</a>). Relies heavily on NumPy (Harris et al. 2020) throughout. |

Core Workflow Stages in the `A15.py` Implementation.

</div>

<figure id="fig-main">
<img src="fig-main.png" style="width:80.0%" />
<figcaption>Composite visualization generated by <code>A15.py</code>
(<code>fig-main.png.txt</code>) showing pyritohedra and tetradecahedra
components alongside stability analysis. This example uses an
<strong>Unstable</strong> output scale <span
class="math inline">\(\epsilon_\delta = 1/96\)</span>. The calculated
non-zero stability difference (<span
class="math inline">\(\epsilon_\Delta \approx \num{0.0104} \neq
0\)</span>) is explicitly annotated in the histogram sidebar (left),
confirming instability relative to the internal integer geometry. The
histogram visually reflects this instability through its wide, gapped
distribution of floating-point denominator exponents (<span
class="math inline">\(n = \log_2 d\)</span>). Main 3D views use
dimensionless coordinates <span class="math inline">\(N_X, N_Y,
N_Z\)</span> relative to the basic unit width <span
class="math inline">\(N_1\)</span>. For this specific unstable scale,
<span class="math inline">\(N_1\)</span> corresponds to 1.0 mm, derived
from applying <span class="math inline">\(\epsilon_\delta=1/96\)</span>
to the internal 96-unit effective lattice baseline. Counts shown refer
to the number of float components analyzed for the
histogram.</figcaption>
</figure>

# Interpretation, Benefits, and Limitations

The investigation detailed in the preceding sections confirms the
distinctive suitability of the *A15* phase structure, when employed
within the numerically stable scaling framework
(<a href="#sec-framework-design" data-reference-type="ref+Label"
data-reference="sec-framework-design">2</a>), for partitioning
interactive 3D spaces. The confluence of its inherent crystallographic
properties (<a href="#sec-introduction" data-reference-type="ref+Label"
data-reference="sec-introduction">1</a>) and the demonstrable numerical
stability achievable through disciplined scaling
($`\epsilon_\Delta = 0`$,
<a href="#subsubsec-stability-diff" data-reference-type="ref+Label"
data-reference="subsubsec-stability-diff">2.4.3</a>) yields a framework
that is both geometrically sophisticated and computationally practical
for demanding immersive applications. This section interprets these
findings, highlights the quantifiable benefits derived directly from the
framework’s mechanics, and candidly discusses practical limitations and
engineering considerations relevant to its implementation.

## Interpretation and Core Findings

This research establishes *A15* encoding as a robust structural
foundation for coordinated spatial partitioning. The core finding is
that by leveraging the specific crystallographic symmetries and
binary-friendly coordinates of the *A15* structure, and by rigorously
adhering to the multi-stage scaling pipeline culminating in a stable
output scale ($`\epsilon_\Delta = 0`$), it is possible to create a
spatial representation that fundamentally eliminates floating-point
representation errors relative to its own discrete grid. This provides a
pathway to achieving verifiable determinism in spatial computations, a
critical enabler for next-generation networked virtual environments and
related spatial computing tasks.

### Privacy and Security of Tracking Data (PII)

The application of this efficient encoding framework to fine-grained
spatial tracking data, especially full-body kinematics derived from
VR/AR systems, carries significant privacy implications that **must** be
addressed with utmost seriousness. Such detailed movement data
constitutes **Personally Identifiable Information (PII)** and may
qualify as sensitive biometric data under various regulations (e.g.,
GDPR (European Parliament and Council of the European Union 2016), CCPA
(California State Legislature 2018)). Handling this data demands
rigorous privacy safeguards and unwavering ethical considerations as a
**non-negotiable** aspect of implementation. Developers and deployers
**must** integrate robust security measures as a foundational
requirement. This includes, at a minimum:

- Employing strong **end-to-end encryption (E2EE)** for all
  *A15*-encoded coordinate streams and associated tracking data during
  network transmission and persistent storage.

- Strict adherence to **data minimization principles** (collecting only
  the data essential for the application’s functionality).

- Implementing transparent **user consent mechanisms** before any
  tracking begins.

- Establishing clearly defined **data retention and deletion policies**.

- Utilizing effective **anonymization or aggregation strategies**
  whenever full individual fidelity is not strictly required (e.g., for
  analytics or heatmaps).

- Ensuring full compliance with all relevant legal and ethical
  regulations.

This data represents individuals and their behavior; it **must** be
treated with the highest degree of care, respect, security, and
transparency. Failure to do so carries significant legal, ethical, and
reputational risks.

It should be noted that while the A15 encoding framework reduces the raw
size of spatial data, potentially making it more efficient to encrypt
and secure, this efficiency gain does not diminish the fundamental
privacy requirements surrounding such sensitive information.

## Quantifiable Benefits

The *A15* encoding framework, when implemented correctly within stable
scaling regimes, offers several key advantages validated by theoretical
analysis and the `A15.py` (Risinger 2024a) reference implementation:

### Guaranteed Determinism via Stable Scaling

Arguably the most significant contribution stems directly from operating
within the rigorously defined **Binary** or **Stable** scaling regimes
($`\epsilon_\Delta = 0`$,
<a href="#subsubsec-stability-regimes" data-reference-type="ref+Label"
data-reference="subsubsec-stability-regimes">2.4.4</a>). By precisely
aligning the chosen output scale ($`\epsilon_\delta`$) with the
structure’s inherent geometric precision requirements ($`\epsilon_N`$,
relative to the 96-unit baseline derived in
<a href="#subsubsec-scaling-baseline" data-reference-type="ref+Label"
data-reference="subsubsec-scaling-baseline">2.3.3</a>), the framework
guarantees that all *A15* lattice points and derived vertex coordinates
map exactly onto hardware binary floating-point formats *relative to
that chosen stable scale*.

It is important to emphasize that this guarantee applies specifically to
the representation of the discrete A15 lattice points themselves; it
ensures that the quantized spatial framework remains consistent across
systems, but does not eliminate all floating-point variability within
applications using this framework. Operations performed between these
stable grid points (such as physics calculations, interpolation, or
velocity updates) may still use floating-point arithmetic and thus be
subject to traditional floating-point variability unless additional
measures are taken.

Nevertheless, this systematic elimination of representation errors for
the spatial partitioning grid itself provides a foundation for
much-improved consistency across disparate machines, platforms, and even
different compilation environments—a fundamental prerequisite for
reliable state synchronization in networked systems, reproducible
physics simulations, efficient network delta compression strategies,
verifiable event sequences and replays, and ultimately, fair and
trustworthy competitive experiences.

### Memory and Bandwidth Efficiency

Mapping continuous coordinates onto compact, often integer-based, *A15*
identifiers provides substantial memory and bandwidth savings compared
to standard floating-point vector representations
(<a href="#subsec-intro-floats" data-reference-type="ref+Label"
data-reference="subsec-intro-floats">1.1</a>). This advantage is
particularly pronounced when quantizing explicitly defined or bounded
volumes where the full dynamic range and mantissa precision of standard
floats represent unnecessary overhead. For a practical example, consider
establishing a coordinate system using a stable scale where the basic
unit width $`N_1`$ (derived from the 96-unit baseline,
<a href="#subsubsec-scaling-baseline,subsubsec-notes-figures"
data-reference-type="ref+Label"
data-reference="subsubsec-scaling-baseline,subsubsec-notes-figures">[subsubsec-scaling-baseline,subsubsec-notes-figures]</a>)
corresponds to 1.5 mm (achieved with $`\epsilon_\delta=1/64`$). A 48
integer identifier could then be structured as follows: 3 bits can
uniquely address the 8 distinct atomic sites within the A15 conventional
cubic cell (2 at Wyckoff 2a, 6 at Wyckoff 6d); allocating 11 bits for
the vertical axis provides
$`2^{11} \times \SI{1.5}{\milli\meter} \approx \SI{3.07}{\meter}`$ of
range, sufficient vertical headroom for human-scale interaction; the
remaining 34 bits, split evenly (17+17) for the horizontal axes, define
a square area of
$`(2^{17} \times \SI{1.5}{\milli\meter})^2 \approx (\SI{196.6}{\meter})^2`$.
This 48-bit structure, capable of encoding a space several times larger
in area than the largest professional stadium fields with sub-millimeter
intra-cell precision, represents a **50% reduction** compared to the 96
typically required for three standard single-precision (32) floats
(<a href="#eq-efficiency-memory" data-reference-type="ref+Label"
data-reference="eq-efficiency-memory">[eq-efficiency-memory]</a>).
``` math
\label{eq-efficiency-memory}
    \text{Savings} = \frac{(96\,\text{bit} - 48\,\text{bit})}{96\,\text{bit}} \times 100\% = 50\%
```
Furthermore, because many application environments (e.g., smaller
arenas, non-square layouts (Risinger 2024b)) require less range, the bit
allocation can often be reduced further, leading to savings
significantly exceeding 50%. This efficiency translates directly into
reduced memory footprints for spatial data structures and significantly
lower network traffic for coordinate updates. Considering just baseline
kinematic tracking data
(<a href="#eq-bandwidth-baseline" data-reference-type="ref+Label"
data-reference="eq-bandwidth-baseline">[eq-bandwidth-baseline]</a>),
this reduction from approximately 67.2   s<sup>−1</sup> down to
34   s<sup>−1</sup> or less offers significant leverage on aggregate
bandwidth, especially in complex scenarios often exceeding
100   s<sup>−1</sup> per user for kinematics alone.
``` math
\label{eq-bandwidth-baseline}
    \text{Rate}_{\text{3xFloat}} = 20\,\text{jts} \times (3_{\text{pos}} + 4_{\text{quat}})\,\tfrac{\text{floats}}{\text{jt}} \times 4\,\tfrac{\text{bytes}}{\text{float}} \times 15\,\text{Hz} = 8400\,\text{byte/s} \approx 67.2\,\text{kbit/s}
```
The *A15* integer encoding effectively reclaims storage and bandwidth
otherwise consumed by unused float exponent bits and excess significand
precision within suitably bounded contexts. Furthermore, this
compactness and structural definition offer benefits for long-term data
storage and archival, providing a potentially more stable and
interpretable format compared to raw floating-point streams. Practical
applications, such as the related `layoutc` project which encodes
physical layouts using similar principles (Risinger 2024b), demonstrate
the utility of such compact, deterministic encodings.

### Geometric Fidelity and Isotropy

Leveraging the intrinsic crystallographic properties of the *A15*
structure (<a href="#subsec-intro-a15" data-reference-type="ref+Label"
data-reference="subsec-intro-a15">1.2</a>) imparts beneficial geometric
qualities to the spatial representation. Its verified high mean
coordination number (13.5) indicates efficient local packing and dense
connectivity between neighboring regions. The crucial $`T_h`$ point
group symmetry provides a high degree of local isotropy by incorporating
near-icosahedral geometric elements within a globally cubic (and thus
perfectly periodic) framework
(<a href="#subsec-intro-a15" data-reference-type="ref+Label"
data-reference="subsec-intro-a15">1.2</a>). This structural integrity,
especially when combined with the topologically matched
Weaire–Phelan Honeycomb local discretization method
(<a href="#subsec-intro-partitioning" data-reference-type="ref+Label"
data-reference="subsec-intro-partitioning">1.3</a>), fosters a spatially
uniform or “fair” representation, helping to minimize the directional
biases or artifacts that can plague simpler grid-based partitioning
schemes.

### Suitability for Distributed Computing

The inherent regularity, crystallographic symmetries, and predictable
neighborhood topology of the A15 lattice provide a structured substrate
well-suited for spatial domain decomposition. This facilitates the
distribution of spatial computations across parallel architectures,
including clusters, GPUs, or edge networks, potentially simplifying load
balancing and data management compared to adaptive or irregular spatial
structures. The deterministic nature of the grid ensures consistent
partitioning across nodes operating under stable scaling regimes.

### Foundation for Enhanced Interoperability

By establishing a common, mathematically precise, and verifiable spatial
structure, the A15 framework offers a robust foundation for enhanced
interoperability. Diverse applications adhering to the same A15
structure, scale, and orientation conventions could potentially
exchange, reference, or merge spatial data with greater semantic
consistency and reduced ambiguity, fostering cohesion across federated
virtual environments or collaborative platforms.

### Validated Implementation Framework

The accompanying `A15.py` script
(<a href="#subsec-implementation-a15py" data-reference-type="ref+Label"
data-reference="subsec-implementation-a15py">2.6</a>) serves as more
than just a visualization aid; it is a validated reference
implementation and analysis framework. It demonstrates the practical
construction of *A15*-based structures, rigorously implements the
multi-stage scaling logic essential for achieving numerical stability
(<a href="#subsec-scaling-pipeline" data-reference-type="ref+Label"
data-reference="subsec-scaling-pipeline">2.3</a>), and provides the
quantitative histogram analysis
(<a href="#subsec-stability-validation" data-reference-type="ref+Label"
data-reference="subsec-stability-validation">2.5</a>,
<a href="#fig-hist" data-reference-type="ref+Label"
data-reference="fig-hist">9</a>) that empirically confirms the existence
and accessibility of the stable scaling regimes
($`\epsilon_\Delta = 0`$). This ensures the reproducibility of the core
findings regarding numerical stability and offers a concrete, verifiable
starting point for developers seeking to implement or explore the *A15*
encoding framework for their own applications.

## Limitations and Engineering Considerations

Despite its significant advantages for achieving determinism and
efficiency, adopting the *A15*-based partitioning approach involves
several practical limitations and engineering challenges that require
careful consideration during implementation. These represent addressable
design aspects rather than fundamental flaws in the underlying concept:

### Complexity of Arbitrary Rotations

Applying arbitrary rotations relative to the lattice axes directly to
the integer lattice coordinates used in *A15* encoding can be intricate,
particularly when mapping points onto valid *A15* basis sites (Wyckoff
2a or 6d) or navigating complex cell boundaries like those in
the Weaire–Phelan Honeycomb. Correct implementation necessitates
carefully calculated transformations aware of the crystal basis and
space group operations (Aroyo 2016).

In practice, this means that when objects rotate in the virtual space,
the mapping of their vertices or bounding volumes to A15 identifiers
requires additional computational steps beyond simple coordinate
transformation. For arbitrary rotations, the system must typically:

1.  Transform object-local coordinates to world space (standard
    transform)

2.  Apply nearest-neighbor or containment tests to determine which A15
    cell contains each point

3.  Map these points to appropriate A15 lattice identifiers

The complexity depends on the nature of the rotation and the chosen
local discretization method; axis-aligned rotations or operations within
simpler partitioning geometries like the Tetrastix Prism may present
fewer challenges. Regardless of the approach, adopting a canonical
orientation convention
(<a href="#subsubsec-limits-handedness" data-reference-type="ref+Label"
data-reference="subsubsec-limits-handedness">3.3.6</a>) is essential for
consistent interpretation and interoperability.

### Quantization Cost and Design Alignment

A practical consideration is the computational cost of
quantization—mapping arbitrary continuous coordinates to the nearest
discrete *A15* identifier. While general-purpose nearest-neighbor
searches in 3D can be computationally intensive, especially with complex
cell boundaries (e.g., Weaire–Phelan Honeycomb), this cost is highly
dependent on implementation strategy and application alignment.

For arbitrary coordinates in arbitrary orientations relative to the A15
grid, the quantization process could require:

- Point-in-cell tests against multiple candidate cells

- Distance calculations to determine the nearest lattice point

- Potentially complex geometric intersection tests for the
  Weaire–Phelan Honeycomb partitioning method

However, as noted
(<a href="#subsec-scaling-framework" data-reference-type="ref+Label"
data-reference="subsec-scaling-framework">2.2</a>), if an application’s
coordinate system is deliberately aligned with a stable *A15* scale
($`\epsilon_\Delta = 0`$), quantization can potentially become a highly
efficient process dominated by integer arithmetic operations (e.g.,
truncation or bit shifts), making performance manageable through
informed design choices rather than being an inherent bottleneck. This
approach represents a co-design strategy where the application’s spatial
system is developed with A15 quantization in mind from the outset,
rather than being applied as an afterthought.

### Integration with Existing Float-Based Systems

Incorporating the A15 framework into existing engines and applications
that rely heavily on floating-point arithmetic presents integration
challenges. Most modern game engines, simulation environments, and
physics systems operate natively with floating-point coordinates and
transformation matrices. Introducing A15 encoding requires careful
attention to the boundaries between:

- Float-based internal physics and rendering systems

- A15-encoded positions for storage and network transmission

- Coordinate and state transformations between these representations

Practical implementation strategies might include:

- Using A15 encoding only for network transmission and state
  synchronization

- Maintaining dual representations (float for local processing, A15 for
  verifiable operations)

- Developing middleware translation layers between engine components

- Implementing custom physics solvers that operate directly on the A15
  grid

Each approach involves trade-offs between implementation complexity,
performance, and the degree of determinism achieved across the entire
application pipeline.

### Scalability and Spatial Federation

The framework naturally defines discrete, structured zones based on the
generated *A15* lattice extents (controlled by the `n` parameter in
`A15.py`). Seamlessly extending this model to create massive, open-world
environments requires robust mechanisms for managing transitions and
maintaining coordinate and state consistency across the boundaries
between independent *A15* zones. Addressing this likely involves
developing standardized inter-zone boundary protocols for coordinate
transformations (potentially involving scale changes), object state
hand-offs, and perhaps authority transfer between simulation domains.
Level-of-Detail (LOD) systems (Luebke et al. 2002) utilizing
multi-resolution *A15* grids (employing coarser quantization for distant
or less critical zones, potentially via relative addressing) might also
play a role in managing complexity at large scales.

### Encoding Efficiency for Irregular Volumes

Using lattice-aligned cuboidal extents for defining *A15* zones, while
straightforward to implement via the `n` parameter, can lead to
inefficient use of the addressable integer coordinate space when
representing environments with highly irregular external boundaries,
complex internal terrain features (like mountains or caves), or
significant voids (e.g., the interior of large, non-rectangular
buildings). This may result in significant portions of the allocated
integer coordinate range being unused or unreachable within the playable
or interactive space. Potential mitigations could include implementing
secondary encoding or metadata schemes (e.g., run-length encoding of
valid zones along axes, hierarchical spatial masks, sparse data
structures like sparse voxel octrees adapted to the A15 lattice),
exploring hybrid approaches that combine the regular *A15* grid with
adaptive structures primarily for managing boundary details or sparse
areas (though this might reintroduce some complexity), or effectively
managing precision and sparsity through hierarchical, multi-scale A15
grids utilizing relative addressing
(<a href="#subsec-scaling-pipeline" data-reference-type="ref+Label"
data-reference="subsec-scaling-pipeline">2.3</a>).

### Requirement for Handedness Convention

While the overall $`Pm\bar{3}n`$ space group is centrosymmetric
(achiral), the specific arrangement of atoms in the *A15* structure’s
basis results in alternating left- and right-handed local coordination
environments around the 6d sites
(<a href="#subsec-intro-a15" data-reference-type="ref+Label"
data-reference="subsec-intro-a15">1.2</a>). This local chirality is
implemented deterministically based on lattice position within the
`A15.py` reference code, and is consistent with a standard
**right-handed coordinate system** interpretation. However, for
interoperability in any practical application, particularly networked
ones, all participating systems **must establish and strictly adhere to
a shared global orientation convention**, independent of any specific
implementation’s internal standard. This convention dictates how local
chiralities and global axes are interpreted, represented, and
transformed, ensuring consistency between different client
implementations and, crucially, when interfacing with host environments
or game engines that may use different native coordinate system
handedness (e.g., left-handed systems common in Unity (Unity
Technologies 2024) and Unreal Engine (Epic Games 2023) versus
right-handed systems standard in physics and mathematics). Without such
a shared convention, mirrored or incorrectly oriented geometry could
easily result from exchanging *A15*-encoded coordinates.

# Immediate Potential and Future Prospects

The *A15*-based spatial partitioning and encoding framework, validated
for its numerical stability and efficiency
(<a href="#sec-discussion" data-reference-type="ref+Label"
data-reference="sec-discussion">3</a>), offers immediate potential for
enhancing current virtual environments and provides a solid foundation
for future advancements in spatial computing. Its inherent structural
and numerical properties lend themselves to broad applicability, while
also opening intriguing avenues for further research and development
aimed at refining and extending its capabilities.

## Immediate Practical Applications

The practical adoption of *A15* partitioning across diverse virtual
platforms is facilitated by several key features, yielding tangible
benefits for various applications available today:

### Cross-Platform Compatibility and Interoperability

The underlying $`Pm\bar{3}n`$ space group of the *A15* structure is
centrosymmetric, lacking inherent global chirality (Aroyo 2016). This
structural property means the fundamental lattice encoding can be
consistently represented within both **left-handed coordinate systems**
(common in game engines like Unity (Unity Technologies 2024) and Unreal
Engine (Epic Games 2023)) and **right-handed systems** (standard in
mathematics and physics) through straightforward external affine
transformations. While this simplifies cross-platform integration,
achieving unambiguous interoperability absolutely requires establishing
and adhering to a shared global orientation convention
(<a href="#subsubsec-limits-handedness" data-reference-type="ref+Label"
data-reference="subsubsec-limits-handedness">3.3.6</a>).

In practical terms, applications implementing the A15 framework should:

- Document the specific handedness convention used (explicitly
  right-handed or left-handed)

- Define the precise alignment of A15 crystallographic axes with
  world-space axes

- Specify the origin point of the A15 grid relative to world space

- Include these specifications in any data exchange protocols or
  standards

Crucially, the shared mathematical foundation provides a pathway toward
application-agnostic spatial understanding, enabling potentially more
seamless data exchange between disparate systems built upon the same A15
conventions.

### Integration Strategies with Existing Engines

For practical implementation in current game engines and simulation
environments, several integration approaches can leverage A15 encoding
while minimizing disruption to existing pipelines:

Dual-Representation Strategy:  
This approach maintains two parallel coordinate representations:

- Native engine coordinates (typically floats) for internal processing,
  physics, and rendering

- A15 identifiers for network transmission, authoritative state storage,
  and verification

The synchronization between these representations occurs at well-defined
boundaries:

1.  Convert native coordinates to A15 identifiers before network
    transmission or state archival

2.  Convert received A15 identifiers back to native coordinates for
    local processing

3.  For critical verifications, operations follow the same path through
    A15 encoding/decoding

This strategy minimizes integration complexity while still gaining the
bandwidth efficiency and deterministic verification benefits of A15
encoding.

Middleware Strategy:  
A specialized middleware layer can handle the translation between native
engine coordinates and A15 identifiers transparently:

- Intercept network communication and convert coordinates on the fly

- Provide verification services for authoritative operations

- Implement comparison operators that account for A15 discretization

This approach allows for incremental adoption without significant engine
modifications, though it may introduce additional computational overhead
at translation boundaries.

Native Implementation Strategy:  
For applications requiring maximal determinism and efficiency, a more
comprehensive approach involves:

- Implementing custom spatial data structures directly using A15
  identifiers

- Developing physics solvers that operate on the A15 grid with
  deterministic integer arithmetic

- Constructing rendering pipelines aware of A15 discretization

This approach requires significant engineering investment but can yield
systems with strong determinism guarantees throughout the entire
application pipeline, not just at network boundaries.

### Alignment with Global Measurement Standards

The explicit output scaling mechanism ($`\epsilon_\delta`$, controlled
via the `-scale` parameter in `A15.py`) allows the dimensionless
internal *A15* lattice coordinates (relative to the 96-unit effective
baseline,
<a href="#subsubsec-scaling-baseline" data-reference-type="ref+Label"
data-reference="subsubsec-scaling-baseline">2.3.3</a>) to map directly
and predictably onto standard physical units, such as SI meters or
Imperial feet. Critically, selecting a scale factor within the Binary or
Stable regimes ($`\epsilon_\Delta=0`$,
<a href="#subsubsec-stability-regimes" data-reference-type="ref+Label"
data-reference="subsubsec-stability-regimes">2.4.4</a>) ensures this
mapping is exact relative to the framework’s chosen resolution. This
facilitates interoperability not only between virtual systems but also
with real-world measurements, enhancing user comprehension and grounding
virtual spaces in familiar metrics. For instance, the recommended
baseline scale of $`\epsilon_\delta = 2^{-6}`$ (`-scale=1/64`) provides
robust sub-millimeter precision while operating within a numerically
stable regime, yielding a basic unit width ($`N_1`$) of 1.5 mm
(<a href="#subsubsec-notes-figures" data-reference-type="ref+Label"
data-reference="subsubsec-notes-figures">6.4.7</a>), suitable for many
human-scale interactions.

### Hierarchical Representation via Relative Addressing

The framework naturally enables multi-scale representations through
relative addressing
(<a href="#subsec-scaling-pipeline" data-reference-type="ref+Label"
data-reference="subsec-scaling-pipeline">2.3</a>). Applications can
utilize this capability to encode coarse global positions efficiently
while representing fine-grained local details (like avatar kinematics or
intricate environmental features) with high precision relative to a
parent coordinate. For example, an avatar’s core position could be
stored on a moderate-resolution A15 grid, while its complex joint
movements are encoded on a much finer A15 grid defined locally relative
to that core position. This approach optimizes the balance between data
size, spatial range, and the level of detail required for different
components of a scene.

### Example Use Cases

These features position *A15* partitioning as a robust foundation
suitable for several demanding applications:

- **Competitive VR Esports:** Where deterministic replays, fairness
  verification, and efficient network synchronization are critical. The
  framework’s guaranteed consistency across different client hardware
  and reduced bandwidth requirements are particularly valuable.

- **Industrial Digital Twins:** Requiring precise spatial fidelity and
  alignment with real-world measurements. The ability to precisely map
  A15 coordinates to physical units enables accurate synchronization
  between physical and virtual counterparts.

- **Distributed Physics Simulations:** Leveraging the regular lattice
  structure for domain decomposition across computing nodes while
  maintaining state consistency. The deterministic representation helps
  ensure reproducible results across different simulation environments.

- **Collaborative Mixed Reality:** Where multiple participants with
  diverse devices need to maintain a shared spatial understanding. The
  compact representation and cross-platform compatibility facilitate
  this shared reference frame.

- **Procedural Content Generation:** Benefiting from deterministic
  spatial addressing to ensure consistent generation results. The
  hierarchical capabilities allow for efficient representation of
  multi-scale structures.

The related `layoutc` project (Risinger 2024b), which applies similar
compact, deterministic encoding principles to physical layouts, provides
a concrete example of these concepts in practice.

## Future Research Directions

The *A15* partitioning framework, while demonstrating immediate utility,
also catalyzes numerous avenues for future research aimed at further
enhancing immersive experiences and extending the capabilities of
spatial computing:

### Performance Characterization and Optimization

A critical direction for future work involves comprehensive empirical
performance analysis of A15 implementations across various hardware and
software environments:

- Measuring quantization costs for different cell geometries (WPH vs.
  TSP) and comparing optimization strategies

- Benchmarking bandwidth consumption in realistic multi-user scenarios
  with full kinematic tracking

- Exploring hardware-accelerated implementations, potentially leveraging
  GPU parallelism or specialized SIMD instructions for efficient
  nearest-neighbor mapping

- Developing optimized data structures specifically tailored to
  A15-based spatial partitioning

Such empirical validation would provide crucial guidance for
implementers regarding performance expectations and optimization
strategies across different hardware targets and use cases.

### Advanced Complex Geometry Representation

Developing robust and efficient methods for representing non-cuboid
volumes or complex geometric features within the A15 framework
(<a href="#subsubsec-limits-complex" data-reference-type="ref+Label"
data-reference="subsubsec-limits-complex">3.3.5</a>) remains a key area
for practical improvement. Research could explore:

- Hybrid approaches combining the base A15 grid with techniques like
  sparse octrees or boundary representations for detail

- Constructive solid geometry (CSG) operations defined relative to A15
  cells

- Multi-resolution *A15* representations (spatial Level-of-Detail, LOD
  (Luebke et al. 2002)) using the framework’s inherent support for
  relative addressing
  (<a href="#subsubsec-apps-relative" data-reference-type="ref+Label"
  data-reference="subsubsec-apps-relative">4.1.4</a>)

- Efficient encoding schemes for representing partially filled or
  complex boundary regions within a coarser A15 grid

These approaches would address the current limitation of efficiently
representing irregular volumes or sparse scenes while maintaining the
core benefits of the A15 framework.

### Deterministic Physics on the A15 Grid

Extending the determinism guarantees beyond just spatial representation
to the actual simulation dynamics represents a promising frontier:

- Developing physics solvers that operate directly on A15 grid points
  using fixed-point or integer arithmetic

- Creating collision detection algorithms specifically optimized for A15
  cell geometries

- Defining temporal integration schemes that maintain determinism across
  hardware platforms

- Exploring trade-offs between simulation fidelity and guaranteed
  reproducibility

Such research could lead to fully deterministic simulation environments
where not just the spatial coordinates but all dynamic interactions
operate with bit-exact consistency across systems.

### Framework for Metaverse Interoperability

Standardized protocols built upon *A15* encoding (or similar
deterministic lattice-based systems) could define universal mechanisms
for:

- Agent state representation across diverse virtual environments

- Coordinate referencing between independently developed worlds

- Interaction semantics and object transformation between spaces

- Capability negotiation and feature discovery between A15-compatible
  systems

This research direction could establish a foundation for a more
coherent, navigable metaverse built on deterministic spatial principles,
where objects and avatars can move between independently developed
environments while maintaining consistent representation.

### AI Integration and Training Efficiency

The structured, deterministic nature of A15-encoded space offers
potential advantages for artificial intelligence systems operating in
virtual environments:

- Investigating whether the discrete, deterministic nature of A15
  representation reduces training noise for spatial AI models

- Exploring potential efficiency gains in reinforcement learning when
  using consistent spatial representations

- Developing AI navigation and pathfinding algorithms optimized for A15
  cell structures

- Creating prediction models that leverage the geometric regularity of
  the A15 lattice

This research could lead to more efficient training methodologies for
spatial AI agents and potentially more robust behavior in deployed
systems.

### Exploration of Alternative Geometric Structures

While A15 offers a compelling balance of properties derived from
crystallography, exploration of alternative partitioning or encoding
schemes derived from related or more exotic mathematical structures
could yield novel insights or properties advantageous for specific
applications. Areas for investigation include:

- **Triply Periodic Minimal Surfaces (TPMS):** Structures like the
  Gyroid (Schoen 1970) possess complex topology useful for flow
  simulation or intricate environment design, likely requiring implicit
  surface representations.

- **Quasicrystalline Patterns:** Non-periodic tilings exhibiting
  symmetries forbidden in periodic crystals (like 5-fold rotation
  (Shechtman et al. 1984)) could enable partitioning schemes with unique
  tiling properties or isotropy characteristics.

- **Optimized Point Sets (e.g., Delone Sets):** Computationally
  generated point distributions balancing criteria like density and
  minimum separation guarantees (Gruber 2007) could tailor partitioning
  more closely to specific application requirements than regular
  lattices.

- **Higher-Dimensional Projections:** Projecting regular polytopes or
  honeycombs (e.g., from 4D (Coxeter 1973)) can generate novel 3D
  structures with useful partitioning properties.

Pursuing these directions promises to expand the capabilities of
structurally informed spatial partitioning. The intersection of
crystallographic principles, computational geometry, numerical analysis,
and real-time interactive systems represents a fertile territory for
extending beyond entertainment into scientific visualization,
collaborative design, distributed simulation, robotics, and the broader
architecture of spatial computing itself.

# Glossary of Terms and Notation

This glossary provides definitions for the specialized terminology and
mathematical notation used throughout this manuscript.

## Mathematical Notation

$`\epsilon_\delta`$  
**Global Output Scaling Factor** — The primary scaling parameter
(specified via `-scale` in `A15.py`) applied to map internal integer
coordinates to output units. Determines the physical size represented by
one unit of the 96-unit baseline dimension. Typically expressed as a
fraction (e.g., $`1/64`$) or power of two (e.g., $`2^{-6}`$).

$`\epsilon_N`$  
**Inherent Base Scale** — The minimum scaling factor, expressible as a
power of two ($`1/2^N`$), required to represent the internal integer
geometry exactly when mapped to a binary number system. Inferred by
analyzing the generated structure’s geometric requirements relative to
the 96-unit baseline dimension.

$`\epsilon_\Delta`$  
**Stability Difference** — Measures the mismatch between
$`\epsilon_\delta`$ and $`\epsilon_N`$. When $`\epsilon_\Delta = 0`$,
the scaling is stable, meaning the output coordinates can be represented
exactly in binary floating-point format relative to the internal grid.
Calculated as the residual after subtracting the nearest integer
multiple of $`\epsilon_N`$ from $`\epsilon_\delta`$.

$`N_1`$  
**Basic Unit Width** — The physical dimension corresponding to the
96-unit effective lattice baseline at the chosen output scale
$`\epsilon_\delta`$. Calculated as $`96 \times \epsilon_\delta`$. For
the recommended scale of $`\epsilon_\delta = 2^{-6}`$, this yields
$`N_1 = 96 \times 2^{-6} = 96/64 = 1.5`$ (typically measured in
millimeters for human-scale applications).

ULP  
**Unit of Least Precision** — The smallest coordinate difference or
spatial distance along a principal axis that can be exactly resolved by
the quantization scheme at a specific scale $`\epsilon_\delta`$. Equal
to the chosen output scale factor, $`\epsilon_\delta`$.

## Technical Terminology

A15  
The crystallographic designation for the $`\beta`$–$`W`$ structure with
space group $`Pm\bar{3}n`$ (No. 223), characterized by high coordination
and near-icosahedral local ordering. Named according to the
Strukturbericht notation system used in crystallography.

$`\beta`$–$`W`$  
Beta-tungsten, the metallurgical name for the A15 crystal structure,
originally observed in tungsten but subsequently found in various
intermetallic compounds with the general formula A$`_3`$B.

C12  
12-coordinated sites within the A15 structure, occupying Wyckoff
position 2a, comprising 25% of the basis sites. The local environment
around these sites resembles a pyritohedron.

C14  
14-coordinated sites within the A15 structure, occupying Wyckoff
position 6d, comprising 75% of the basis sites. The local environment
around these sites resembles a tetradecahedron.

Determinism  
Guarantee that identical operations produce bit-identical results across
different computing environments, regardless of hardware architecture,
operating system, or compiler optimizations.

Isotropy  
Uniformity of properties across different directions, indicating how
"fair" or unbiased a spatial representation is. Higher isotropy means
minimal directional artifacts or biases.

Lattice  
A periodic arrangement of points in space, defined by translational
symmetry. In crystallography, a lattice is characterized by its Bravais
type and basis vectors.

$`Pm\bar{3}n`$  
Hermann-Mauguin notation for the space group (No. 223) of the A15
structure, describing its complete symmetry operations. The notation
indicates: P (primitive unit cell), m (mirror plane), $`\bar{3}`$
(3-fold rotoinversion axis), n (glide plane).

Quantization  
The process of mapping continuous coordinates to discrete A15 lattice
identifiers. This transformation discretizes space according to the
chosen partitioning method (WPH or TSP) and scaling factor.

$`T_h`$  
Point group symmetry ($`m\bar{3}`$, order 24) describing the local
symmetry around sites in the A15 structure. This group is a maximal
subgroup common to both cubic symmetry ($`O_h`$) and icosahedral
symmetry ($`I_h`$), enabling the structure’s high local isotropy within
a periodic lattice.

TSP  
Tetrastix Prism, a simplified partitioning method using axis-aligned
planar faces to divide space into cubic blocks centered on A15 lattice
sites. Offers computational efficiency at some cost to isotropy.

WPH  
Weaire-Phelan Honeycomb, a partitioning method using two distinct
polyhedra (pyritohedra and tetradecahedra) in a 1:3 ratio, corresponding
to the Voronoi decomposition of A15 lattice points. Provides high
isotropy but with more complex cell boundaries.

Wyckoff Position  
Standardized designation for sets of equivalent points within a crystal
structure, identified by multiplicity and site symmetry. Named after
Ralph W.G. Wyckoff, who systematized the classification of crystal
structures.

96-Unit Baseline  
The fundamental period or effective unit dimension of the comprehensive
internal integer grid required to represent the complete A15 structure
without loss of precision. Derived as $`4 \times 24`$ to accommodate
basis site offsets of $`1/4`$ within the lattice spacing factor of
$`24`$.

## Implementation Parameters in `A15.py`

`prescale`  
Internal integer multiplier applied to shape primitives, establishing
primitive resolution (defaults: 20 for WPH, 24 for TSP). Converts base
fractional coordinates to integer vertices relative to shape center.

`rescale`  
Optional power-of-two scaling applied before generation, adjusting size
relative to the 96-unit baseline. Specified via `-rescale` flag or
implicitly through ‘+/-’ suffixes appended to shape names (e.g.,
`pyritohedra++` implies a rescale factor of $`2^2=4`$).

`scale`  
Equivalent to $`\epsilon_\delta`$, the final global scaling factor
applied after generation to map internal integers to output units.
Specified via the `-scale=<value>` command-line option.

`n`  
Controls lattice extent: integer specifies cuboid dimensions (number of
lattice cells along each axis); float specifies spherical radius cutoff
based on lattice vector magnitude relative to the origin.

`-bars`  
Visualization option enabling histogram analysis of floating-point
denominators to validate numerical stability. The resulting histogram
pattern directly indicates the stability regime ($`\epsilon_\Delta = 0`$
or $`\epsilon_\Delta \neq 0`$).

`-stix`  
Configuration option activating the TSP partitioning method instead of
the default WPH geometry. Trades some isotropy for computational
efficiency in point-in-cell tests.

`-pop`  
Option to display the visualization in a pop-up window. Can be specified
as `-pop` or with a specific command like `-pop=open`.

`-savefig`  
Option to save the visualization to a file. Default filename is
`savefig.png` unless specified otherwise.

## Integration Strategies

Dual-Representation  
Integration approach maintaining two parallel coordinate
representations: native engine coordinates (typically floats) for
internal processing and A15 identifiers for network transmission and
verification. Minimizes integration complexity while providing bandwidth
and determinism benefits at interface boundaries.

Middleware  
Integration approach implementing a translation layer between an
existing engine and A15 encoding. Transparently converts coordinates
without requiring significant engine modifications, facilitating
incremental adoption.

Native Implementation  
Integration approach building spatial data structures and physics
directly on A15 identifiers, maximizing determinism throughout the
entire pipeline. Requires more engineering investment but provides the
strongest guarantees.

Relative Addressing  
Technique using hierarchical A15 grids at different scales, with
fine-grained coordinates defined relative to a parent object’s position.
Enables efficient multi-scale representation while maintaining precision
where needed.

This glossary provides a reference for the specialized terms and
mathematical notation used throughout the manuscript. For implementation
details and additional technical specifications, refer to the `A15.py`
reference implementation
(<a href="#subsec-implementation-a15py" data-reference-type="ref+Label"
data-reference="subsec-implementation-a15py">2.6</a>) and the
supplementary information
(<a href="#sec-supplementary" data-reference-type="ref+Label"
data-reference="sec-supplementary">6</a>).

# Supplementary Information

This section provides guidelines for replicating the results presented
using the accompanying code, details supplementary resources available
online, and discusses additional technical considerations relevant to
the implementation and interpretation of the *A15* encoding framework.

## Code Availability and Replication Protocols

The figures and structural data presented in this research were
generated using the accompanying Python script, `A15.py` (Risinger
2024a), which serves as the reference implementation. To ensure
reproducibility, users should have Python 3 installed along with the
standard scientific libraries NumPy (Harris et al. 2020), SciPy
(Virtanen et al. 2020), and Matplotlib (Hunter 2007). Understanding the
script’s execution flow
(<a href="#tab-a15py-workflow" data-reference-type="ref+Label"
data-reference="tab-a15py-workflow">5</a>) and key parameters,
particularly those governing the multi-stage scaling framework
(<a href="#subsec-scaling-pipeline" data-reference-type="ref+Label"
data-reference="subsec-scaling-pipeline">2.3</a>) and numerical
stability validation
(<a href="#subsec-stability-validation" data-reference-type="ref+Label"
data-reference="subsec-stability-validation">2.5</a>), is essential for
proper use and interpretation of results.

**Note on performance:** `A15.py` is a clarity-first reference
implementation intended for correctness and pedagogical value. It does
not reflect the performance achievable through optimized
implementations. Developers targeting real-time applications (e.g.,
networked VR, physics engines) are advised to implement accelerated
versions using fixed-point or SIMD-based quantization pipelines
alongside other integration strategies outlined in
<a href="#subsubsec-apps-integration" data-reference-type="ref+Label"
data-reference="subsubsec-apps-integration">4.1.2</a> and design
optimized implementations tailored to their application domains.

The core Python script (`A15.py`), configuration files (`*.png.txt`)
used for figure generation, the LaTeX source for this manuscript (or a
version thereof), and extended documentation are publicly available
within the Infima Labs `space` repository on GitHub (Infima Labs 2023):

<div class="center">

<https://github.com/infimalabs/space/>

</div>

A project overview and supplementary materials may also be found at the
project’s homepage:

<div class="center">

<https://infima.space/A15/>

</div>

The related project `layoutc`, applying similar compact encoding
concepts to paintball field layouts (Risinger 2024b), is also available
via associated repositories:

<div class="center">

<https://github.com/infimalabs/layoutc/>

</div>

Achieving deterministic results, a core goal of this framework, relies
crucially on operating within the **Binary** or **Stable** scaling
regimes ($`\epsilon_\Delta = 0`$,
<a href="#subsubsec-stability-regimes" data-reference-type="ref+Label"
data-reference="subsubsec-stability-regimes">2.4.4</a>). The recommended
baseline scale of $`\epsilon_\delta = 2^{-6}`$ (`-scale=1/64`) generally
provides a practical balance for human-scale interactions, offering high
precision (approximately 1.5 mm basic unit width $`N_1`$, see
<a href="#subsubsec-apps-measurement" data-reference-type="ref+Label"
data-reference="subsubsec-apps-measurement">4.1.3</a> and
<a href="#subsubsec-notes-figures" data-reference-type="ref+Label"
data-reference="subsubsec-notes-figures">6.4.7</a>) while ensuring exact
floating-point representability relative to the internal grid for
typical configurations. Users are strongly encouraged to employ the
`-bars` analysis feature
(<a href="#subsec-stability-validation" data-reference-type="ref+Label"
data-reference="subsec-stability-validation">2.5</a>) to explicitly
verify the stability ($`\epsilon_\Delta = 0`$) of any custom
configurations before deployment in applications where determinism is
critical.

Furthermore, while `A15.py` deterministically implements alternating
handedness for local coordination environments based on lattice position
(<a href="#subsubsec-limits-handedness" data-reference-type="ref+Label"
data-reference="subsubsec-limits-handedness">3.3.6</a>), networked
applications or systems exchanging *A15*-encoded data **must** establish
and consistently apply a shared global orientation convention to ensure
interoperability and prevent geometric mirroring between different
clients or system components.

### Example Replication Commands

The primary figures presented in this manuscript can be regenerated
using the `A15.py` script and the corresponding configuration files
(typically named `fig-`*`name`*`.png.txt`) provided in the supplementary
materials repository
(<a href="#subsec-replication" data-reference-type="ref+Label"
data-reference="subsec-replication">6.1</a>). Ensure the script and
configuration files are accessible in the execution environment. Use the
`-i` (or `-pop`) option for interactive viewing (requires a graphical
display environment):

<div class="description">

Nested *A15* (Multiple Scales)  
`python3 A15.py -i fig-intro.png.txt`

Left-Handed $`^1/_2`$ Unit Cell  
`python3 A15.py -i fig-cell2.png.txt`

Weaire–Phelan Honeycomb and Tetrastix Prism  
`python3 A15.py -i fig-wp.png.txt`  
`python3 A15.py -i fig-ts.png.txt`

Representative Example (Unstable Scale)  
`python3 A15.py -i fig-main.png.txt`

Binary, Stable, and Unstable Scales  
`python3 A15.py -i fig-histb.png.txt`  
`python3 A15.py -i fig-hists.png.txt`  
`python3 A15.py -i fig-histu.png.txt`

</div>

For a detailed explanation of all command-line options, parameters,
configuration file syntax, and advanced usage, refer to
`python3 A15.py –help`.

## Middleware Integration Example

To illustrate a simple partial-integration strategy, consider a
middleware function that translates floating-point positions into A15
lattice identifiers for network transmission, and vice versa on receipt:

``` python
from A15 import encode, decode

# Application float-space position
world_pos = [1.25, 0.75, 2.5]  # meters

# Encode to compact A15 ID for transmission
a15_id = encode(world_pos, scale=1/64)
send_to_network(a15_id)

# On receiver: decode back to float
recv_pos = decode(a15_id, scale=1/64)
render_object_at(recv_pos)
```

This pattern allows developers to gain determinism and bandwidth
benefits without deeply modifying the internal coordinate systems of a
host engine. The encoded form remains compact and platform-agnostic
during transmission, with reconstruction yielding consistent results.

## Practical Implementation Considerations

Beyond the theoretical foundation and validation tools provided in this
research, several practical considerations merit attention for
developers seeking to implement the A15 framework in real-world
applications:

### Performance Optimization Strategies

While comprehensive benchmarking remains an area for future work
(<a href="#subsubsec-outlook-performance" data-reference-type="ref+Label"
data-reference="subsubsec-outlook-performance">4.2.1</a>), preliminary
experience suggests several approaches to optimize A15 implementation
performance:

- **Aligned Coordinate Systems:** Designing application coordinate
  systems to align with the A15 grid allows quantization to be
  implemented as simple integer truncation rather than complex geometric
  tests.

- **Caching A15 Identifiers:** For static objects or slowly changing
  environments, pre-computing and caching A15 identifiers can amortize
  quantization costs.

- **Hierarchical Spatial Partitioning:** Combining the A15 grid with
  higher-level partitioning schemes (e.g., spatial hashing) can
  accelerate nearest-neighbor searches for quantization.

- **Simpler Cell Geometry for Performance-Critical Paths:** Using the
  TSP partitioning method
  (<a href="#subsec-intro-partitioning" data-reference-type="ref+Label"
  data-reference="subsec-intro-partitioning">1.3</a>) for paths
  requiring frequent quantization trading some isotropy for
  computational efficiency.

### Incremental Adoption Path

Rather than attempting a complete transition to A15-based spatial
representation, many applications will benefit from an incremental
adoption strategy:

1.  **Network Transmission Only** - Use A15 encoding solely for
    transmitting position updates over the network, while maintaining
    traditional floating-point coordinates for all internal processing.

2.  **Authoritative State Recording** - Extend A15 usage to
    authoritative state storage, enabling deterministic replay and
    verification capabilities.

3.  **Critical System Integration** - Selectively implement A15-aware
    subsystems for components where determinism is most crucial (e.g.,
    collision detection).

4.  **Comprehensive Integration** - Develop fully A15-native physics and
    simulation systems if warranted by application requirements.

This phased approach allows developers to gain immediate benefits from
A15 encoding (bandwidth reduction, improved network consistency) while
deferring more complex integration challenges until warranted by
specific application needs.

### Testing and Verification

Implementing a system that guarantees determinism requires rigorous
testing approaches:

- **Cross-Platform Verification:** Test A15-encoded operations across
  different hardware, operating systems, and compiler settings to verify
  bit-exact results.

- **Stability Regime Validation:** For any chosen scale, verify that
  $`\epsilon_\Delta = 0`$ using techniques similar to the histogram
  analysis in `A15.py`
  (<a href="#subsec-stability-validation" data-reference-type="ref+Label"
  data-reference="subsec-stability-validation">2.5</a>).

- **Boundary Case Testing:** Thoroughly test edge cases near cell
  boundaries where quantization decisions might be sensitive to
  implementation details.

- **Deterministic Replay:** Validate that identical initial conditions
  and inputs reliably produce identical event sequences when using A15
  encoding.

These testing approaches help ensure that the theoretical determinism
guarantees of the A15 framework translate into practical consistency in
deployed applications.

## Supporting Notes and Clarifications

Further details, data, and clarifications related to this research are
provided below.

### Bandwidth Calculation Basis

The discussion regarding network bandwidth requirements
(<a href="#subsubsec-benefits-efficiency" data-reference-type="ref+Label"
data-reference="subsubsec-benefits-efficiency">3.2.2</a>) utilizes a
baseline calculation
(<a href="#eq-bandwidth-baseline" data-reference-type="ref+Label"
data-reference="eq-bandwidth-baseline">[eq-bandwidth-baseline]</a>)
assuming a single avatar with 20 tracked joints, each transmitting 3D
position (3x 32 floats) and an orientation quaternion (4x 32 floats) at
a 15 Hz update rate (67.2   s<sup>−1</sup>). This kinematic component
typically represents only a fraction of the total network traffic in
complex interactive applications.

### Memory Efficiency Context

The estimate of **50% or more** memory and bandwidth savings
(<a href="#eq-efficiency-memory" data-reference-type="ref+Label"
data-reference="eq-efficiency-memory">[eq-efficiency-memory]</a>,
<a href="#subsubsec-benefits-efficiency" data-reference-type="ref+Label"
data-reference="subsubsec-benefits-efficiency">3.2.2</a>) compares
storing 3D coordinates using compact integer *A15* identifiers versus
raw 32 floating-point vectors. A representative scenario assumes a 48
integer representation per 3D point for the *A15* identifier (balancing
range and precision) compared to
$`3 \times 32\,\text{bit} = 96\,\text{bit}`$ for three standard floats.
This efficiency gain is most pronounced when quantizing bounded volumes
where the extreme dynamic range and full mantissa precision of floats
are not required. The exact saving achieved depends on the application’s
required spatial extent, desired intra-cell resolution (potentially
managed via relative addressing,
<a href="#subsubsec-apps-relative" data-reference-type="ref+Label"
data-reference="subsubsec-apps-relative">4.1.4</a>), and the chosen bit
depth for the *A15* identifier.

### Geometric Data Availability

Tables containing the precise internal integer vertex coordinates
(relative to shape centers at the relevant `prescale` value) for the
fundamental polyhedra (pyritohedra with various $`h`$ parameters,
tetradecahedra) generated by `A15.py` functions are available within the
code repository
(<a href="#subsec-replication" data-reference-type="ref+Label"
data-reference="subsec-replication">6.1</a>), allowing independent
verification of geometric constructions.

### Pyritohedra Parameter for Weaire–Phelan Honeycomb Geometry

As implemented in `A15.py`, invoking the `pyritohedron()` function with
the specific height parameter $`h=7/5`$ yields internal integer
coordinates that, after appropriate scaling and placement by
`lattice()`, correspond precisely to the vertex coordinates defining the
pyritohedral cells within the geometric Weaire–Phelan Honeycomb
partition used in this framework
(<a href="#subsec-intro-partitioning" data-reference-type="ref+Label"
data-reference="subsec-intro-partitioning">1.3</a>).

### Unit of Least Precision (ULP) Definition

Within this framework, when operating at a specific Binary or Stable
output scale $`\epsilon_\delta`$
(<a href="#subsubsec-stability-regimes" data-reference-type="ref+Label"
data-reference="subsubsec-stability-regimes">2.4.4</a>), the **Unit of
Least Precision (ULP)** represents the smallest coordinate difference or
spatial distance along a principal axis that can be exactly resolved by
the quantization scheme. This ULP corresponds directly to an integer
difference of 1 in the underlying *internal integer* coordinate system
(relative to the 96-unit baseline,
<a href="#subsubsec-scaling-baseline" data-reference-type="ref+Label"
data-reference="subsubsec-scaling-baseline">2.3.3</a>) before the final
scaling by $`\epsilon_\delta`$ is applied. Therefore, the physical size
of the ULP is precisely equal to the chosen output scale factor,
$`\epsilon_\delta`$. Features, movements, or discrepancies smaller than
$`\epsilon_\delta`$ cannot be distinctly represented by the encoding at
that scale. Selecting an appropriate $`\epsilon_\delta`$ involves
balancing the desired spatial resolution (ULP) against the overall
spatial range achievable within a fixed-bit integer representation
chosen for the *A15* identifier, potentially leveraging relative
addressing
(<a href="#subsubsec-apps-relative" data-reference-type="ref+Label"
data-reference="subsubsec-apps-relative">4.1.4</a>) to manage this
trade-off across different parts of a scene.

### Framework Adaptability

While this research focuses intensely on the *A15* phase structure due
to its compelling combination of advantageous properties for
deterministic spatial encoding, the underlying `A15.py` software
framework possesses inherent adaptability. Key components, particularly
the `lattice()` function for replicating geometric units according to
symmetry rules and the visualization tools within the `figure()`
function, could be modified or extended to generate and visualize other
crystal lattice types or different space-filling structures, offering a
versatile platform for broader geometric exploration, albeit likely
requiring non-trivial adaptation of the core geometric and scaling
logic.

### Interpretation of Figure Annotations

Annotations visible in figures generated by `A15.py` with the `-bars`
option (e.g., <a href="#fig-main" data-reference-type="ref+Label"
data-reference="fig-main">10</a>) directly illustrate key stability
concepts discussed in
<a href="#subsec-stability" data-reference-type="ref+Label"
data-reference="subsec-stability">2.4</a>. The calculated value $`N_1`$
(labeled “basic unit width” or similar in some outputs) shown in the
histogram sidebar represents the physical dimension corresponding to the
chosen output scale $`\epsilon_\delta`$ applied to the fundamental
internal lattice dimension (the 96-unit baseline,
<a href="#subsubsec-scaling-baseline" data-reference-type="ref+Label"
data-reference="subsubsec-scaling-baseline">2.3.3</a>). For the unstable
scale $`\epsilon_\delta=1/96`$ shown in
<a href="#fig-main" data-reference-type="ref+Label"
data-reference="fig-main">10</a>, this results in
$`N_1 = 96 \times (1/96) = 1.0`$ (assuming millimeters as the base unit,
thus 1.0 mm). In contrast, the recommended stable scale
$`\epsilon_\delta = 1/64`$ yields $`N_1 = 96 \times (1/64) = 1.5`$ (or
1.5 mm,
<a href="#subsubsec-apps-measurement" data-reference-type="ref+Label"
data-reference="subsubsec-apps-measurement">4.1.3</a>). The displayed
$`\epsilon_\Delta`$ value explicitly confirms the calculated stability
difference for the configuration (e.g.,
$`\epsilon_\Delta \approx \num{0.0104} \neq 0`$ for the unstable case in
<a href="#fig-main" data-reference-type="ref+Label"
data-reference="fig-main">10</a>, confirming $`\epsilon_\delta`$ is not
an integer multiple of $`\epsilon_N`$), providing direct numerical
validation alongside the visual histogram representation.

<div class="acknowledgements">

This work benefited significantly from the democratization of knowledge,
particularly Wikipedia’s extensive collection of mathematical and
physical concepts. The accessibility of such resources proved
invaluable.

Unabashed credit is given to the transformative role of artificial
intelligence assistants in accelerating the development and articulation
of these ideas. Their capabilities enabled rapid iteration and
refinement of concepts and language.

Special appreciation is extended to the Infima Labs team for their
confidence and trust throughout this research. Their shared vision for
advancing spatial computing was—and continues to be—a wellspring of
inspiration.

Above all, deep gratitude is expressed to the author’s wife and family
for their unwavering support and patience throughout this work’s
development. Their encouragement and understanding were instrumental in
bringing these ideas to fruition.

</div>

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-ITCVolumeA2016" class="csl-entry">

Aroyo, Mois Ilia, ed. 2016. *<span class="nocase">International Tables
for Crystallography, Volume A: Space-Group Symmetry</span>*. 6th ed.
Chester, UK: International Union of Crystallography, Wiley.

</div>

<div id="ref-AshcroftMermin1976" class="csl-entry">

Ashcroft, Neil William, and N. David Mermin. 1976. *Solid State
Physics*. New York: Holt, Rinehart; Winston.

</div>

<div id="ref-Bentley1975" class="csl-entry">

Bentley, Jon Louis. 1975. “Multidimensional Binary Search Trees Used for
Associative Searching.” *Communications of the ACM* 18 (9): 509–17.
<https://doi.org/10.1145/361002.361007>.

</div>

<div id="ref-CCPA-2018" class="csl-entry">

California State Legislature. 2018. “<span class="nocase">California
Consumer Privacy Act of 2018 (CCPA)</span>.”
<https://leginfo.legislature.ca.gov/faces/codes_displayText.xhtml?division=3.&part=4.&lawCode=CIV&title=1.81.5>.

</div>

<div id="ref-Claypool2006" class="csl-entry">

Claypool, Mark, and Kajal Claypool. 2006. “Latency and Player Actions in
Online Games.” *Communications of the ACM* 49 (11): 40–45.
<https://doi.org/10.1145/1167838.1167860>.

</div>

<div id="ref-ConwaySloane1999" class="csl-entry">

Conway, John Horton, and Neil James Alexander Sloane. 1999. *Sphere
Packings, Lattices and Groups*. 3rd ed. Vol. 290. Grundlehren Der
Mathematischen Wissenschaften. New York: Springer-Verlag.
<https://doi.org/10.1007/978-1-4757-6568-7>.

</div>

<div id="ref-Coxeter1973" class="csl-entry">

Coxeter, Harold Scott MacDonald. 1973. *Regular Polytopes*. 3rd ed. New
York: Dover Publications.

</div>

<div id="ref-CoxeterMoser1972" class="csl-entry">

Coxeter, Harold Scott MacDonald, and William Oscar Jules Moser. 1972.
*Generators and Relations for Discrete Groups*. 3rd ed. Vol. 14.
Ergebnisse Der Mathematik Und Ihrer Grenzgebiete. Berlin:
Springer-Verlag. <https://doi.org/10.1007/978-3-662-21943-0>.

</div>

<div id="ref-UnrealCoords" class="csl-entry">

Epic Games. 2023. “Coordinate Space Terminology.”
<https://docs.unrealengine.com/5.3/en-US/coordinate-space-terminology-in-unreal-engine/>.

</div>

<div id="ref-GDPR-2016" class="csl-entry">

European Parliament and Council of the European Union. 2016.
“<span class="nocase">Regulation (EU) 2016/679 of the European
Parliament and of the Council of 27 April 2016 on the protection of
natural persons with regard to the processing of personal data and on
the free movement of such data (General Data Protection
Regulation)</span>.” <https://eur-lex.europa.eu/eli/reg/2016/679/oj>.

</div>

<div id="ref-Finkel1974" class="csl-entry">

Finkel, Raphael A., and Jon Louis Bentley. 1974. “Quad Trees: A Data
Structure for Retrieval on Composite Keys.” *Acta Informatica* 4 (1):
1–9. <https://doi.org/10.1007/BF00288933>.

</div>

<div id="ref-FrankKasper1958" class="csl-entry">

Frank, Frederick Charles, and John Samuel Kasper. 1958. “Complex Alloy
Structures Regarded as Sphere Packings. I. Definitions and Basic
Principles.” *Acta Crystallographica* 11 (3): 184–90.
<https://doi.org/10.1107/s0365110x5800048x>.

</div>

<div id="ref-FrankKasper1959" class="csl-entry">

———. 1959. “Complex Alloy Structures Regarded as Sphere Packings. II.
Analysis and Classification of Representative Structures.” *Acta
Crystallographica* 12 (7): 483–99.
<https://doi.org/10.1107/S0365110X5900149X>.

</div>

<div id="ref-Goldberg1991" class="csl-entry">

Goldberg, David. 1991. “What Every Computer Scientist Should Know about
Floating-Point Arithmetic.” *ACM Computing Surveys* 23 (1): 5–48.
<https://doi.org/10.1145/103162.103163>.

</div>

<div id="ref-Gruber2007" class="csl-entry">

Gruber, Peter Michael. 2007. *Convex and Discrete Geometry*. Vol. 336.
Grundlehren Der Mathematischen Wissenschaften. Berlin: Springer.
<https://doi.org/10.1007/978-3-540-71133-9>.

</div>

<div id="ref-Guttman1984" class="csl-entry">

Guttman, Antonin. 1984. “R-Trees: A Dynamic Index Structure for Spatial
Searching.” In *Proceedings of the 1984 ACM SIGMOD International
Conference on Management of Data*, 47–57. SIGMOD ’84. Boston,
Massachusetts, USA: ACM. <https://doi.org/10.1145/602259.602266>.

</div>

<div id="ref-Harris2020" class="csl-entry">

Harris, Charles R., K. Jarrod Millman, Stéfan J. van der Walt, Ralf
Gommers, Pauli Virtanen, David Cournapeau, Eric Wieser, et al. 2020.
“Array Programming with NumPy.” *Nature* 585 (7825): 357–62.
<https://doi.org/10.1038/s41586-020-2649-2>.

</div>

<div id="ref-Hunter2007" class="csl-entry">

Hunter, John D. 2007. “Matplotlib: A 2D Graphics Environment.”
*Computing in Science & Engineering* 9 (3): 90–95.
<https://doi.org/10.1109/MCSE.2007.55>.

</div>

<div id="ref-IEEE754-2019" class="csl-entry">

IEEE. 2019. “<span class="nocase">IEEE Standard for Floating-Point
Arithmetic</span>.” 754-2019. New York, NY, USA: Institute of Electrical
and Electronics Engineers; IEEE.
<https://doi.org/10.1109/IEEESTD.2019.8766229>.

</div>

<div id="ref-InfimaSpace" class="csl-entry">

Infima Labs. 2023. “Infima Space Repository.” <https://infima.space>.

</div>

<div id="ref-Kusner1996" class="csl-entry">

Kusner, Robert, and John Morgan Sullivan. 1996. “Comparing the
Weaire-Phelan Equal-Volume Foam to Kelvin’s Foam.” *Forma* 11 (3):
233–42.

</div>

<div id="ref-Luebke2002" class="csl-entry">

Luebke, David, Martin Reddy, Jonathan David Cohen, Amitabh Varshney,
Benjamin Watson, and Robert Huebner. 2002. *Level of Detail for 3D
Graphics*. San Francisco, CA: Morgan Kaufmann.

</div>

<div id="ref-PythonDocsFloatRatio" class="csl-entry">

Python Software Foundation. 2025. “<span class="nocase">Built-in Types —
Python 3 documentation: float.as_integer_ratio</span>.”
<https://docs.python.org/3/library/stdtypes.html#float.as_integer_ratio>.

</div>

<div id="ref-Risinger2024A15" class="csl-entry">

Risinger, C. Anthony. 2024a. “<span class="nocase">A15.py: A15 Phase
Visualizer</span>.”
<https://github.com/infimalabs/space/blob/main/A15/A15.py>.

</div>

<div id="ref-Risinger2024Layoutc" class="csl-entry">

———. 2024b. “<span class="nocase">layoutc: Paintball Field Layout
Compressor</span>.” <https://github.com/infimalabs/layoutc>.

</div>

<div id="ref-Samet1990" class="csl-entry">

Samet, Hanan. 1990. *The Design and Analysis of Spatial Data
Structures*. Reading, MA: Addison-Wesley.

</div>

<div id="ref-Schoen1970" class="csl-entry">

Schoen, Alan Howard. 1970. “Infinite Periodic Minimal Surfaces Without
Self-Intersections.” NASA TN D-5541. NASA.
<https://ntrs.nasa.gov/citations/19700020214>.

</div>

<div id="ref-Shechtman1984" class="csl-entry">

Shechtman, Dan, Ilan Blech, Denis Gratias, and John Werner Cahn. 1984.
“Metallic Phase with Long-Range Orientational Order and No Translational
Symmetry.” *Physical Review Letters* 53 (20): 1951–53.
<https://doi.org/10.1103/PhysRevLett.53.1951>.

</div>

<div id="ref-Kelvin1887" class="csl-entry">

Thomson, William. 1887. “On the Division of Space with Minimum
Partitional Area.” *Philosophical Magazine* 24 (151): 503–14.
<https://doi.org/10.1080/14786448708628135>.

</div>

<div id="ref-UnityCoords" class="csl-entry">

Unity Technologies. 2024. “<span class="nocase">Understanding Vector
Arithmetic in Unity</span>.”
<https://docs.unity3d.com/Manual/UnderstandingVectorArithmetic.html>.

</div>

<div id="ref-Virtanen2020" class="csl-entry">

Virtanen, Pauli, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler
Reddy, David Cournapeau, Evgeni Burovski, et al. 2020. “SciPy 1.0:
Fundamental Algorithms for Scientific Computing in Python.” *Nature
Methods* 17: 261–72. <https://doi.org/10.1038/s41592-019-0686-2>.

</div>

<div id="ref-WeaireHutzler2001" class="csl-entry">

Weaire, Denis, and Stefan Hutzler. 2001. *The Physics of Foams*. Oxford:
Oxford University Press.

</div>

<div id="ref-WeairePhelan1994" class="csl-entry">

Weaire, Denis, and Robert Phelan. 1994. “A Counter-Example to Kelvin’s
Conjecture on Minimal Surfaces.” *Philosophical Magazine Letters* 69
(2): 107–10. <https://doi.org/10.1080/09500839408241577>.

</div>

</div>

[^1]: Floats are abundant in software yet maddeningly fickle; among the
    first one-hundred simple reciprocals ($`1/n`$), ninety-three require
    approximation in standard binary float formats (Goldberg 1991).

[^2]: This simplified view omits details such as exponent bias,
    subnormal numbers, infinities, and NaNs, rigorously defined in (IEEE
    2019) but not central to the core issue of approximation.
