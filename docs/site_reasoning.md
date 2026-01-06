# Interstitial Site Analysis & Reasoning

Based on Voronoi tessellation and local environment analysis.

## Methodology
1.  **Generation**: Voronoi tessellation was used to identify the largest voids in the lattice compatible with Cr insertion.
2.  **Analysis**: For each distinct site, `CrystalNN` was used to determine the coordination environment and nearest neighbors.
3.  **Classification**: Sites are classified by their coordination number and the chemical identity of the neighbors.

## Identified Sites
### Site 0: Distorted-7-coord
- **Label**: `Distorted-7-coord_CN7-Te4-Sb3`
- **Coordinates (Frac)**: (0.667, 0.333, 0.418)
- **Coordination**: 7
- **Environment**: CN7-Te4-Sb3 (Neighbors: ['Sb', 'Sb', 'Sb', 'Te', 'Te', 'Te', 'Te'])
- **Nearest Neighbor Distance**: 2.583 Å
- **Reason for Selection**: Identified as a distinct Voronoi node representing a potential local minimum for intercalation. Distorted-7-coord geometry suggests stability for transition metals.

### Site 1: Tetrahedral-like
- **Label**: `Tetrahedral-like_CN4-Te4`
- **Coordinates (Frac)**: (0.667, 0.333, 0.470)
- **Coordination**: 4
- **Environment**: CN4-Te4 (Neighbors: ['Te', 'Te', 'Te', 'Te'])
- **Nearest Neighbor Distance**: 2.570 Å
- **Reason for Selection**: Identified as a distinct Voronoi node representing a potential local minimum for intercalation. Tetrahedral-like geometry suggests stability for transition metals.

### Site 2: Distorted-7-coord
- **Label**: `Distorted-7-coord_CN7-Te6-Sb1`
- **Coordinates (Frac)**: (0.000, 0.000, 0.520)
- **Coordination**: 7
- **Environment**: CN7-Te6-Sb1 (Neighbors: ['Te', 'Te', 'Te', 'Sb', 'Te', 'Te', 'Te'])
- **Nearest Neighbor Distance**: 2.668 Å
- **Reason for Selection**: Identified as a distinct Voronoi node representing a potential local minimum for intercalation. Distorted-7-coord geometry suggests stability for transition metals.

### Site 3: Distorted-8-coord
- **Label**: `Distorted-8-coord_CN8-Te4-Sb4`
- **Coordinates (Frac)**: (0.333, 0.667, 0.973)
- **Coordination**: 8
- **Environment**: CN8-Te4-Sb4 (Neighbors: ['Te', 'Te', 'Te', 'Sb', 'Sb', 'Sb', 'Te', 'Sb'])
- **Nearest Neighbor Distance**: 2.628 Å
- **Reason for Selection**: Identified as a distinct Voronoi node representing a potential local minimum for intercalation. Distorted-8-coord geometry suggests stability for transition metals.

