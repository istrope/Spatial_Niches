# Neighborhood Differential Expression Utility

This repository provides a `Neighborhood` class for performing differential expression analysis in spatial transcriptomics datasets (e.g., Xenium, Visium) based on gene expression and neighboring cell-type filters.

## Features

* Create binary masks of cells expressing a target gene and belonging to specified sender cell types.
* Filter cells by the proportion or absolute number of neighboring cells of a given receiver cell type.
* Optionally filter by expression of a neighboring gene.
* Support for both "sender" mode (comparing gene‑expressing vs non‑expressing cells within a cell type) and "receiver" mode (comparing neighbors of expressing cells vs non‑neighbors).
* Bootstrap or direct Wilcoxon differential expression testing via `scanpy.tl.rank_genes_groups`.

## Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/<your-username>/neighborhood-de.git
   cd neighborhood-de
   ```

2. Create a virtual environment and install dependencies:

   ```bash
   python -m venv venv
   source venv/bin/activate
   pip install --upgrade pip
   pip install -r requirements.txt
   ```
> **Required**: Python 3.8+
> **Note:** Your AnnData object should contain a connectivity graph in `adata.obsp['connectivities']`.

## Usage

```python
import scanpy as sc
from neighborhood import Neighborhood  # path/to/your_module

# Load or prepare your AnnData
adata = sc.read_h5ad('my_spatial_data.h5ad')

# Initialize Neighborhood utility
nb = Neighborhood(
    gene='CTLA4',
    neighbor_gene='PDCD1',           # optional
    celltype_col='celltype_sender',  # obs column for sender filter
    celltype_vals=['T cell'],
    neighbor_col='celltype_neighbor',# obs column for neighbor filter
    neighbor_vals='Macrophage',      # single or list
    neighbor_perc=0.2,               # fraction of neighbors
    adata=adata,
    sender=True,                     # sender mode
    receiver=False                   # receiver mode
)

# Get indices of group1 (gene+) and group2 (gene-)
group1_idx, group2_idx = nb.get_indices()

# Create a filtered AnnData with a new 'group' obs column
adata_groups = nb.create_adata()

# Run differential expression (no bootstrapping)
b_deg = nb.de_test(nperm=0)
print(b_deg.head())
```

## API Reference

### `Neighborhood` Constructor

```python
Neighborhood(
    gene: str,
    neighbor_gene: Optional[str] = None,
    celltype_col: Optional[str] = None,
    celltype_vals: Union[str, List[str]] = None,
    neighbor_col: Optional[str] = None,
    neighbor_vals: Union[str, List[str]] = None,
    num_neighbors: int = 0,
    adata: anndata.AnnData = None,
    neighbor_perc: float = 0.0,
    sender: bool = True,
    receiver: bool = False
)
```

* **gene**: Target gene for expression mask.
* **neighbor\_gene**: Optional neighboring gene for expression filtering.
* **celltype\_col** / **celltype\_vals**: Obs column and values for sender cell-type filtering.
* **neighbor\_col** / **neighbor\_vals**: Obs column and values for neighbor-cell-type filtering.
* **num\_neighbors**: Minimum number of neighbors required.
* **neighbor\_perc**: Minimum fraction of neighbors required.
* **sender**: If `True`, runs sender-mode filtering; otherwise receiver-mode.
* **receiver**: If `True`, runs receiver-mode filtering.

### Key Methods

* `get_indices() -> Tuple[np.ndarray, np.ndarray]`: Returns `(group1_indices, group2_indices)`.
* `create_adata() -> anndata.AnnData`: Returns a new AnnData with a `group` obs column indicating each cell’s group.
* `de_test(nperm: int, sort: bool = True, test: str = 'wilcoxon') -> pandas.DataFrame`: Runs DE testing (bootstrap if `nperm>0`).

## License

This project is licensed under the [MIT License](LICENSE).
