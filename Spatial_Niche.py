import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
from typing import Union, List
from random import sample

class Neighborhood:
    def __init__(
        self,
        gene: str,
        neighbor_gene: str = None,
        # NEW: one obs‐column & values for sender (celltype) filter
        celltype_col: str = None,
        celltype_vals: Union[str, List[str]] = None,
        # NEW: one obs‐column & values for receiver (neighbor) filter
        neighbor_col: str = None,
        neighbor_vals: Union[str, List[str]] = None,
        num_neighbors: int = 0,
        adata: ad.AnnData = None,
        neighbor_perc: float = 0.0,
        sender: bool = True,
        receiver: bool = False
    ):
        self.gene = gene
        self.neighbor_gene = neighbor_gene
        self.celltype_col = celltype_col
        self.celltype_vals = (
            [celltype_vals] if isinstance(celltype_vals, str) else celltype_vals
        )
        self.neighbor_col = neighbor_col
        self.neighbor_vals = (
            [neighbor_vals] if isinstance(neighbor_vals, str) else neighbor_vals
        )
        self.num_neighbors = num_neighbors
        self.adata = adata
        self.neighbor_perc = neighbor_perc
        self.sender = sender
        self.receiver = receiver

        # pick whichever obsp you like for your neighborhood graph
        self.neighbor_matrix = adata.obsp['connectivities']

    def _make_binary_mask(
        self,
        vals: List[str],
        col: str,
        adata: ad.AnnData
    ) -> np.ndarray:
        """Return an int‐mask of obs[col] ∈ vals; if vals is None, return all 1s."""
        if vals is None or col is None:
            return np.ones(adata.n_obs, dtype=int)
        if col not in adata.obs.columns:
            raise ValueError(f"Column '{col}' not found in adata.obs")
        return adata.obs[col].isin(vals).astype(int).values

    def create_masks(self, gene_mask: np.ndarray, cell_mask: np.ndarray):
        """
        Combine gene_mask (cells expressing gene) with a cell‐type mask.
        Returns two binary vectors of length n_obs.
        """
        # gene_mask is already a boolean array of length n_obs
        # cell_mask is int array of length n_obs
        return gene_mask.astype(bool), (cell_mask.astype(bool))

    def neighbor_filter(self, cell_mask):
        adj = (self.neighbor_matrix > 0).astype(int)
        neighbor_counts = np.array(adj.sum(axis=1)).flatten()
        target_counts = np.array(self.neighbor_matrix.dot(cell_mask)).flatten()

        if self.neighbor_perc > 0:
            frac = np.zeros_like(neighbor_counts, dtype=float)
            valid = neighbor_counts > 0
            frac[valid] = target_counts[valid] / neighbor_counts[valid]
            return frac >= self.neighbor_perc

        if self.num_neighbors > 0:
            return target_counts >= self.num_neighbors

        raise ValueError(
            "Set either neighbor_perc>0 or num_neighbors>0 for neighbor filtering"
        )

    def get_indices(self):
        adata = self.adata

        # 1) gene expression mask
        if self.gene not in adata.var_names:
            raise ValueError(f"Gene '{self.gene}' not in adata.var_names")
        gene_mask = (adata[:, self.gene].X.toarray() > 0).flatten()

        # 2) sender cell‐type mask
        sender_cell_mask = self._make_binary_mask(
            self.celltype_vals, self.celltype_col, adata
        )

        # combine for group1 (gene+ & in celltype) and group2 (gene− & in celltype)
        if self.sender:
            mask_expr   = gene_mask & sender_cell_mask.astype(bool)
            mask_noexpr = (~gene_mask) & sender_cell_mask.astype(bool)

            # optional neighbor‐cell‐type filter
            if self.neighbor_perc > 0 or self.num_neighbors > 0:
                if self.neighbor_col is None or self.neighbor_vals is None:
                    raise ValueError("Must set neighbor_col & neighbor_vals for neighbor filter")
                recv_mask = self._make_binary_mask(
                    self.neighbor_vals, self.neighbor_col, adata
                )
                pass_filter = self.neighbor_filter(recv_mask)
                mask_expr   &= pass_filter
                mask_noexpr &= pass_filter

            # optional neighbor‐gene expression filter
            if self.neighbor_gene:
                neigh_expr = (adata[:, self.neighbor_gene].X.toarray() > 0).flatten()
                neigh_counts = self.neighbor_matrix.dot(neigh_expr)
                neigh_pass = neigh_counts.A1 if hasattr(neigh_counts, 'A1') else neigh_counts
                neigh_pass = neigh_pass > 0
                mask_expr   &= neigh_pass
                mask_noexpr &= neigh_pass

        else:
            # receiver mode: look at neighborhood of gene‐expressing cells
            neigh_counts = self.neighbor_matrix.dot(gene_mask)
            mask_expr   = neigh_counts > 0
            mask_noexpr = neigh_counts == 0

            if self.neighbor_col and self.neighbor_vals:
                recv_mask = self._make_binary_mask(
                    self.neighbor_vals, self.neighbor_col, adata
                )
                mask_expr   &= recv_mask
                mask_noexpr &= recv_mask

        g1 = np.where(mask_expr)[0]
        g2 = np.where(mask_noexpr)[0]
        if len(g1) < 50 or len(g2) < 50:
            raise ValueError("Not enough cells in one of the groups")
        return g1, g2

    # … rest of your methods (create_adata, de_bootstrap, de_test) remain the same,
    #    but wherever you referenced celltype_sender['celltype'] or
    #    celltype_receiver['rctd_celltype'], switch to using your
    #    self.celltype_col, self.celltype_vals, self.neighbor_col, self.neighbor_vals …
