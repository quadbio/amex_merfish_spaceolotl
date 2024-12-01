import numpy as np
import scanpy as sc

from _logging import _setup_logger
logger = _setup_logger()

def corr_ann(
    adata: sc.AnnData,
    basis: str,
    obs_keys: list[str] = ['nCount_RNA', 'nFeature_RNA'],
    components: list[int] = [1, 2]
    ) -> None:
    
    """Compute and log correlations between specified observation keys and embedding components.

    Args:
        adata (sc.AnnData): The AnnData object containing the data.
        basis (str): The basis to check in adata.obsm e.g. X_pca, X_umap
        obs_keys (list[str], optional): List of observation keys to correlate. Defaults to ['nCount_RNA', 'nFeature_RNA'].
        components (list[int], optional): List of components to correlate with. Defaults to [1, 2].

    Raises:
        ValueError: If the basis or any obs_key does not exist.
    """

    if basis not in adata.obsm.keys():
        raise ValueError(f'The basis "{basis}" has not been computed yet.')

    for key in obs_keys:
        if key not in adata.obs.keys():
            raise ValueError(f'The key "{key}" does not exist in adata.obs.')

    X_em = adata.obsm[basis]
    max_dim = X_em.shape[1]

    for key in obs_keys:
        for comp in components:
            if comp > max_dim:
                logger.warning(f"Component {comp} exceeds available components ({max_dim}) in basis '{basis}'. Skipping this component.")
                continue
            
            coordinates = X_em[:, comp - 1]
            ann = adata.obs[key]
            corr = np.corrcoef(coordinates, ann)[0, 1]
            logger.info(f'Correlation between "{key}" and component {comp} of basis "{basis}" is {corr:.2f}.')

def std_sc_wf(
    adata: sc.AnnData,
    min_counts: int = 10,
    min_cells: int = 0,
    normalize_to: float = 1e4,
    scale: bool = True,
    regress_out: str | None = None,
    scale_clip: float = 10,
    n_pcs: int = 50,
    n_umap: int = 2,
    n_neighbors: int = 15,
    resolutions: list[float] | None = [1.0], 
    de = False,
    random_state: int = 0,
    copy: bool = True,
    **kwargs
    ) -> sc.AnnData | None:
    
    """ 
    Going through all standard motions of sc analysis

    Args:
        adata (sc.AnnData)                : The AnnData object containing the data.
        min_counts (int, optional)        : Minimum counts per cell to retain. Defaults to 10.
        min_cells (int, optional)         : Minimum cells per gene to retain. Defaults to 0.
        normalize_to (float, optional)    : Target sum for normalization. Defaults to 1e4.
        scale (bool, optional)            : Whether to scale the data. Defaults to True.
        regress_out (str | None, optional): Variable(s) to regress out. Defaults to None.
        scale_clip (float, optional)      : Maximum value for scaling. Defaults to 10.
        n_pcs (int, optional)             : Number of principal components to compute. Defaults to 50.
        n_umap (int, optional)            : Number of UMAP components to compute. Defaults to 2.
        n_neighbors (int, optional)       : Number of neighbors for UMAP. Defaults to 15.
        resoltions (list | None, optional): List of clustering resolutions. Defaults to 1.0
        random_state (int, optional)      : Random state for reproducibility. Defaults to 0.
        copy (bool, optional)             : If True, return a copy of the AnnData object. Defaults to True.
        **kwargs (optional)               : passed to corr_ann.

    Returns:
        sc.AnnData | None: The processed AnnData object or None if copy is False.
    """

    if copy:
        adata = adata.copy()

    # Filtering
    cells_before = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=min_counts)
    cells_after = adata.n_obs

    genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    genes_after= adata.n_vars

    logger.info(f"Removed {cells_before - cells_after} cells and {genes_before - genes_after} genes. Object now has shape {adata.n_obs} cells x {adata.n_vars} genes.")

    # Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=normalize_to)
    sc.pp.log1p(adata)

    adata.raw = adata.copy()

    if regress_out is not None:
        sc.pp.regress_out(adata, regress_out)

    if scale:
        sc.pp.scale(adata, max_value=scale_clip)

    sc.pp.pca(adata, svd_solver='arpack', random_state = random_state)
    logger.info(f"PCA completed with {n_pcs} components.")

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)
    logger.info(f"Neighbors computed with {n_neighbors} neighbors.")

    if resolutions:
        for res in resolutions:
            sc.tl.leiden(adata,
                         resolution = res,
                         random_state = random_state,
                         flavor="igraph",
                         n_iterations=2,
                         key_added = "leiden_{}".format(res))
            logger.info(f"Clustering at resolution {str(res)}.")

            if de:
                key = "leiden_{}".format(res)
                sc.tl.dendrogram(adata, groupby = key, key_added = key)
                sc.tl.rank_genes_groups(adata, groupby = key, key_added = key)
                logger.info(f"Detecting DE genes for {key}.")


    sc.tl.umap(adata, random_state = random_state, n_components = n_umap)
    logger.info(f"UMAP computed with {n_umap} dimensions.")

    if kwargs:
        corr_ann(adata, **kwargs)

    return adata if copy else None
